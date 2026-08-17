## ---------------------------------------------------------------------------
## preprocess.R -- Stage 1.
##
## THE PREP CHAIN
## --------------
## sesame's default prep string is "QCDPB":
##     Q  qualityMask            flag design-problem probes
##     C  inferInfiniumIChannel  correct the declared colour channel
##     D  dyeBiasNL              non-linear dye bias correction
##     P  pOOBAH                 detection p-values -> mask
##     B  noob                   background correction
##
## methylQC v3 runs, per sample:
##
##     read IDAT
##       -> stats_channel()   -> Infinium-I channel switch counts (BEFORE "C")
##       -> prepSesame "C"    -> measure intensity / dye bias / probe counts
##       -> prepSesame "DB"   -> dye bias correction, then noob
##       -> ELBAR(return.pval = TRUE)                           [ONE call]
##       -> getBetas(mask = FALSE)
##
## Four deliberate differences from v2, each traced to sesame's source:
##
## 1. ELBAR replaces pOOBAH. ELBAR derives its background empirically from the
##    data's own intensity profile rather than assuming out-of-band signal is
##    background.
##
## 2. Detection runs ONCE, after noob, and its p-values are kept. v2 ran
##    pOOBAH inside openSesame (which set the mask) and then ran it AGAIN in
##    extract_detP() on the finished object, so the stored p-values did not
##    match the stored mask.
##
##    The ordering deserves care, because the right answer differs by method.
##    pOOBAH DEFINES background as the out-of-band signal, and noob rewrites
##    MG/UG/MR/UR (sesame R/background.R:112-118) -- exactly those columns --
##    so pOOBAH after noob compares corrected foreground against corrected
##    background and is incoherent. ELBAR does not: it re-derives background
##    empirically from whatever data it is handed, so running it after noob is
##    self-consistent, and running it BEFORE noob breaks it. methylQC 3.0.1
##    carried the pOOBAH-specific argument across to ELBAR and put detection
##    before noob; the consequence is documented at the call site below.
##
## 3. "Q" is dropped. qualityMask() writes only to sdf$mask
##    (sesame R/mask.R:176-186 -> addMask at :15-21) and nothing downstream in
##    this chain reads that column: dyeBiasNL defaults to mask = TRUE, which
##    means "use all Infinium-I probes INCLUDING masked ones"
##    (R/dye_bias.R:123-124); ELBAR calls signalMU(mask = FALSE); noob filters
##    by nonuniqMask(platform) directly (R/background.R:86); and
##    getBetas(mask = FALSE) skips masking entirely (R/sesame.R:209).
##    "CDB" and "QCDB" therefore produce identical numbers. The design mask is
##    stored separately instead, as a vector, so it stays independently
##    addressable -- which matters when a clock probe sits in the design mask
##    and you want its value anyway.
##
## 4. "I" is not in the string because we call ELBAR() ourselves after B.
##    Running "I" as well would compute the same thing a second time and throw
##    away the p-values.
##
## WHERE THE STATISTICS ARE MEASURED
## ---------------------------------
## Each group is taken at the last point before the step that would destroy it,
## verified on EPIC.1.SigDF rather than assumed:
##
##   channel switches  BEFORE "C".  "C" resolves the very disagreement the
##                     statistic counts, so afterwards R2G and G2R are 0 by
##                     construction and R2R/G2G merely repeat num_probes_IR and
##                     num_probes_IG. Before "C" the same array reports 65 and
##                     904 switches. v3.0.1 measured these after "C" and so
##                     recorded four uninformative columns for every sample.
##
##   dye bias          AFTER "C", BEFORE "D".  "D" corrects dye bias by
##                     definition: RGratio goes 1.512 -> 0.999 across it.
##
##   out-of-band       AFTER "C", BEFORE "B".  These summarise background, and
##                     noob subtracts background: mean_oob_red goes 545 -> 278.
##                     Measuring after noob measures what background remains
##                     after background was removed.
##
##   mean intensity    AFTER "C", BEFORE "D"/"B", though this one matters less
##                     than the others. Scaling an array to 75/50/25/10% of its
##                     signal is recovered as 0.75/0.50/0.25/0.10 before noob
##                     and 0.751/0.503/0.254/0.105 after, so the cohort-relative
##                     MAD rule would work either way. What noob does change is
##                     the absolute level (about 10% lower), and `intfloor` is
##                     an absolute threshold, so the uncorrected scale is the
##                     one it can be calibrated against.
## ---------------------------------------------------------------------------

## Thresholds at which per-probe failure counts are cached, so that qcplots()
## can retune `detp` without re-reading the detection matrix.
.MQC_PGRID <- c(0.01, 0.02, 0.05, 0.10, 0.20)

## Bin count for the cached per-sample p-value histogram. 1000 bins puts the
## call rate at any threshold within 0.001 of the exact value.
.MQC_PBINS <- 1000L


#' Process a single IDAT pair
#'
#' Never throws: any failure is caught and returned in the \code{error} field
#' with fail-closed placeholders, so one corrupt array cannot abort a cohort.
#'
#' @param pfx IDAT basename
#' @param platform platform string
#' @param addr annotation from \code{\link{mqcanno}}
#' @param keep_sdf retain the preprocessed SigDF
#' @param want_design also return the design mask (only needed for one sample)
#' @return a list with \code{betas}, \code{pvals}, \code{stats}, \code{phist},
#'   optionally \code{sdf} and \code{design}, and \code{error}.
#' @keywords internal
#' @noRd
process_one <- function(pfx, platform, addr, keep_sdf = FALSE,
                        want_design = FALSE) {
  res <- list(betas = NULL, pvals = NULL, stats = NULL, phist = NULL,
              sdf = NULL, design = NULL, error = NA_character_,
              warnings = character(0))

  ## ELBAR warns when its assumptions fail (dichotomous background, no
  ## variation), and those warnings are the reason to distrust a sample's
  ## detection p-values. They are raised inside a forked worker, where R's
  ## warning machinery does not reach the parent, so collect them here and
  ## carry them back for the parent to log.
  withCallingHandlers(tryCatch({
    sdf <- read_idat(pfx, platform, addr)

    ## ---- channel switches: measured BEFORE C -----------------------------
    ## inferInfiniumIChannel(summary = TRUE) reports how many Infinium-I probes
    ## disagree with their declared channel. C is what resolves that
    ## disagreement, so after C the off-diagonal counts are zero by
    ## construction and the diagonal ones are just the type I probe counts
    ## again. Measured before C they say how many probes the array mis-declared,
    ## which is a real chemistry diagnostic.
    ch <- stats_channel(sdf)

    ## ---- stage C, then the remaining instrument diagnostics --------------
    ## Sex-chromosome intensity is measured HERE, in the worker, rather than in
    ## Stage 2 from a retained SigDF. Two numbers per sample cost nothing to
    ## carry, and it keeps the sex check working when savesdf is downgraded by
    ## the memory guard -- which happens on exactly the large cohorts where a
    ## sample swap is most likely.
    sdf <- sesame::prepSesame(sdf, "C")
    st <- stats_raw(sdf)
    for (k in names(ch)) st[[k]] <- ch[[k]]
    si <- sexintensity(sdf)
    st$sex_chrX_intensity <- si$chrX
    st$sex_chrY_intensity <- si$chrY
    res$stats <- st
    if (isTRUE(want_design)) res$design <- designmask(sdf)

    ## ---- stages D and B: dye bias, then noob -----------------------------
    sdf <- sesame::prepSesame(sdf, "DB")

    ## ---- detection: ONE ELBAR call, after noob ---------------------------
    ## ELBAR must run AFTER noob. It locates background as the low-intensity
    ## population whose beta values are undifferentiated around 0.5, and
    ## without noob the two channels still carry different additive offsets, so
    ## dark type II probes (M green, U red) are pushed towards 0 and 1 while
    ## dark type I probes (both alleles one channel) stay near 0.5. The pooled
    ## distribution is then bimodal, ELBAR reports "Background signal is
    ## dichotomous", abandons the search and defines background from the ten
    ## dimmest probes alone. Measured on EPIC.1.SigDF: beta spread across the
    ## 500 dimmest probes is 0.758 before noob against 0.083 after, and the
    ## fallback leaves eight distinct p-values for the whole array with every
    ## probe at p = 0 -- detection silently switched off.
    pv <- sesame::ELBAR(sdf, return.pval = TRUE)
    pv[is.na(pv)] <- 1.0                    # fail closed
    res$pvals <- pv
    res$phist <- phist_one(pv)

    ## ---- unmasked betas --------------------------------------------------
    b <- sesame::getBetas(sdf, mask = FALSE)
    res$betas <- b

    if (isTRUE(keep_sdf)) res$sdf <- sdf
  }, error = function(e) {
    res$error <<- conditionMessage(e)
  }), warning = function(w) {
    res$warnings <<- c(res$warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

  res
}

#' Instrument and background diagnostics, measured before noob
#'
#' @param sdf a SigDF that has had "C" applied but not "D" or "B"
#' @return a one-row data.frame
#' @keywords internal
#' @noRd
stats_raw <- function(sdf) {
  qc <- sesame::sesameQC_calcStats(
    sdf, c("intensity", "numProbes", "dyeBias"))
  s <- sesame::sesameQC_getStats(qc, drop = FALSE)

  num <- function(k) {
    v <- s[[k]]
    if (is.null(v) || !length(v)) NA_real_ else as.numeric(v[1])
  }
  data.frame(
    mean_intensity_raw = num("mean_intensity"),
    mean_intensity_mu  = num("mean_intensity_MU"),
    mean_inb_grn       = num("mean_inb_grn"),
    mean_inb_red       = num("mean_inb_red"),
    mean_oob_grn       = num("mean_oob_grn"),
    mean_oob_red       = num("mean_oob_red"),
    rg_distort         = num("RGdistort"),
    rg_ratio           = num("RGratio"),
    n_probes           = num("num_probes"),
    n_probes_ii        = num("num_probes_II"),
    n_probes_ir        = num("num_probes_IR"),
    n_probes_ig        = num("num_probes_IG"),
    na_intensity_m     = num("na_intensity_M"),
    na_intensity_u     = num("na_intensity_U"),
    stringsAsFactors = FALSE)
}

#' Infinium-I channel switch counts, measured before channel inference
#'
#' \code{inferInfiniumIChannel(summary = TRUE)} reports how many Infinium-I
#' probes were read in a channel other than the one the manifest declares. That
#' disagreement is exactly what prep code \code{"C"} resolves, so the statistic
#' has to be taken beforehand. Measured after \code{"C"} on a real EPIC array,
#' \code{R2G} and \code{G2R} are 0 by construction and \code{R2R}/\code{G2G}
#' reproduce \code{num_probes_IR}/\code{num_probes_IG} exactly -- four columns
#' carrying no information. Measured before \code{"C"} the same array reports
#' 65 red-to-green and 904 green-to-red switches, which is a genuine chemistry
#' diagnostic.
#'
#' @param sdf a \code{SigDF} that has NOT had \code{"C"} applied.
#' @return a one-row data.frame of the four switch counts.
#' @keywords internal
#' @noRd
stats_channel <- function(sdf) {
  s <- tryCatch({
    qc <- sesame::sesameQC_calcStats(sdf, "channel")
    sesame::sesameQC_getStats(qc, drop = FALSE)
  }, error = function(e) list())
  num <- function(k) {
    v <- s[[k]]
    if (is.null(v) || !length(v)) NA_real_ else as.numeric(v[1])
  }
  data.frame(
    inf1_r2r = num("InfI_switch_R2R"), inf1_g2g = num("InfI_switch_G2G"),
    inf1_r2g = num("InfI_switch_R2G"), inf1_g2r = num("InfI_switch_G2R"),
    stringsAsFactors = FALSE)
}

#' Histogram of one sample's detection p-values
#' @param pv numeric p-values
#' @return integer vector of length \code{.MQC_PBINS}
#' @keywords internal
#' @noRd
phist_one <- function(pv) {
  ## Bins are RIGHT-closed -- bin k covers ((k-1)/nb, k/nb] -- so that summing
  ## bins 1..k answers "how many p <= k/nb", exactly the comparison callrate()
  ## makes on the matrix directly. With left-closed bins the two paths
  ## disagreed for any p sitting exactly on the threshold.
  pv[is.na(pv)] <- 1                      # fail closed, as process_one does
  b <- findInterval(pv, seq(0, 1, length.out = .MQC_PBINS + 1L), left.open = TRUE)
  b[b < 1L] <- 1L                         # p == 0 belongs in the first bin
  b[b > .MQC_PBINS] <- .MQC_PBINS
  tabulate(b, nbins = .MQC_PBINS)
}


## ---------------------------------------------------------------------------
## Stage 1 driver
## ---------------------------------------------------------------------------

#' Preprocess a cohort of IDAT pairs
#'
#' Runs the three-stage prep chain over every sample, assembling the beta and
#' detection p-value matrices plus the design mask and per-sample diagnostics.
#'
#' @param ss sample sheet from \code{\link{discover}}; needs \code{sample_id}
#'   and \code{Basename}.
#' @param platform platform string.
#' @param outdir output directory (used for per-batch files and checkpoints).
#' @param workers number of forked workers; \code{NULL} inherits from the plan.
#' @param batch samples per task; \code{NULL} inherits from the plan.
#' @param extreme extreme memory conservation: write per-batch files and never
#'   hold a full matrix.
#' @param savesdf retain preprocessed SigDF objects.
#' @param memcap abort if resident memory exceeds this many bytes;
#'   \code{NULL} disables the guard.
#' @param logger optional logger.
#' @return a list with \code{betas}, \code{detp}, \code{design}, \code{stats},
#'   \code{failed}, \code{phist}, \code{pfail} and \code{probe_ids}. Under
#'   \code{extreme = TRUE} the matrices are \code{NULL} and batch files are on
#'   disk instead.
#' @export
#' @examples
#' \dontrun{
#' ss <- discover("~/idats")
#' out <- runsesame(ss, "EPIC", "~/qcout")
#' }
runsesame <- function(ss, platform, outdir,
                      workers = NULL, batch = NULL, extreme = NULL,
                      savesdf = NULL, memcap = NULL, logger = NULL) {

  logger <- logger %||% nulllog()
  cfg <- mqcopts()
  extreme <- .opt(extreme, "extreme", cfg)
  savesdf <- .opt(savesdf, "savesdf", cfg)

  stopifnot(all(c("sample_id", "Basename") %in% names(ss)))
  n <- nrow(ss)
  if (n < 1L) stop("no samples to process", call. = FALSE)
  if (anyDuplicated(ss$sample_id))
    stop("sample_id values must be unique", call. = FALSE)

  ## Resolve annotation ONCE, here in the parent, and pre-warm the cache.
  ## Everything below passes `addr` explicitly so no worker ever calls
  ## sesameDataGet(). See annotation.R for why.
  mqcprewarm(logger)
  addr <- mqcanno(platform, logger)

  workers <- .resolve_workers(workers, cfg, n, logger)
  batch <- .resolve_batch(batch, cfg, n, workers)
  logger$log("stage1", sprintf("%d samples, %d worker(s), batch size %d%s",
                               n, workers, batch,
                               if (isTRUE(extreme)) ", extreme memory mode" else ""))

  chunks <- split(seq_len(n), ceiling(seq_len(n) / batch))

  ## ---- sample 1 alone, to establish the probe set and the design mask ----
  first <- process_one(ss$Basename[1], platform, addr,
                       keep_sdf = savesdf, want_design = TRUE)
  if (!is.na(first$error))
    stop("the first sample failed to process, so the probe set cannot be ",
         "established: ", first$error, call. = FALSE)

  probe_ids <- names(first$betas)
  np <- length(probe_ids)
  design <- align_to(first$design, probe_ids, "design mask")
  logger$log("stage1", sprintf("%d probes; %d masked by design (%.1f%%)",
                               np, sum(design), 100 * mean(design)))

  ## ---- allocate ----------------------------------------------------------
  keep_full <- !isTRUE(extreme)
  betas <- if (keep_full) matrix(NA_real_, np, n,
                                 dimnames = list(probe_ids, ss$sample_id)) else NULL
  detp  <- if (keep_full) matrix(NA_real_, np, n,
                                 dimnames = list(probe_ids, ss$sample_id)) else NULL
  sdfs  <- if (isTRUE(savesdf)) vector("list", n) else NULL
  if (isTRUE(savesdf)) names(sdfs) <- ss$sample_id

  phist <- matrix(0L, .MQC_PBINS, n, dimnames = list(NULL, ss$sample_id))
  pfail <- matrix(0L, np, length(.MQC_PGRID),
                  dimnames = list(probe_ids, sprintf("p%g", .MQC_PGRID)))
  stats_l <- vector("list", n)
  failed <- data.frame(sample_id = character(0), basename = character(0),
                       error = character(0), stringsAsFactors = FALSE)

  if (isTRUE(extreme)) mqcmakedirs(outdir, batches = TRUE)

  mem_track <- data.frame(done = integer(0), rss = numeric(0))
  t0 <- Sys.time()

  ## ---- absorb one worker result into the accumulators --------------------
  absorb <- function(r, i, bcol) {
    sid <- ss$sample_id[i]
    for (w in unique(r$warnings %||% character(0)))
      logger$log("stage1", sprintf("sample %s: %s", sid, w))
    if (!is.na(r$error)) {
      failed <<- rbind(failed, data.frame(
        sample_id = sid, basename = ss$Basename[i], error = r$error,
        stringsAsFactors = FALSE))
      logger$log("stage1", sprintf("sample %s FAILED: %s", sid, r$error),
                 warn = TRUE)
      ## fail closed: absent data must look maximally suspect, never perfect
      stats_l[[i]] <<- .empty_stats()
      phist[, i] <<- c(rep(0L, .MQC_PBINS - 1L), np)
      pfail <<- pfail + 1L
      return(invisible(NULL))
    }
    b <- align_to(r$betas, probe_ids, paste0("betas for ", sid))
    p <- align_to(r$pvals, probe_ids, paste0("detection p-values for ", sid))

    if (keep_full) {
      betas[, i] <<- b
      detp[, i]  <<- p
    } else {
      ## bcol is an environment, so ordinary assignment mutates it in place;
      ## `<<-` here would skip the local frame and look for `bcol` outside.
      bcol$b[[length(bcol$b) + 1L]] <- b
      bcol$p[[length(bcol$p) + 1L]] <- p
      bcol$s <- c(bcol$s, sid)
    }
    phist[, i] <<- r$phist
    for (g in seq_along(.MQC_PGRID))
      pfail[, g] <<- pfail[, g] + as.integer(p > .MQC_PGRID[g])
    stats_l[[i]] <<- r$stats
    if (isTRUE(savesdf)) sdfs[[i]] <<- r$sdf
    invisible(NULL)
  }

  ## ---- iterate over chunks ----------------------------------------------
  bpp <- .make_bpparam(workers)
  done <- 0L
  aborted <- FALSE
  warned_mem <- FALSE

  for (ci in seq_along(chunks)) {
    idx <- chunks[[ci]]
    bcol <- new.env(parent = emptyenv())
    bcol$b <- list(); bcol$p <- list(); bcol$s <- character(0)

    idx_run <- idx
    if (ci == 1L) {
      absorb(first, idx[1], bcol)
      idx_run <- idx[-1]
      first <- NULL          # release the retained SigDF before the next chunk
    }

    if (length(idx_run)) {
      res <- .run_chunk(ss$Basename[idx_run], platform, addr, savesdf, bpp, logger)
      for (k in seq_along(idx_run)) absorb(res[[k]], idx_run[k], bcol)
      rm(res)
    }

    if (isTRUE(extreme) && length(bcol$s)) .write_batch(bcol, probe_ids, ci, outdir)

    done <- done + length(idx)
    rss <- procmem()
    mem_track <- rbind(mem_track, data.frame(done = done, rss = rss))
    gc(verbose = FALSE)

    logger$log("stage1", sprintf("batch %d/%d complete (%d/%d samples, RSS %s)",
                                 ci, length(chunks), done, n, fmtbytes(rss)))

    ## memcap is a PARENT-process budget: qcplan() sizes the total envelope and
    ## subtracts the worker allowance before handing it over, because procmem()
    ## can only see this process. Comparing parent RSS against the whole-run
    ## budget would let the fleet's memory grow unseen.
    ## Crossing the cap is a trigger to LOOK, not a reason to stop. Resident
    ## memory includes R's high-water mark and pages freed but not returned to
    ## the OS, so it drifts up even when live data is flat, and the cap is
    ## necessarily an estimate. Abort only when the growth actually extrapolates
    ## past the cap for the full cohort. 3.0.1 aborted on the instantaneous
    ## reading and killed runs whose own projection showed a zero shortfall.
    if (!is.null(memcap) && is.finite(memcap) && !is.na(rss) && rss > memcap) {
      proj <- .project_mem(mem_track, n)
      if (is.finite(proj) && proj <= memcap) {
        if (!warned_mem) {
          logger$log("memory", sprintf(
            paste("resident memory %s is over the cap of %s, but growth",
                  "projects to %s for all %d samples, which still fits;",
                  "continuing"),
            fmtbytes(rss), fmtbytes(memcap), fmtbytes(proj), n), warn = TRUE)
          warned_mem <- TRUE
        }
      } else {
        .memory_abort(mem_track, n, memcap, outdir, betas, detp, design,
                      probe_ids, ss, stats_l, failed, logger)
        aborted <- TRUE
        break
      }
    }
  }
  .stop_bpparam(bpp)

  if (aborted)
    stop("run aborted at the memory cap after ", done, " of ", n,
         " samples. Partial results were written; see the log for the ",
         "projected requirement and a suggested batch size.", call. = FALSE)

  stats_df <- .bind_stats(stats_l, ss$sample_id)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  logger$log("stage1", sprintf("Stage 1 complete in %s; %d failed sample(s)",
                               fmtdur(elapsed), nrow(failed)))

  list(betas = betas, detp = detp, design = design, stats = stats_df,
       failed = failed, phist = phist, pfail = pfail, pgrid = .MQC_PGRID,
       probe_ids = probe_ids, sdfs = sdfs, elapsed = elapsed,
       platform = platform, n_batches = length(chunks))
}


## ---------------------------------------------------------------------------
## Internals
## ---------------------------------------------------------------------------

## stop.on.error = FALSE does NOT make bplapply fault-tolerant: it controls
## when the error is raised, not whether. A worker that dies outright (OOM
## killer, fork failure) still propagates, and because runsesame() accumulates
## into matrices held in the parent, that discarded every completed batch.
## process_one() already catches its own errors, so any throw reaching here is
## infrastructural and a serial redo of that one chunk is the cheap repair.
.run_chunk <- function(prefixes, platform, addr, savesdf, bpp, logger = NULL) {
  logger <- logger %||% nulllog()
  fn <- function(pfx) process_one(pfx, platform, addr, keep_sdf = savesdf)
  if (is.null(bpp)) return(lapply(prefixes, fn))

  res <- tryCatch(BiocParallel::bplapply(prefixes, fn, BPPARAM = bpp),
                  error = function(e) {
                    logger$log("stage1", paste(
                      "parallel dispatch failed, re-running this batch serially:",
                      conditionMessage(e)), warn = TRUE)
                    NULL
                  })
  if (is.null(res) || length(res) != length(prefixes))
    return(lapply(prefixes, fn))

  ## Repair anything that came back in the wrong shape rather than trusting it.
  bad <- vapply(res, function(r) !is.list(r) || is.null(r[["error"]]), logical(1))
  for (i in which(bad))
    res[[i]] <- tryCatch(fn(prefixes[i]),
                         error = function(e) .failed_result(conditionMessage(e)))
  res
}

## The shape process_one() promises, for the paths that cannot call it.
.failed_result <- function(msg) {
  list(betas = NULL, pvals = NULL, stats = NULL, phist = NULL,
       sdf = NULL, design = NULL, error = as.character(msg))
}

.make_bpparam <- function(workers) {
  if (workers <= 1L || !can_fork()) return(NULL)
  ## stop.on.error = FALSE is belt and braces: process_one() already catches
  ## everything itself, so a worker should never throw.
  BiocParallel::MulticoreParam(workers = workers, stop.on.error = FALSE,
                               progressbar = FALSE)
}

.stop_bpparam <- function(bpp) {
  if (is.null(bpp)) return(invisible(NULL))
  try(BiocParallel::bpstop(bpp), silent = TRUE)
  invisible(NULL)
}

.resolve_workers <- function(workers, cfg, n, logger) {
  w <- .opt(workers, "workers", cfg)
  if (is.null(w)) {
    pl <- tryCatch(qcplan_quiet(n), error = function(e) NULL)
    w <- pl$workers %||% 1L
  }
  w <- max(1L, min(as.integer(w), n))
  if (!can_fork() && w > 1L) {
    logger$log("stage1",
               "forked parallelism is unavailable on this platform; running serially",
               warn = TRUE)
    w <- 1L
  }
  w
}

.resolve_batch <- function(batch, cfg, n, workers) {
  b <- .opt(batch, "batch", cfg)
  if (is.null(b)) {
    ## Keep at least 4 tasks per worker so one slow sample cannot strand the
    ## rest at the tail. Batch size is a memory dial; per-task overhead is
    ## milliseconds against seconds of per-sample work, so time is flat.
    b <- max(1L, min(50L, floor(n / max(1L, 4L * workers))))
    if (b < 1L) b <- 1L
  }
  max(1L, min(as.integer(b), n))
}

.empty_stats <- function() {
  z <- stats_raw_template()
  z[] <- NA_real_
  z
}

stats_raw_template <- function() {
  data.frame(
    mean_intensity_raw = NA_real_, mean_intensity_mu = NA_real_,
    mean_inb_grn = NA_real_, mean_inb_red = NA_real_,
    mean_oob_grn = NA_real_, mean_oob_red = NA_real_,
    rg_distort = NA_real_, rg_ratio = NA_real_,
    n_probes = NA_real_, n_probes_ii = NA_real_,
    n_probes_ir = NA_real_, n_probes_ig = NA_real_,
    inf1_r2r = NA_real_, inf1_g2g = NA_real_,
    inf1_r2g = NA_real_, inf1_g2r = NA_real_,
    na_intensity_m = NA_real_, na_intensity_u = NA_real_,
    sex_chrX_intensity = NA_real_, sex_chrY_intensity = NA_real_,
    stringsAsFactors = FALSE)
}

.bind_stats <- function(lst, sample_ids) {
  tmpl <- stats_raw_template()
  lst <- lapply(lst, function(x) {
    if (is.null(x)) return(tmpl)
    miss <- setdiff(names(tmpl), names(x))
    for (m in miss) x[[m]] <- NA_real_
    x[, names(tmpl), drop = FALSE]
  })
  out <- do.call(rbind, lst)
  out <- cbind(sample_id = sample_ids, out, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

.write_batch <- function(bcol, probe_ids, ci, outdir) {
  bm <- do.call(cbind, bcol$b); colnames(bm) <- bcol$s; rownames(bm) <- probe_ids
  pm <- do.call(cbind, bcol$p); colnames(pm) <- bcol$s; rownames(pm) <- probe_ids
  d <- mqcpath(outdir, "batches", create = TRUE)
  save_rds_atomic(bm, file.path(d, sprintf("betas_%04d.rds", ci)))
  save_rds_atomic(pm, file.path(d, sprintf("detP_%04d.rds", ci)))
  invisible(NULL)
}

## Extrapolate resident memory to the full cohort. Memory grows linearly in
## samples because the accumulating matrices and SigDF list do, so a straight
## line through the completed batches gives bytes per sample as the slope and
## fixed overhead as the intercept.
.project_mem <- function(mem_track, n) {
  mt <- mem_track[!is.na(mem_track$rss), , drop = FALSE]
  if (nrow(mt) >= 3L) {
    fit <- stats::lm(rss ~ done, data = mt)
    cf <- stats::coef(fit)
    return(unname(cf[1] + cf[2] * n))
  }
  if (nrow(mt) >= 1L)
    return(mt$rss[nrow(mt)] / mt$done[nrow(mt)] * n)
  NA_real_
}

## Warn, checkpoint, extrapolate. A relative memory model is fitted across
## completed batches: resident memory grows linearly in samples because the
## matrices do, so slope is bytes per sample and intercept is fixed overhead.
.memory_abort <- function(mem_track, n, memcap, outdir, betas, detp, design,
                          probe_ids, ss, stats_l, failed, logger) {
  logger$log("memory", sprintf("resident memory %s exceeded the cap of %s",
                               fmtbytes(mem_track$rss[nrow(mem_track)]),
                               fmtbytes(memcap)), warn = TRUE)

  ## Checkpoint first -- aborting without saving throws away hours -- but NOT
  ## onto the canonical paths. 3.0.0 wrote the partially-filled matrix to
  ## betas_all.rds, where mqccheckdir() and prep() then accepted it as a
  ## complete Stage 1 and silently analysed a cohort whose unprocessed tail was
  ## all NA. Partial output goes to its own names and prep() ignores it.
  ok <- tryCatch({
    mqcmakedirs(outdir)
    d <- mqcpath(outdir, "matrices", create = TRUE)
    if (!is.null(betas)) save_rds_atomic(betas, file.path(d, "betas_partial.rds"))
    if (!is.null(detp))  save_rds_atomic(detp,  file.path(d, "detP_partial.rds"))
    save_rds_atomic(design, file.path(d, "design_mask_partial.rds"))
    TRUE
  }, error = function(e) { logger$log("memory",
      paste("checkpoint failed:", conditionMessage(e)), warn = TRUE); FALSE })
  if (ok)
    logger$log("memory", paste(
      "partial matrices checkpointed as *_partial.rds. They are NOT a complete",
      "Stage 1 and prep() will not read them; re-run with a smaller batch or",
      "extreme = TRUE."))

  mt <- mem_track[!is.na(mem_track$rss), , drop = FALSE]
  if (nrow(mt) >= 3L) {
    fit <- stats::lm(rss ~ done, data = mt)
    b0 <- unname(stats::coef(fit)[1]); b1 <- unname(stats::coef(fit)[2])
  } else if (nrow(mt) >= 1L) {
    b1 <- mt$rss[nrow(mt)] / mt$done[nrow(mt)]; b0 <- 0
  } else {
    b1 <- NA_real_; b0 <- NA_real_
  }
  total <- b0 + b1 * n
  avail <- sysmem()$available
  suggest <- if (is.na(b1) || b1 <= 0) NA_real_ else
    floor((memcap - b0) / b1)

  logger$log("memory", sprintf(
    "measured %s per sample plus %s fixed; projected %s for all %d samples",
    fmtbytes(b1), fmtbytes(b0), fmtbytes(total), n))
  if (!is.na(avail))
    logger$log("memory", sprintf("available memory is %s; shortfall %s",
                                 fmtbytes(avail), fmtbytes(max(0, total - avail))))
  if (!is.na(suggest) && suggest >= 1)
    logger$log("memory", sprintf(
      "process in runs of about %d samples, or set extreme = TRUE",
      max(1L, as.integer(suggest))))
  invisible(NULL)
}
