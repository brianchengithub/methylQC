## ---------------------------------------------------------------------------
## pipeline.R -- the user-facing entry points.
##
## Goal:    take a directory of IDATs to a finished QC report in one call.
## Approach: two stages behind one front door, with every stage also callable
##          on its own for re-runs.
## Inputs:  a directory of IDAT pairs, optionally a sample sheet.
## Outputs: an output directory (see mqcpath()/utils.R for the layout).
## Usage:   pipeline("~/idats", "~/qcout")
##
##   pipeline() one call: preflight, Stage 1, Stage 2, summary
##   qcplan()   preflight only: memory, workers, batch size, SigDF retention
##   qc()       Stage 1 (+ Stage 2 unless stage2 = FALSE)
##   prep()     Stage 2 alone, against an existing output directory
##   qcplots()  regenerate the threshold-dependent panels from cache
##   makemat()  reassemble per-batch matrices written under extreme mode
##
## THE SAMPLE SHEET
## ----------------
## There is exactly ONE sample_sheet.csv. Stage 1 creates it with discovery
## information and per-sample diagnostics; Stage 2 adds flag, sex, age and
## cell composition columns to that same file, for every sample. It is never
## subset. Columns are assigned rather than joined, so re-running is
## idempotent and never produces .x/.y duplicates.
##
## RUN FACTS
## ---------
## Stage 1 writes run_info.rds. Without it a standalone prep() had no idea what
## Stage 1 did and rewrote METHODS.txt with NA placeholders, a 0.0% design-mask
## figure and no EPICv2 section.
## ---------------------------------------------------------------------------

#' Run the whole pipeline in one call
#'
#' Preflights the run, processes every IDAT pair, derives quality metrics,
#' calls sex, estimates cell composition and epigenetic age, writes the QC
#' report, and prints a short summary of what was produced.
#'
#' This is the entry point to reach for. \code{\link{qc}}, \code{\link{prep}}
#' and \code{\link{qcplots}} remain available for re-running one stage.
#'
#' @param indir directory containing IDAT files.
#' @param outdir destination directory; created if absent.
#' @param ... passed to \code{\link{qc}}, e.g. \code{workers}, \code{batch},
#'   \code{extreme}, \code{savesdf}, \code{platform}, \code{sheet}.
#' @param quiet suppress the printed summary.
#' @return the output directory, invisibly.
#' @export
#' @examples
#' \dontrun{
#' pipeline("~/idats", "~/qcout")
#' pipeline("~/idats", "~/qcout", workers = 4)
#' }
pipeline <- function(indir, outdir, ..., quiet = FALSE) {
  qc(indir = indir, outdir = outdir, ...)
  if (!isTRUE(quiet)) cat(.run_summary(outdir), sep = "\n")
  invisible(outdir)
}

#' Run the methylQC pipeline
#'
#' Stage 1 reads and preprocesses every IDAT pair; Stage 2 derives quality
#' metrics, calls sex, estimates cell composition and epigenetic age, and
#' writes the QC report.
#'
#' When \code{workers}, \code{batch} and \code{savesdf} are not supplied they
#' are inherited from \code{\link{qcplan}}, which measures actual memory use on
#' this machine.
#'
#' @param indir directory containing IDAT files.
#' @param outdir destination directory; created if absent.
#' @param platform platform string; inferred when \code{NULL}.
#' @param sheet optional explicit sample sheet path.
#' @param workers forked worker processes; \code{NULL} inherits from the plan.
#' @param batch samples per task; \code{NULL} inherits from the plan.
#' @param extreme extreme memory conservation mode.
#' @param savesdf retain preprocessed SigDF objects. The preflight withdraws
#'   this when the cohort is too large to hold them; nothing else depends on
#'   them, so the run is unaffected.
#' @param collapse resolve EPICv2 replicate probes; \code{NA} decides
#'   automatically from the probe identifiers.
#' @param stage2 also run Stage 2.
#' @param plan an existing plan from \code{\link{qcplan}}, to avoid measuring
#'   twice.
#' @return the output directory, invisibly.
#' @export
#' @examples
#' \dontrun{
#' qcplan("~/idats")
#' qc("~/idats", "~/qcout")
#' }
qc <- function(indir, outdir, platform = NULL, sheet = NULL,
               workers = NULL, batch = NULL, extreme = NULL, savesdf = NULL,
               collapse = NULL, stage2 = TRUE, plan = NULL) {

  cfg <- mqcopts()
  extreme <- .opt(extreme, "extreme", cfg)

  mqcmakedirs(outdir, batches = isTRUE(extreme))
  lg <- makelog(mqcpath(outdir, "log", create = TRUE))
  on.exit(lg$close(), add = TRUE)

  lg$log("qc", sprintf("methylQC %s starting",
                       utils::packageVersion("methylQC")))
  checkdeps(stop_on_missing = TRUE)

  ## ---- discovery ---------------------------------------------------------
  ss <- discover(indir, sheet = sheet, logger = lg)
  seen <- unique(stats::na.omit(ss$detected_platform))
  if (length(seen) > 1L)
    stop("this batch mixes array platforms (", paste(seen, collapse = ", "),
         "); process each platform separately.", call. = FALSE)
  platform <- platform %||% seen[1] %||% NA_character_
  if (is.na(platform))
    stop("could not determine the array platform; pass platform= explicitly",
         call. = FALSE)
  lg$log("qc", sprintf("platform: %s", platform))

  cols <- checkmeta(ss, logger = lg)

  ## ---- plan --------------------------------------------------------------
  ## Always run the preflight, even when workers and batch are given, because
  ## it is what decides whether SigDFs can be retained and what the memory
  ## guard's parent budget is.
  if (is.null(plan)) {
    plan <- tryCatch(qcplan(ss = ss, platform = platform, extreme = extreme,
                            savesdf = savesdf, outdir = outdir, verbose = FALSE),
                     error = function(e) {
                       lg$log("qc", paste("preflight failed:", conditionMessage(e)),
                              warn = TRUE)
                       NULL
                     })
    if (!is.null(plan)) for (l in plan$text) lg$log("plan", l)
  }
  workers <- workers %||% plan$workers %||% 1L
  batch   <- batch   %||% plan$batch   %||% 25L
  savesdf <- plan$savesdf %||% .opt(savesdf, "savesdf", cfg)
  parentcap <- (plan$parentcap %|NA|% NULL) %||% .default_memcap(cfg)
  if (isTRUE(plan$savesdf_downgraded))
    for (l in plan$sdf_note) lg$log("memory", sub("^\\s+", "", l), warn = FALSE)
  if (!is.null(plan) && !isTRUE(plan$fits))
    lg$log("qc", paste("the preflight projects that this run will not fit in",
                       "memory; consider extreme = TRUE"), warn = TRUE)

  ## ---- Stage 1 -----------------------------------------------------------
  t0 <- Sys.time()
  s1 <- runsesame(ss, platform, outdir, workers = workers, batch = batch,
                  extreme = extreme, savesdf = savesdf, memcap = parentcap,
                  logger = lg)

  n_pre <- length(s1$probe_ids)
  collapsed <- FALSE
  do_collapse <- .decide_collapse(collapse, cfg, s1$probe_ids, lg)
  if (isTRUE(do_collapse) && isTRUE(extreme)) {
    lg$log("qc", paste("EPICv2 replicate collapsing needs a full matrix and is",
                       "skipped under extreme = TRUE. Run makemat() then",
                       "prep(collapse = TRUE)."), warn = TRUE)
  } else if (isTRUE(do_collapse)) {
    cv <- collapsev2(s1$betas, s1$detp, s1$design,
                     method = cfg$collapsemethod, logger = lg)
    s1$betas <- cv$betas; s1$detp <- cv$detp; s1$design <- cv$design
    s1$probe_ids <- rownames(cv$betas)
    ## Both derived summaries are rebuilt. 3.0.0 rebuilt pfail but not phist,
    ## so prep() (which reads the collapsed matrix) and qcplots() (which reads
    ## the cached pre-collapse histogram) reported different call rates for the
    ## same samples, and qcplots() rewrote the sheet with the second one.
    s1$pfail <- .pfail_from_detp(s1$detp, s1$pgrid)
    s1$phist <- .phist_from_detp(s1$detp)
    collapsed <- TRUE
    attr(s1, "collapsemethod") <- cv$method
  }

  ## ---- persist Stage 1 ---------------------------------------------------
  prov <- list(methylqc_version = as.character(utils::packageVersion("methylQC")),
               sesame_version = as.character(utils::packageVersion("sesame")),
               platform = platform, prep = "C|DB|ELBAR",
               created = as.character(Sys.time()))
  if (!is.null(s1$betas)) {
    attr(s1$betas, "methylqc") <- prov
    save_rds_atomic(s1$betas, mqcpath(outdir, "betas", create = TRUE))
    save_rds_atomic(s1$detp, mqcpath(outdir, "detp"))
    snp <- snpbetas(s1$betas)
    if (!is.null(snp)) save_rds_atomic(snp, mqcpath(outdir, "snpbetas"))
  }
  save_rds_atomic(s1$design, mqcpath(outdir, "designmask", create = TRUE))
  if (isTRUE(savesdf) && !is.null(s1$sdfs))
    save_rds_atomic(s1$sdfs, mqcpath(outdir, "sdfs"))
  write_csv_atomic(s1$failed, mqcpath(outdir, "failedsample", create = TRUE))

  ## ---- Stage 1 sample sheet ---------------------------------------------
  ss <- .assign_cols(ss, s1$stats, "sample_id")
  write_csv_atomic(ss, mqcpath(outdir, "sheet", create = TRUE))
  lg$log("qc", sprintf("Stage 1 written to %s", mqcpath(outdir, "matrices")))

  info <- list(
    pkg_version = prov$methylqc_version, sesame_version = prov$sesame_version,
    date = as.character(Sys.Date()), platform = platform, indir = indir,
    n_samples = nrow(ss), n_batches = length(unique(ss$batch_folder)),
    n_probes = length(s1$probe_ids), n_probes_precollapse = n_pre,
    n_failed = nrow(s1$failed), n_design_masked = sum(s1$design),
    collapsed = collapsed,
    collapsemethod = attr(s1, "collapsemethod") %||% cfg$collapsemethod,
    maskuse = cfg$maskuse, detp = cfg$detp, samplemin = cfg$samplemin,
    intmad = cfg$intmad, intfloor = cfg$intfloor,
    workers = workers, batch = batch, extreme = isTRUE(extreme),
    savesdf = isTRUE(savesdf), savesdf_downgraded = isTRUE(plan$savesdf_downgraded),
    stage1_time = fmtdur(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
    dishref = cfg$dishref, dishmethod = cfg$dishmethod)
  save_rds_atomic(info, mqcpath(outdir, "runinfo", create = TRUE))

  if (!isTRUE(stage2)) {
    writemethods(outdir, info)
    lg$log("qc", "Stage 2 skipped by request")
    return(invisible(outdir))
  }

  if (isTRUE(extreme)) {
    lg$log("qc", paste("extreme mode: Stage 2 needs a full matrix.",
                       "Run makemat() then prep()."), warn = TRUE)
    writemethods(outdir, info)
    return(invisible(outdir))
  }

  prep(outdir, s1 = s1, cols = cols, info = info, logger = lg)
  invisible(outdir)
}


#' Stage 2: metrics, flags, sex, cell composition, age and the QC report
#'
#' @param dir an output directory produced by \code{\link{qc}}.
#' @param s1 in-memory Stage 1 results; read from disk when \code{NULL}.
#' @param cols resolved metadata column names.
#' @param info run facts for METHODS.txt; read from \code{run_info.rds} when
#'   \code{NULL}.
#' @param collapse collapse EPICv2 replicates before Stage 2; \code{NA}
#'   decides automatically. Use this after \code{makemat()} on an
#'   \code{extreme = TRUE} run.
#' @param logger optional logger.
#' @return \code{dir}, invisibly.
#' @export
#' @examples
#' \dontrun{
#' prep("~/qcout")
#' }
prep <- function(dir, s1 = NULL, cols = NULL, info = NULL, collapse = NULL,
                 logger = NULL) {

  own <- is.null(logger)
  lg <- logger %||% makelog(mqcpath(dir, "log", create = TRUE))
  if (own) on.exit(lg$close(), add = TRUE)
  cfg <- mqcopts()

  if (is.null(s1)) {
    mqccheckdir(dir)
    s1 <- .load_stage1(dir, lg)
    ## An extreme-mode run reaches Stage 2 through makemat(), so the matrix on
    ## disk may still carry EPICv2 replicate suffixes.
    if (isTRUE(.decide_collapse(collapse, cfg, s1$probe_ids, lg))) {
      cv <- collapsev2(s1$betas, s1$detp, s1$design,
                       method = cfg$collapsemethod, logger = lg)
      s1$betas <- cv$betas; s1$detp <- cv$detp; s1$design <- cv$design
      s1$probe_ids <- rownames(cv$betas)
      s1$phist <- NULL; s1$pfail <- NULL      # rebuilt by build_cache()
      save_rds_atomic(s1$betas, mqcpath(dir, "betas"))
      save_rds_atomic(s1$detp, mqcpath(dir, "detp"))
      save_rds_atomic(s1$design, mqcpath(dir, "designmask"))
    }
  }
  info <- info %||% .load_runinfo(dir, lg)

  ss <- utils::read.csv(mqcpath(dir, "sheet"), stringsAsFactors = FALSE,
                        check.names = FALSE)
  cols <- cols %||% checkmeta(ss, logger = lg)
  platform <- unique(stats::na.omit(ss$detected_platform))[1] %||% NA_character_

  ## ---- per-sample metrics ------------------------------------------------
  qcm <- qcmetrics(s1$stats, detp = s1$detp, pthresh = cfg$detp,
                   samplemin = cfg$samplemin, intmad = cfg$intmad,
                   intfloor = cfg$intfloor, logger = lg)

  ## Attach the metadata columns BEFORE flagging. qcm is derived from the
  ## Stage 1 diagnostics and carries no sample sheet columns, so calling
  ## flagsamples() first would leave reported sex and age as NA and make
  ## sex_mismatch silently always FALSE.
  qcm <- .carry_meta(qcm, ss, cols)

  ## ---- frozen panels -----------------------------------------------------
  ## Computed BEFORE flagging, because the MDS outlier is one of the flags.
  sexp <- if (!is.na(platform)) sexprobes(platform, lg) else character(0)
  ## The matrix may have been collapsed to bare identifiers while sexprobes()
  ## returns whatever the manifest carries, so align the two before use.
  if (length(sexp) && !is.null(s1$probe_ids) && !has_v2_suffix(s1$probe_ids))
    sexp <- unique(stripv2(sexp))
  bt <- pcainput(s1$betas, s1$detp, s1$design, sexp, lg)
  pca <- runpca(bt, lg)
  mds <- runmds(bt, lg)
  dens <- betadensity(s1$betas)
  rm(bt); gc(verbose = FALSE)

  flags <- flagsamples(qcm, sexcol = cols$sex, agecol = cols$age,
                       mdsflag = mdsoutlier(mds), logger = lg)
  sex_threshold <- attr(flags, "sex_threshold")

  ## ---- epigenetic age (no masking) ---------------------------------------
  ages <- predictages(s1$betas, logger = lg)
  if (!is.null(ages)) flags <- .assign_cols(flags, ages, "sample_id")

  ## ---- cell composition ---------------------------------------------------
  props <- rundish(s1$betas, logger = lg)
  if (!is.null(props)) flags <- .assign_cols(flags, props, "sample_id")

  ## ---- cache -------------------------------------------------------------
  cc <- build_cache(s1, platform, pca = pca, mds = mds, density = dens)
  cc$design <- s1$design
  cc$sex_threshold <- sex_threshold
  save_rds_atomic(cc, mqcpath(dir, "cache", create = TRUE))

  ## ---- SNP identity ------------------------------------------------------
  snp <- snpbetas(s1$betas)
  donors <- if (!is.na(cols$donor) && cols$donor %in% names(ss))
    stats::setNames(as.character(ss[[cols$donor]]), ss$sample_id) else NULL
  sc <- snpcheck(snp, donors = donors, logger = lg)
  if (!is.null(sc))
    write_csv_atomic(sc, mqcpath(dir, "snpconc", create = TRUE))

  ## ---- report ------------------------------------------------------------
  pf <- probefail(s1$detp, cfg$detp)
  qcreport(flags, cc, pf, mqcpath(dir, "plots", create = TRUE),
           mqcpath(dir, "failedprobe", create = TRUE),
           mqcpath(dir, "pcscores"), cfg = cfg, logger = lg)

  ## ---- single sample sheet, updated in place -----------------------------
  ss <- .assign_cols(ss, flags, "sample_id")
  write_csv_atomic(ss, mqcpath(dir, "sheet"))

  info <- .extend_info(info, ss, flags, pca, props, ages, cfg, platform)
  save_rds_atomic(info, mqcpath(dir, "runinfo", create = TRUE))
  writemethods(dir, info)
  reportflags(flags, lg)
  lg$log("prep", "Stage 2 complete")
  invisible(dir)
}


#' Regenerate the threshold-dependent QC panels
#'
#' Recomputes only what a threshold change actually invalidates. PCA, MDS and
#' the beta density panels are frozen, so they never need recomputing; the
#' cached per-sample p-value histogram serves call rates at any threshold it
#' can resolve, and the cached per-probe grid counts serve probe failure rates
#' at any threshold on the grid. Only an off-grid detection threshold requires
#' re-reading the detection matrix.
#'
#' Without \code{suffix} the sample sheet is updated to match, so the report
#' and the sheet never disagree about how many samples failed. With
#' \code{suffix} nothing canonical is touched, so a side-by-side comparison
#' does not disturb the primary result.
#'
#' @param dir output directory.
#' @param detp,samplemin,failmin,intmad,intfloor,inclqual thresholds to change;
#'   \code{NULL} keeps the cached value.
#' @param suffix write \code{qc_plots_<suffix>.pdf} and matching sidecar files
#'   instead of overwriting the canonical ones.
#' @param dry print the plan and return without doing the work.
#' @return the PDF path, invisibly.
#' @export
#' @examples
#' \dontrun{
#' qcplots("~/qcout", detp = 0.01, dry = TRUE)
#' qcplots("~/qcout", detp = 0.01, suffix = "detp01")
#' }
qcplots <- function(dir, detp = NULL, samplemin = NULL, failmin = NULL,
                    intmad = NULL, intfloor = NULL, inclqual = NULL,
                    suffix = NULL, dry = FALSE) {

  mqccheckdir(dir)
  lg <- makelog(mqcpath(dir, "log", create = TRUE))
  on.exit(lg$close(), add = TRUE)

  cc <- loadcache(dir, lg)
  cfg <- mqcopts()
  new <- list(detp = detp, samplemin = samplemin, failmin = failmin,
              intmad = intmad, intfloor = intfloor, inclqual = inclqual)
  new <- new[!vapply(new, is.null, logical(1))]
  for (k in names(new)) cfg[[k]] <- new[[k]]

  pl <- plan_replot(cc, cfg$detp)
  msg <- c(sprintf("qcplots plan - %s", dir),
           if (length(new))
             sprintf("  changed: %s",
                     paste(sprintf("%s -> %s", names(new), unlist(new)),
                           collapse = ", "))
           else "  changed: nothing",
           sprintf("  %s", pl$note),
           "  frozen  : PCA, MDS, beta density",
           sprintf("  recompute: %s",
                   if (pl$needs_matrix) "detection matrix scan" else "nothing (cache only)"))
  for (m in msg) lg$log("qcplots", m)
  if (isTRUE(dry)) return(invisible(NULL))

  if (is.null(cc))
    stop("no usable QC cache in ", dir, "; run prep() first.", call. = FALSE)

  ## ---- per-sample metrics -------------------------------------------------
  ss <- utils::read.csv(mqcpath(dir, "sheet"), stringsAsFactors = FALSE,
                        check.names = FALSE)
  cols <- checkmeta(ss, logger = lg)

  ## The cached histogram cannot resolve a threshold finer than one bin, so
  ## fall back to the matrix rather than reporting a call rate of zero.
  cr_ok <- !is.null(callrate_from_hist(cc$phist, cfg$detp))
  dp <- NULL
  if (!cr_ok) {
    lg$log("qcplots", sprintf(
      "detp %g is finer than the cached histogram resolution; reading the detection matrix",
      cfg$detp))
    dp <- readRDS(mqcpath(dir, "detp"))
  }
  qcm <- qcmetrics(cc$stats, detp = dp, phist = if (cr_ok) cc$phist else NULL,
                   pthresh = cfg$detp, samplemin = cfg$samplemin,
                   intmad = cfg$intmad, intfloor = cfg$intfloor, logger = lg)
  qcm <- .carry_meta(qcm, ss, cols)
  flags <- flagsamples(qcm, sexcol = cols$sex, agecol = cols$age,
                       mdsflag = mdsoutlier(cc$mds), logger = lg)
  keep <- setdiff(names(ss), names(flags))
  for (k in keep) flags[[k]] <- ss[[k]][match(flags$sample_id, ss$sample_id)]

  ## ---- per-probe failure rate --------------------------------------------
  pf <- probefail_from_grid(cc$pfail, cc$n_samples, cfg$detp, cc$pgrid)
  if (is.null(pf)) {
    lg$log("qcplots", "reading the detection matrix for an off-grid threshold")
    dp <- dp %||% readRDS(mqcpath(dir, "detp"))
    pf <- probefail(dp, cfg$detp)
  }
  rm(dp); gc(verbose = FALSE)

  ## ---- destinations -------------------------------------------------------
  ## With a suffix, EVERY artefact is suffixed. 3.0.0 wrote a separate PDF but
  ## still clobbered failed_probes.csv and sample_sheet.csv, so the canonical
  ## report and its sidecar files described different thresholds.
  tag <- function(path, ext) if (is.null(suffix)) path else
    sub(paste0("\\", ext, "$"), paste0("_", suffix, ext), path)
  out <- if (is.null(suffix)) mqcpath(dir, "plots", create = TRUE) else
    file.path(mqcpath(dir, "qc", create = TRUE),
              sprintf("qc_plots_%s.pdf", suffix))

  qcreport(flags, cc, pf, out,
           tag(mqcpath(dir, "failedprobe", create = TRUE), ".csv"),
           tag(mqcpath(dir, "pcscores"), ".csv"), cfg = cfg, logger = lg)

  reportflags(flags, lg)
  if (is.null(suffix)) {
    ss <- .assign_cols(ss, flags, "sample_id")
    write_csv_atomic(ss, mqcpath(dir, "sheet"))
    lg$log("qcplots", "sample sheet updated to match the new thresholds")
  } else {
    lg$log("qcplots", sprintf(
      "suffixed run: %s written; the canonical report and sample sheet are untouched",
      basename(out)))
  }
  invisible(out)
}


#' Reassemble per-batch matrices written under extreme memory mode
#'
#' @param dir output directory.
#' @param what \code{"betas"} or \code{"detp"}.
#' @param write also write the assembled matrix to \code{data/matrices/}.
#' @return the assembled matrix.
#' @export
#' @examples
#' \dontrun{
#' b <- makemat("~/qcout", "betas", write = TRUE)
#' }
makemat <- function(dir, what = c("betas", "detp"), write = FALSE) {
  what <- match.arg(what)
  bd <- mqcpath(dir, "batches")
  if (!dir.exists(bd)) stop("no batch directory at ", bd, call. = FALSE)
  pat <- if (identical(what, "betas")) "^betas_[0-9]+\\.rds$" else "^detP_[0-9]+\\.rds$"
  files <- sort(list.files(bd, pattern = pat, full.names = TRUE))
  if (!length(files)) stop("no ", what, " batch files in ", bd, call. = FALSE)

  ## Only accept batch files belonging to the CURRENT Stage 1. A shorter re-run
  ## into the same directory leaves higher-numbered files from the previous
  ## one, and 3.0.0 cbind()ed them in, producing a matrix with phantom columns
  ## from samples that were not in this cohort at all.
  ids <- .expected_samples(dir)
  if (!is.null(ids)) {
    ok <- vapply(files, function(f) {
      cn <- tryCatch(colnames(readRDS(f)), error = function(e) NULL)
      !is.null(cn) && all(cn %in% ids)
    }, logical(1))
    if (any(!ok)) {
      warning(sum(!ok), " batch file(s) in ", bd, " contain samples that are ",
              "not in this run's sample sheet and were ignored as stale: ",
              paste(basename(files[!ok]), collapse = ", "), call. = FALSE)
      files <- files[ok]
    }
    if (!length(files))
      stop("every batch file in ", bd, " is stale; re-run qc().", call. = FALSE)
  }

  first <- readRDS(files[1])
  rid <- rownames(first)
  parts <- vector("list", length(files))
  parts[[1]] <- first
  for (k in seq_along(files)[-1]) {
    m <- readRDS(files[k])
    if (!identical(rownames(m), rid)) m <- m[match(rid, rownames(m)), , drop = FALSE]
    parts[[k]] <- m
  }
  out <- do.call(cbind, parts)
  rownames(out) <- rid
  if (anyDuplicated(colnames(out)))
    stop("duplicate sample columns after reassembly; the batch directory ",
         "mixes runs. Clear ", bd, " and re-run qc().", call. = FALSE)
  if (isTRUE(write))
    save_rds_atomic(out, mqcpath(dir, if (identical(what, "betas")) "betas" else "detp",
                                 create = TRUE))
  out
}


## ---------------------------------------------------------------------------
## Internals
## ---------------------------------------------------------------------------

## The sample IDs this output directory is supposed to contain.
.expected_samples <- function(dir) {
  p <- mqcpath(dir, "sheet")
  if (!file.exists(p)) return(NULL)
  ss <- tryCatch(utils::read.csv(p, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(ss) || !"sample_id" %in% names(ss)) return(NULL)
  as.character(ss$sample_id)
}

## Copy the resolved metadata columns onto the metrics frame, so that
## flagsamples() can see reported sex and age.
.carry_meta <- function(qcm, ss, cols) {
  ## The resolved metadata columns, plus the chip-position columns discover()
  ## parsed from the Sentrix barcode and anything that looks like a plate or
  ## well. These all end up on the PC/metadata pages, which cannot show a
  ## variable that never reached this frame.
  want <- stats::na.omit(unlist(cols[c("sex", "age", "batch", "cell", "donor")]))
  extra <- c(.MQC_POSITION_COLS, "Plate", "Sample_Plate", "Sample_Group",
             "Well", "Sample_Well", "Well_Position", "Sentrix_Position",
             "Slide", "Chip", "Sentrix_ID")
  want <- c(as.character(want), intersect(names(ss), extra))
  want <- intersect(unique(want), names(ss))
  want <- setdiff(want, names(qcm))
  if (!length(want)) return(qcm)
  .assign_cols(qcm, ss[, c("sample_id", want), drop = FALSE], "sample_id")
}

## A memory cap derived from the system even when no plan was run.
.default_memcap <- function(cfg) {
  a <- sysmem()$available
  if (is.na(a)) NULL else a * cfg$memfrac
}

## Assign columns from `src` onto `dst`, matched on `key`. Assignment rather
## than a join, so re-running never produces .x/.y duplicates.
.assign_cols <- function(dst, src, key) {
  if (is.null(src) || !nrow(src)) return(dst)
  if (!key %in% names(dst) || !key %in% names(src))
    stop("both frames need a '", key, "' column", call. = FALSE)
  m <- match(dst[[key]], src[[key]])
  for (cl in setdiff(names(src), key)) dst[[cl]] <- src[[cl]][m]
  dst
}

.decide_collapse <- function(collapse, cfg, probe_ids, lg) {
  want <- if (!is.null(collapse)) collapse else cfg$collapse
  auto <- has_v2_suffix(probe_ids)
  if (length(want) != 1L || is.na(want)) {
    if (auto) lg$log("qc", "EPICv2 replicate suffixes detected; collapsing")
    return(auto)
  }
  if (isTRUE(want) && !auto) {
    lg$log("qc", "collapse requested but no replicate suffixes are present; skipping")
    return(FALSE)
  }
  if (!isTRUE(want) && auto)
    lg$log("qc", paste("EPICv2 replicate suffixes are present but collapse is",
                       "disabled. Clock and cell-type references key on bare",
                       "cg identifiers and will not match."), warn = TRUE)
  isTRUE(want)
}

.load_runinfo <- function(dir, lg) {
  p <- mqcpath(dir, "runinfo")
  if (!file.exists(p)) {
    lg$log("prep", paste("no run_info.rds; METHODS.txt will omit the Stage 1",
                         "facts. Re-run qc() to regenerate them."), warn = TRUE)
    return(NULL)
  }
  tryCatch(readRDS(p), error = function(e) NULL)
}

.load_stage1 <- function(dir, lg) {
  betas <- readRDS(mqcpath(dir, "betas"))
  detp <- readRDS(mqcpath(dir, "detp"))
  design <- readRDS(mqcpath(dir, "designmask"))
  sdfs <- if (file.exists(mqcpath(dir, "sdfs"))) readRDS(mqcpath(dir, "sdfs")) else NULL
  cc <- loadcache(dir, lg)
  ss <- utils::read.csv(mqcpath(dir, "sheet"), stringsAsFactors = FALSE,
                        check.names = FALSE)
  stats <- cc$stats
  if (is.null(stats)) {
    keepc <- intersect(names(ss), c("sample_id", names(stats_raw_template())))
    stats <- ss[, keepc, drop = FALSE]
  }
  detp <- conform_to(detp, betas, "detP matrix")
  lg$log("prep", sprintf("loaded Stage 1: %d probes x %d samples",
                         nrow(betas), ncol(betas)))
  list(betas = betas, detp = detp, design = design, stats = stats,
       sdfs = sdfs, probe_ids = rownames(betas),
       phist = cc$phist, pfail = cc$pfail, pgrid = cc$pgrid %||% .MQC_PGRID)
}

.extend_info <- function(info, ss, flags, pca, props, ages, cfg, platform) {
  info <- info %||% list()
  info$platform <- info$platform %||% platform
  info$n_samples <- nrow(ss)
  info$n_low_callrate <- sum(flags$low_callrate, na.rm = TRUE)
  info$n_low_intensity <- sum(flags$low_intensity, na.rm = TRUE)
  info$n_sex_mismatch <- sum(flags$sex_mismatch, na.rm = TRUE)
  info$n_sex_unclear <- sum(flags$sex_unclear, na.rm = TRUE)
  info$sex_called <- any(!is.na(flags$inferred_sex))
  info$int_cutoff <- sprintf("%.0f", unique(stats::na.omit(flags$intensity_cutoff))[1] %||% NA_real_)
  info$int_median <- sprintf("%.0f", stats::median(flags$mean_intensity_raw, na.rm = TRUE))
  info$cohort_low <- isTRUE(stats::median(flags$mean_intensity_raw, na.rm = TRUE) < cfg$intfloor)
  info$n_pca_probes <- pca$n_probes %||% NA
  info$cells_ok <- !is.null(props)
  ## These are counts, not objects. 3.0.0 stored the whole data.frame and then
  ## interpolated it into sprintf("%s"), which emitted one copy of the
  ## paragraph per column, and reported the number of cell TYPES under the
  ## label "reference probes matched the array".
  info$n_cell_types <- if (is.null(props)) NA else
    (attr(props, "n_types") %||% (ncol(props) - 1L))
  info$n_dish_probes <- if (is.null(props)) NA else (attr(props, "n_probes") %||% NA)
  info$clock_ok <- !is.null(ages)
  info$clock <- if (is.null(ages)) NA_character_ else "Horvath353"
  if (!is.null(ages)) {
    np <- grep("nprobes$", names(ages), value = TRUE)[1]
    nm <- grep("nmodel$", names(ages), value = TRUE)[1]
    info$n_clock_present <- if (is.na(np)) NA else max(ages[[np]], na.rm = TRUE)
    info$n_clock_model <- if (is.na(nm)) NA else ages[[nm]][1]
  }
  info$maskuse <- cfg$maskuse; info$detp <- cfg$detp
  info$samplemin <- cfg$samplemin; info$intmad <- cfg$intmad
  info$intfloor <- cfg$intfloor
  info$dishref <- cfg$dishref; info$dishmethod <- cfg$dishmethod
  info
}

## The end-of-run summary pipeline() prints.
.run_summary <- function(dir) {
  L <- character(0)
  add <- function(...) L <<- c(L, sprintf(...))
  info <- tryCatch(readRDS(mqcpath(dir, "runinfo")), error = function(e) NULL)
  ss <- tryCatch(utils::read.csv(mqcpath(dir, "sheet"), stringsAsFactors = FALSE),
                 error = function(e) NULL)
  add("")
  add("methylQC complete - %s", dir)
  add("")
  if (!is.null(info)) {
    add("  %-26s: %s, %s probes", "Platform", info$platform %||% "?",
        format(info$n_probes %||% NA, big.mark = ","))
    add("  %-26s: %s processed, %s failed", "Samples",
        format(info$n_samples %||% NA, big.mark = ","), info$n_failed %||% 0)
    add("  %-26s: %s", "Stage 1 time", info$stage1_time %||% "?")
  }
  if (!is.null(ss)) {
    f <- function(k) if (k %in% names(ss)) sum(ss[[k]], na.rm = TRUE) else NA
    add("  %-26s: %s", "Flagged samples", f("flagged"))
    add("      low call rate         : %s", f("low_callrate"))
    add("      low intensity         : %s", f("low_intensity"))
    add("      sex mismatch          : %s", f("sex_mismatch"))
    add("      MDS outlier           : %s", f("mds_outlier"))
    if (isFALSE(info$sex_called))
      add("      (no sex was called: cohort is not bimodal in chrY intensity)")
  }
  add("")
  add("  Report    : %s", mqcpath(dir, "plots"))
  add("  Sheet     : %s", mqcpath(dir, "sheet"))
  add("  Methods   : %s", mqcpath(dir, "methods"))
  add("  Log       : %s", mqcpath(dir, "log"))
  add("")
  add("  Retune a threshold without reprocessing:")
  add("    qcplots(\"%s\", detp = 0.01)", dir)
  add("")
  L
}
