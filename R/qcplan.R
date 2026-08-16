## ---------------------------------------------------------------------------
## qcplan.R -- preflight sizing.
##
## Goal:    decide worker count, batch size and whether SigDFs can be retained,
##          before a single IDAT is read.
## Approach: the parent-side requirement is arithmetic, not estimation --
##          n_probes * n_samples * 8 bytes per retained double matrix, and both
##          terms are known in advance. Only the per-sample worker cost needs
##          measuring, and it is measured by processing one sample (which the
##          real run has to do anyway) bracketed by gc(reset = TRUE).
## Inputs:  a discovered sample sheet, or a directory to discover.
## Outputs: a plan list, printed as a short report.
## Usage:   plan <- qcplan("~/idats")
## ---------------------------------------------------------------------------

## Bytes per probe for one retained SigDF. A SigDF is a data.frame of
## Probe_ID (character, 8 bytes of pointer once the string pool is shared
## across samples), MG/MR/UG/UR (4 x 8 bytes double), col (factor, 4) and mask
## (logical, 4): 8 + 32 + 4 + 4 = 48. Measured marginal cost at 866,553 EPIC
## probes is 34 MB, i.e. 41 bytes/probe; 48 is used so the guard errs towards
## disabling retention rather than towards an out-of-memory abort.
.MQC_SIGDF_BYTES_PER_PROBE <- 48

## Rough compression ratio of a SigDF under saveRDS()'s default gzip, measured
## at 20 MB written per 34 MB resident.
.MQC_SIGDF_DISK_RATIO <- 0.6

#' Bytes one retained SigDF occupies
#'
#' @param n_probes probes on the array.
#' @return bytes.
#' @keywords internal
#' @noRd
sigdfbytes <- function(n_probes) as.numeric(n_probes) * .MQC_SIGDF_BYTES_PER_PROBE

#' Preflight a methylQC run
#'
#' Measures the actual per-sample memory cost on this machine with this array,
#' reads available system memory, and recommends a worker count, a batch size
#' and whether the preprocessed \code{SigDF} objects can be retained.
#'
#' \code{qc()} calls this automatically when \code{workers} and \code{batch}
#' are not given, so the defaults are self-derived. Any argument you pass to
#' \code{qc()} overrides the corresponding recommendation.
#'
#' @param indir directory of IDAT files, or \code{NULL} if \code{ss} is given.
#' @param ss an already-discovered sample sheet, to avoid re-scanning.
#' @param platform platform string; inferred when \code{NULL}.
#' @param memfrac fraction of available memory the run may use.
#' @param extreme plan for extreme memory conservation mode.
#' @param savesdf whether SigDF retention is wanted; the plan may withdraw it.
#' @param outdir destination directory, used only to check free disk space.
#' @param measure if \code{FALSE}, skip processing a sample and use a
#'   conservative default for per-sample cost.
#' @param verbose print the report.
#' @return a list with \code{workers}, \code{batch}, \code{savesdf},
#'   \code{n_samples}, \code{n_probes}, \code{per_sample}, \code{parent_bytes},
#'   \code{projected_peak}, \code{memcap}, \code{parentcap}, \code{fits} and
#'   \code{text}.
#' @export
#' @examples
#' \dontrun{
#' plan <- qcplan("~/idats")
#' qc("~/idats", "~/qcout", workers = plan$workers)
#' }
qcplan <- function(indir = NULL, ss = NULL, platform = NULL,
                   memfrac = NULL, extreme = NULL, savesdf = NULL,
                   outdir = NULL, measure = TRUE, verbose = TRUE) {

  cfg <- mqcopts()
  memfrac <- .opt(memfrac, "memfrac", cfg)
  extreme <- .opt(extreme, "extreme", cfg)
  savesdf <- .opt(savesdf, "savesdf", cfg)

  if (is.null(ss)) {
    if (is.null(indir)) stop("supply either indir or ss", call. = FALSE)
    ss <- discover(indir, logger = nulllog())
  }
  n <- nrow(ss)
  platform <- platform %||%
    unique(stats::na.omit(ss$detected_platform))[1] %||% NA_character_
  if (is.na(platform))
    stop("could not determine the array platform; pass platform= explicitly",
         call. = FALSE)

  mem <- sysmem()
  avail <- mem$available
  budget <- if (is.na(avail)) NA_real_ else avail * memfrac

  ## ---- measure one sample -------------------------------------------------
  per_sample <- NA_real_; np <- NA_integer_; t_sample <- NA_real_
  if (isTRUE(measure)) {
    m <- tryCatch(.measure_one(ss$Basename[1], platform),
                  error = function(e) NULL)
    if (!is.null(m)) {
      per_sample <- m$bytes; np <- m$n_probes; t_sample <- m$seconds
    }
  }
  measured <- !is.na(per_sample)
  if (is.na(np)) np <- .probe_guess(platform)
  if (is.na(per_sample)) per_sample <- 400 * 2^20     # conservative fallback
  if (is.na(t_sample)) t_sample <- 6

  ## ---- parent-side matrices: exact arithmetic -----------------------------
  n_double <- if (isTRUE(extreme)) 0L else 2L         # betas + detP
  fixed <- as.numeric(np) * n * 8 * n_double
  fixed <- fixed + as.numeric(np) * 4                          # design mask
  fixed <- fixed + as.numeric(np) * length(.MQC_PGRID) * 4     # grid counts
  fixed <- fixed + .MQC_PBINS * n * 4                          # p-value hist

  ## The SigDF list is the single largest parent allocation and, unlike worker
  ## memory, it grows monotonically and is never released between chunks.
  sdf_bytes <- if (isTRUE(savesdf)) sigdfbytes(np) * n else 0
  parent_bytes <- fixed + sdf_bytes

  ## ---- withdraw SigDF retention if it does not fit ------------------------
  sdf_note <- character(0)
  savesdf_downgraded <- FALSE
  if (isTRUE(savesdf) && !is.na(budget)) {
    room_needed <- fixed + sdf_bytes + per_sample     # at least one worker
    if (room_needed > budget) {
      savesdf <- FALSE
      savesdf_downgraded <- TRUE
      parent_bytes <- fixed
      sdf_note <- .sdf_downgrade_text(sdf_bytes, room_needed, budget, n, np)
    }
  }

  ## ---- disk guard for the sdfs_all.rds write ------------------------------
  if (isTRUE(savesdf) && !is.null(outdir)) {
    need <- sdf_bytes * .MQC_SIGDF_DISK_RATIO
    free <- diskfree(outdir)
    if (!is.na(free) && need > free) {
      savesdf <- FALSE
      savesdf_downgraded <- TRUE
      parent_bytes <- fixed
      sdf_note <- c(sprintf(
        "  SigDF retention disabled: sdfs_all.rds needs about %s but only %s is free on the output volume.",
        fmtbytes(need), fmtbytes(free)),
        "  Nothing else is affected; sex, age, cell composition and every QC panel are unchanged.")
    }
  }

  ## ---- workers ------------------------------------------------------------
  ## detectCores() reports failure as NA, not NULL, so %||% let NA through into
  ## the arithmetic and workers came out NA.
  cores <- max(1L, as.integer(parallel::detectCores(logical = FALSE) %|NA|% 1L))
  ceiling_workers <- max(1L, cores - 1L)
  if (is.na(budget)) {
    workers <- 1L
  } else {
    room <- budget - parent_bytes
    workers <- if (room <= 0) 1L else
      max(1L, min(ceiling_workers, floor(room / per_sample)))
  }
  if (isTRUE(extreme)) workers <- min(workers, 2L)
  if (!can_fork()) workers <- 1L
  workers <- max(1L, as.integer(workers))

  ## ---- batch size ---------------------------------------------------------
  ## Memory scales linearly with workers * batch; time is flat because
  ## per-task overhead is milliseconds against seconds of per-sample work.
  ## Keep >= 4 tasks per worker so a slow sample cannot strand the tail.
  payload <- as.numeric(np) * 8 * 2                     # betas + pvals, 1 sample
  by_mem <- if (is.na(budget)) 25 else
    floor(max(1, (budget - parent_bytes) * 0.5) / (workers * payload))
  by_tail <- floor(n / max(1L, 4L * workers))
  batch <- max(1L, min(50L, by_mem, if (by_tail >= 1) by_tail else 1L, n))

  projected <- parent_bytes + workers * per_sample
  fits <- !is.na(budget) && projected <= budget

  ## The guard inside runsesame() can only see this process, so it is given the
  ## budget with the worker allowance already removed.
  parentcap <- if (is.na(budget)) NA_real_ else budget - workers * per_sample

  plan <- list(
    platform = platform, n_samples = n, n_probes = np,
    workers = workers, batch = as.integer(batch),
    savesdf = isTRUE(savesdf), savesdf_downgraded = savesdf_downgraded,
    sigdf_bytes = sdf_bytes,
    per_sample = per_sample, measured = measured,
    parent_bytes = parent_bytes,
    projected_peak = projected, memcap = budget, parentcap = parentcap,
    total_mem = mem$total, avail_mem = avail, memfrac = memfrac,
    cores = cores, extreme = isTRUE(extreme),
    seconds_per_sample = t_sample,
    projected_seconds = t_sample * n / max(1L, workers),
    fits = fits, sdf_note = sdf_note)

  plan$text <- .plan_text(plan)
  if (isTRUE(verbose)) cat(plan$text, sep = "\n")
  invisible(plan)
}

## A quiet, cheap variant used internally when qc() needs defaults but the
## caller did not run qcplan() explicitly.
#' @keywords internal
#' @noRd
qcplan_quiet <- function(n_samples) {
  mem <- sysmem()
  cores <- max(1L, as.integer(parallel::detectCores(logical = FALSE) %|NA|% 1L))
  w <- if (is.na(mem$available)) 1L else
    max(1L, min(cores - 1L, floor(mem$available * 0.5 / (400 * 2^20))))
  if (!can_fork()) w <- 1L
  list(workers = as.integer(min(w, max(1L, n_samples))))
}

.measure_one <- function(pfx, platform) {
  mqcprewarm(nulllog())
  addr <- mqcanno(platform, nulllog())
  gc(reset = TRUE, verbose = FALSE)
  t0 <- Sys.time()
  r <- process_one(pfx, platform, addr, keep_sdf = FALSE, want_design = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.na(r$error)) stop(r$error, call. = FALSE)
  list(bytes = .gc_peak_bytes(), n_probes = length(r$betas), seconds = secs)
}

## gc()'s column layout is not fixed: R inserts a "limit (Mb)" column when a
## memory limit is in force, which on this build makes "max used" column 6 and
## its megabyte figure column 7. 3.0.0 hard-coded column 6 and multiplied raw
## CELL COUNTS by 2^20, overstating per-sample cost by about six orders of
## magnitude -- enough that workers was pinned to 1 and every plan reported
## INSUFFICIENT MEMORY. Locate the column by name instead.
.gc_peak_bytes <- function() {
  g <- gc(verbose = FALSE)
  j <- which(colnames(g) == "max used")
  if (!length(j) || (j[1] + 1L) > ncol(g)) return(NA_real_)
  sum(g[, j[1] + 1L], na.rm = TRUE) * 2^20
}

.probe_guess <- function(platform) {
  switch(as.character(platform),
         HM450 = 485512L, EPIC = 866553L, EPICv2 = 937055L,
         MSA = 269094L, MM285 = 296070L, 866553L)
}

## What the user is told when SigDF retention is withdrawn. It has to say what
## was disabled, why, by how much, and what to do -- and that nothing else
## changes, because after 3.0.1 the sex check no longer depends on SigDFs.
.sdf_downgrade_text <- function(sdf_bytes, need, budget, n, np) {
  c(sprintf("  SigDF retention disabled: keeping them for %d samples needs %s,",
            n, fmtbytes(sdf_bytes)),
    sprintf("  which takes the run to %s against a budget of %s (short by %s).",
            fmtbytes(need), fmtbytes(budget), fmtbytes(need - budget)),
    "  Nothing else is affected: sex calling, epigenetic age, cell composition",
    "  and every QC panel are computed without them.",
    "  To keep them anyway, do one of:",
    sprintf("    - process in runs of about %d samples",
            max(1L, floor((budget * 0.8) / max(1, sigdfbytes(np))))),
    "    - raise the share of RAM methylQC may use, e.g. mqcset(memfrac = 0.9)",
    "    - free memory, then re-run")
}

.plan_text <- function(p) {
  L <- character(0)
  add <- function(...) L <<- c(L, sprintf(...))

  add("")
  add("methylQC plan - %s, %d samples, %s probes",
      p$platform, p$n_samples, format(p$n_probes, big.mark = ","))
  add("")
  add("  %-30s: %s", if (isTRUE(p$measured)) "Measured per-sample peak"
                     else "Assumed per-sample peak (unmeasured)",
      fmtbytes(p$per_sample))
  add("  %-30s: %s (of %s)", "Available memory",
      fmtbytes(p$avail_mem), fmtbytes(p$total_mem))
  add("  %-30s: %s", "Parent matrices", fmtbytes(p$parent_bytes))
  add("  %-30s: %s", sprintf("Budget (memfrac %.2f)", p$memfrac),
      fmtbytes(p$memcap))
  add("")
  add("  %-30s: %d   (of %d physical cores)", "Recommended workers",
      p$workers, p$cores)
  add("  %-30s: %d", "Recommended batch size", p$batch)
  add("  %-30s: %s", "Retain SigDFs (savesdf)",
      if (isTRUE(p$savesdf)) "yes" else
        if (isTRUE(p$savesdf_downgraded)) "no  (disabled by the memory check)" else "no")
  add("  %-30s: %s   %s", "Projected peak",
      fmtbytes(p$projected_peak), if (isTRUE(p$fits)) "OK" else "EXCEEDS BUDGET")
  add("  %-30s: ~%s", "Projected Stage 1 time", fmtdur(p$projected_seconds))
  add("")
  if (length(p$sdf_note)) { L <- c(L, p$sdf_note, "") }

  if (isTRUE(p$fits)) {
    add("  Run with:  pipeline(indir, outdir)            # inherits the above")
    add("             pipeline(indir, outdir, workers = %d)   # override one field",
        max(1L, p$workers - 1L))
  } else {
    short <- p$projected_peak - (p$memcap %|NA|% 0)
    nruns <- max(2L, ceiling(p$projected_peak / max(1, p$memcap %|NA|% 1)))
    add("  INSUFFICIENT MEMORY. Options, cheapest first:")
    add("    1. pipeline(..., extreme = TRUE)")
    add("    2. Process in %d runs of about %d samples",
        nruns, max(1L, floor(p$n_samples / nruns)))
    add("    3. Add about %s of RAM", fmtbytes(short))
  }
  add("")
  L
}
