## ---------------------------------------------------------------------------
## metrics.R -- sample-level quality metrics and thresholds.
##
## INTENSITY THRESHOLDS
## --------------------
## v2 used a fixed cutoff (intmin = 1300) on mean intensity. A fixed cutoff
## has two problems: it depends on scanner, chemistry, DNA input and array
## version, and it was calibrated against intensity measured after the full
## QCDPB chain, whereas v3 measures intensity before dye bias and noob. Those
## are different quantities.
##
## v3 flags per-sample outliers relative to the cohort, using the median
## absolute deviation of log2 intensity. Log first because intensity is
## right-skewed; one-sided because only LOW intensity indicates failure.
##
## A relative rule cannot detect whole-cohort failure -- if every array ran
## low, the median moves with them and nothing is flagged. Simulation at
## n = 200: with three genuinely bad arrays in a healthy cohort, the fixed
## rule and the MAD rule both flagged exactly those three; with the entire
## cohort degraded, the fixed rule flagged 176 samples and the MAD rule
## flagged none. Neither answer alone is right, so v3 keeps an absolute floor
## as a COHORT-LEVEL WARNING (is this whole run degraded?) while individual
## flags come from the MAD rule.
## ---------------------------------------------------------------------------

#' Flag low-intensity outliers within a cohort
#'
#' @param intensity numeric vector of per-sample mean intensity.
#' @param nmad MAD multiplier; 3 is the conventional "extreme" line.
#' @param floor absolute value below which the whole cohort is reported as
#'   degraded. Set \code{NA} to disable.
#' @return a list with \code{flag} (logical vector), \code{cutoff} (on the
#'   original scale), \code{cohort_low} and \code{median}.
#' @export
#' @examples
#' x <- c(rlnorm(50, log(4000), 0.2), 700)
#' intflag(x)$flag[51]
intflag <- function(intensity, nmad = NULL, floor = NULL) {
  cfg <- mqcopts()
  nmad <- .opt(nmad, "intmad", cfg)
  floor <- if (missing(floor)) cfg$intfloor else floor

  x <- as.numeric(intensity)
  ok <- is.finite(x) & x > 0
  out <- list(flag = rep(FALSE, length(x)), cutoff = NA_real_,
              cohort_low = FALSE, median = NA_real_, n_usable = sum(ok))

  ## Too few usable samples to estimate a spread. The relative rule is off, but
  ## a sample with no measurable intensity is still a failure -- 3.0.0 returned
  ## all-FALSE here, passing unmeasurable arrays silently.
  if (sum(ok) < 5L) {
    out$flag <- !ok
    return(out)
  }

  lx <- log2(x[ok])
  med <- stats::median(lx)
  s <- stats::mad(lx)
  out$median <- 2^med

  if (is.finite(s) && s > 0) {
    cut_log <- med - nmad * s
    out$cutoff <- 2^cut_log
    f <- rep(FALSE, length(x))
    f[ok] <- lx < cut_log
    f[!ok] <- TRUE                       # unmeasurable intensity is a failure
    out$flag <- f
  } else {
    out$flag <- !ok
  }

  if (!is.na(floor) && is.finite(floor) && out$median < floor)
    out$cohort_low <- TRUE
  out
}

#' Assemble per-sample QC metrics
#'
#' Joins the Stage 1 diagnostics with detection-derived metrics and applies
#' the sample-level thresholds.
#'
#' @param stats per-sample diagnostics from \code{\link{runsesame}}.
#' @param detp detection p-value matrix, or \code{NULL} when \code{phist} is
#'   supplied instead.
#' @param phist cached per-sample p-value histogram (bins x samples).
#' @param pthresh detection p-value threshold.
#' @param samplemin minimum acceptable call rate.
#' @param intmad MAD multiplier for the intensity outlier rule.
#' @param intfloor absolute cohort-level intensity floor.
#' @param logger optional logger.
#' @return a data.frame, one row per sample.
#' @export
#' @examples
#' \dontrun{
#' qm <- qcmetrics(out$stats, out$detp)
#' }
qcmetrics <- function(stats, detp = NULL, phist = NULL,
                      pthresh = NULL, samplemin = NULL,
                      intmad = NULL, intfloor = NULL, logger = NULL) {

  logger <- logger %||% nulllog()
  cfg <- mqcopts()
  pthresh <- .opt(pthresh, "detp", cfg)
  samplemin <- .opt(samplemin, "samplemin", cfg)
  intmad <- .opt(intmad, "intmad", cfg)
  intfloor <- .opt(intfloor, "intfloor", cfg)

  qc <- stats
  if (!"sample_id" %in% names(qc))
    stop("stats must carry a sample_id column", call. = FALSE)

  cr <- if (!is.null(detp)) callrate(detp, pthresh) else
          callrate_from_hist(phist, pthresh)
  if (is.null(cr))
    stop("cannot compute call rate: supply detp, or a phist whose bin width ",
         "(", if (is.null(phist)) "no histogram" else 1 / nrow(phist),
         ") resolves a threshold of ", pthresh, call. = FALSE)
  qc$call_rate <- unname(cr[match(qc$sample_id, names(cr))])
  qc$detp_threshold <- pthresh

  ## Fail closed: a sample whose call rate could not be computed is failed,
  ## not silently passed. v2's `cr < samplemin` produced NA for such samples,
  ## and the NA then propagated into `out[out$flagged, ]`, emitting a phantom
  ## all-NA row into the CSV while never flagging the real sample.
  qc$low_callrate <- is.na(qc$call_rate) | qc$call_rate < samplemin

  ints <- intflag(qc$mean_intensity_raw, nmad = intmad, floor = intfloor)
  qc$low_intensity <- ints$flag
  qc$intensity_cutoff <- ints$cutoff

  if (isTRUE(ints$cohort_low))
    logger$log("qc", sprintf(
      paste("COHORT WARNING: median intensity %.0f is below the floor of %.0f.",
            "The whole run looks degraded; per-sample flags are relative to",
            "this cohort and will not detect it."),
      ints$median, intfloor), warn = TRUE)

  logger$log("qc", sprintf(
    "call rate: median %.4f, %d below %.3f | intensity: median %.0f, cutoff %.0f, %d flagged",
    stats::median(qc$call_rate, na.rm = TRUE), sum(qc$low_callrate),
    samplemin, ints$median %||% NA_real_, ints$cutoff %||% NA_real_,
    sum(qc$low_intensity)))

  qc
}

#' Call rate recovered from a cached p-value histogram
#'
#' The histogram is built once during Stage 1 at about 0.1 s per sample,
#' inside the worker where the p-values are already in hand. It makes the call
#' rate available at any threshold without re-reading the detection matrix,
#' which at 1500 EPIC samples is a 9.7 GiB file.
#'
#' @param phist bins-by-samples integer matrix.
#' @param pthresh threshold.
#' @return a named numeric vector of call rates.
#' @export
#' @examples
#' h <- matrix(c(rep(0L, 999), 100L, rep(50L, 1000)), 1000, 2,
#'             dimnames = list(NULL, c("a", "b")))
#' callrate_from_hist(h, 0.5)
callrate_from_hist <- function(phist, pthresh = NULL) {
  pthresh <- .opt(pthresh, "detp")
  if (is.null(phist)) return(NULL)
  nb <- nrow(phist)
  ## Bin k covers ((k-1)/nb, k/nb], so summing bins 1..k counts exactly the
  ## p <= k/nb that callrate() counts on the matrix.
  k <- floor(pthresh * nb + 1e-9)
  ## A threshold finer than one bin cannot be answered from the histogram at
  ## all. 3.0.0 returned a call rate of 0 for every sample here and wrote that
  ## to the sample sheet; returning NULL makes the caller re-read the matrix,
  ## which is what probefail_from_grid() already does off-grid.
  if (k < 1L) return(NULL)
  k <- min(as.integer(k), nb)
  tot <- colSums(phist)
  hit <- colSums(phist[seq_len(k), , drop = FALSE])
  stats::setNames(ifelse(tot > 0, hit / tot, NA_real_), colnames(phist))
}

#' Per-probe failure rate recovered from cached grid counts
#'
#' @param pfail probes-by-grid integer matrix of counts of samples exceeding
#'   each cached threshold.
#' @param n_samples number of samples.
#' @param pthresh requested threshold.
#' @param grid the cached threshold grid.
#' @return a named numeric vector, or \code{NULL} if the threshold is not on
#'   the grid (the caller must then re-read the detection matrix).
#' @export
#' @examples
#' pf <- matrix(c(0L, 10L), 2, 1, dimnames = list(c("cg1", "cg2"), "p0.05"))
#' probefail_from_grid(pf, 10, 0.05, 0.05)
probefail_from_grid <- function(pfail, n_samples, pthresh, grid) {
  j <- which(abs(grid - pthresh) < 1e-12)
  if (!length(j)) return(NULL)
  stats::setNames(pfail[, j[1]] / n_samples, rownames(pfail))
}
