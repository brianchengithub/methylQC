## ---------------------------------------------------------------------------
## cache.R -- derived summaries that let qcplots() retune thresholds without
## re-reading the matrices or recomputing anything expensive.
##
## Measured at EPIC scale (866,000 probes), extrapolated to 1500 samples:
##
##   sample call rate   colMeans      ~9 s      linear in samples
##   probe fail rate    rowMeans     ~11 s      linear in samples
##   prcomp, 100k probes            ~15 min     QUADRATIC in samples
##
## plus reading two 9.7 GiB matrices from disk. PCA is roughly 98% of the
## work, so the cache is designed around removing it from the loop entirely
## rather than around the matrix scans.
##
## What is threshold-dependent and therefore regenerated:
##     call rate vs sample (detp, samplemin), intensity (intmad, intfloor),
##     probe failure (detp, failmin), sex, age
##
## What is threshold-independent and therefore frozen:
##     PCA, MDS, beta density
##
## PCA is frozen legitimately rather than merely reused, because its input
## probe set is defined without reference to any user-tunable threshold: see
## pcainput() in plots.R.
##
## Note on storage: the prcomp object's `rotation` component is
## n_probes x n_PCs, which at 100,000 probes and 1500 samples is 1.2 GiB.
## The plots need only the scores, the standard deviations and the probe
## count, so only those are kept.
## ---------------------------------------------------------------------------

.MQC_CACHE_VERSION <- 3L

#' Build the QC cache
#'
#' @param s1 the list returned by \code{\link{runsesame}}.
#' @param platform platform string.
#' @param pca list with \code{sdev} and \code{n_probes}, or \code{NULL}.
#' @param mds matrix of MDS coordinates, or \code{NULL}.
#' @param density list of per-sample density curves, or \code{NULL}.
#' @return a cache list.
#' @keywords internal
#' @noRd
build_cache <- function(s1, platform, pca = NULL, mds = NULL, density = NULL) {
  ## phist and pfail come free from Stage 1, but only when Stage 1 ran in this
  ## process. prep(dir) loads Stage 1 from disk and gets them from the cache --
  ## which does not exist the first time prep() runs -- so 3.0.0 wrote a cache
  ## with phist = NULL, and every later qcplots() died with "supply either detp
  ## or phist". Worse, re-running prep() on a good directory replaced a working
  ## cache with a broken one. Derive them from the detection matrix instead
  ## when they are absent; it is one pass over a matrix already in memory.
  ph <- s1$phist %||% .phist_from_detp(s1$detp)
  pf <- s1$pfail %||% .pfail_from_detp(s1$detp, s1$pgrid %||% .MQC_PGRID)
  ns <- if (!is.null(ph)) ncol(ph) else
    if (!is.null(s1$betas)) ncol(s1$betas) else NA_integer_

  list(
    cache_version = .MQC_CACHE_VERSION,
    pkg_version   = as.character(utils::packageVersion("methylQC")),
    created       = as.character(Sys.time()),
    platform      = platform,
    n_probes      = length(s1$probe_ids),
    n_samples     = ns,
    sample_ids    = if (!is.null(ph)) colnames(ph) else NULL,
    probe_ids     = s1$probe_ids,
    phist         = ph,
    pbins         = .MQC_PBINS,
    pfail         = pf,
    pgrid         = s1$pgrid %||% .MQC_PGRID,
    stats         = s1$stats,
    pca           = pca,
    mds           = mds,
    density       = density)
}

## Rebuild the per-sample p-value histogram from a detection matrix.
#' @keywords internal
#' @noRd
.phist_from_detp <- function(detp) {
  if (is.null(detp) || !ncol(detp)) return(NULL)
  out <- vapply(seq_len(ncol(detp)),
                function(j) phist_one(detp[, j]), integer(.MQC_PBINS))
  colnames(out) <- colnames(detp)
  out
}

## Rebuild the per-probe failure counts at each cached threshold.
#' @keywords internal
#' @noRd
.pfail_from_detp <- function(detp, grid) {
  if (is.null(detp) || !ncol(detp)) return(NULL)
  out <- matrix(0L, nrow(detp), length(grid),
                dimnames = list(rownames(detp), sprintf("p%g", grid)))
  for (g in seq_along(grid))
    out[, g] <- as.integer(rowSums(is.na(detp) | detp > grid[g]))
  out
}

#' Load the QC cache from an output directory
#'
#' @param dir output directory.
#' @param logger optional logger.
#' @return the cache list, or \code{NULL} if it is absent or unusable.
#' @export
#' @examples
#' \dontrun{
#' cc <- loadcache("~/qcout")
#' }
loadcache <- function(dir, logger = NULL) {
  logger <- logger %||% nulllog()
  p <- mqcpath(dir, "cache")
  if (!file.exists(p)) {
    logger$log("cache", "no QC cache found; falling back to the matrices")
    return(NULL)
  }
  cc <- tryCatch(readRDS(p), error = function(e) NULL)
  if (is.null(cc) || !identical(cc$cache_version, .MQC_CACHE_VERSION)) {
    logger$log("cache",
               "QC cache is missing or from a different package version; ignoring",
               warn = TRUE)
    return(NULL)
  }
  cc
}

#' Is a detection threshold available from the cached grid?
#'
#' @param cc cache list
#' @param pthresh requested threshold
#' @return \code{TRUE} when probe-level failure rates can be served from cache
#' @keywords internal
#' @noRd
cache_has_threshold <- function(cc, pthresh) {
  !is.null(cc) && !is.null(cc$pgrid) &&
    any(abs(cc$pgrid - pthresh) < 1e-12)
}

#' Summarise what a threshold change would cost
#'
#' @param cc cache list, possibly \code{NULL}
#' @param pthresh requested detection threshold
#' @return a list with \code{tier}, \code{needs_matrix} and \code{note}
#' @keywords internal
#' @noRd
plan_replot <- function(cc, pthresh) {
  if (is.null(cc))
    return(list(tier = "matrix", needs_matrix = TRUE,
                note = "no cache; reading the matrices"))
  if (cache_has_threshold(cc, pthresh))
    return(list(tier = "cache", needs_matrix = FALSE,
                note = sprintf("detp %g is on the cached grid", pthresh)))
  list(tier = "matrix", needs_matrix = TRUE,
       note = sprintf(paste("detp %g is off the cached grid (%s); the",
                            "detection matrix must be re-read"),
                      pthresh, paste(cc$pgrid, collapse = ", ")))
}
