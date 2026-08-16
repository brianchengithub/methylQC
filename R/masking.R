## ---------------------------------------------------------------------------
## masking.R -- applying masks, on demand, under a policy the caller chooses.
##
## v3 never stores a masked matrix. betas_all.rds holds every probe for every
## sample with no NA introduced by masking; the design mask and the detection
## p-values are stored separately. Masking is a view, applied here.
##
## Two mask types, deliberately kept apart:
##
##   design     -- probes with known design problems (multi-mapping, SNP
##                 overlap, cross-hybridising). A platform property, identical
##                 for every sample, stored as a named logical vector.
##
##   detection  -- probes whose signal was indistinguishable from background in
##                 THIS sample. Derived by thresholding the stored ELBAR
##                 p-values, so the threshold is a decision made at use time,
##                 not baked into Stage 1.
##
## v2 fused them into one logical matrix via sesame's addMask(), which only
## ever sets TRUE and never records why. plots.R then had to reverse-engineer
## the split as `mask & (detP <= pthresh)`, which silently misclassifies
## whenever the mask and the p-values disagree -- and in v2 they always did,
## because the p-values were recomputed after noob had rewritten the signal
## they were derived from.
## ---------------------------------------------------------------------------

#' Apply a masking policy to a beta matrix
#'
#' @param betas probes-by-samples numeric matrix (unmasked, as stored).
#' @param detp matching detection p-value matrix, or \code{NULL}.
#' @param design named logical design mask, or \code{NULL}.
#' @param maskuse which masks to apply: \code{"both"}, \code{"detection"},
#'   \code{"design"} or \code{"none"}.
#' @param pthresh detection p-value threshold.
#' @param logger optional logger.
#' @return \code{betas} with masked cells set to \code{NA}.
#' @export
#' @examples
#' b <- matrix(runif(6), 3, 2,
#'             dimnames = list(c("cg1", "cg2", "cg3"), c("s1", "s2")))
#' d <- matrix(c(0.001, 0.9, 0.001, 0.001, 0.001, 0.001), 3, 2,
#'             dimnames = dimnames(b))
#' cleanmat(b, detp = d, maskuse = "detection")
cleanmat <- function(betas, detp = NULL, design = NULL,
                     maskuse = NULL, pthresh = NULL, logger = NULL) {

  logger <- logger %||% nulllog()
  cfg <- mqcopts()
  maskuse <- .opt(maskuse, "maskuse", cfg)
  pthresh <- .opt(pthresh, "detp", cfg)
  maskuse <- match.arg(maskuse, c("both", "detection", "design", "none"))

  if (!is.matrix(betas)) stop("betas must be a matrix", call. = FALSE)
  if (identical(maskuse, "none")) return(betas)

  n_before <- sum(is.na(betas))
  use_det <- maskuse %in% c("both", "detection")
  use_des <- maskuse %in% c("both", "design")

  if (use_det) {
    if (is.null(detp)) {
      logger$log("mask", "maskuse requests detection masking but no detP was supplied",
                 warn = TRUE)
    } else {
      detp <- conform_to(detp, betas, "detP matrix")
      ## NA is treated as a failure, not as a pass. In v2 the index
      ## `detp > pthresh` carried NAs, and R leaves NA-indexed elements
      ## untouched, so a probe whose detection could not be computed was
      ## silently treated as having passed.
      bad <- is.na(detp) | detp > pthresh
      n_na <- sum(is.na(detp))
      betas[bad] <- NA_real_
      if (n_na > 0)
        logger$log("mask", sprintf(
          "%d detection p-value(s) were NA and were treated as failures", n_na))
    }
  }

  if (use_des) {
    if (is.null(design)) {
      logger$log("mask", "maskuse requests design masking but no design mask was supplied",
                 warn = TRUE)
    } else {
      dm <- align_to(design, rownames(betas), "design mask", fill = TRUE)
      dm[is.na(dm)] <- TRUE                     # unknown means masked
      betas[dm, ] <- NA_real_
    }
  }

  n_after <- sum(is.na(betas))
  logger$log("mask", sprintf("maskuse='%s', pthresh=%.4g: masked %s additional cell(s)",
                             maskuse, pthresh,
                             format(n_after - n_before, big.mark = ",")))
  betas
}

#' Convert beta values to M-values
#'
#' M-values are the logit of beta. A small offset keeps the result finite at
#' beta values of exactly 0 or 1, bounding the range at roughly +/- 20.
#'
#' v2 materialised a full M-value matrix in Stage 1 and wrote it to disk,
#' which at 1500 EPIC samples is 9.7 GiB of storage for a quantity that is a
#' one-line transform of a matrix already on disk. v3 provides this helper
#' instead.
#'
#' @param betas a beta value matrix or vector.
#' @param offset numeric offset guarding against 0 and 1.
#' @return M-values, same shape as the input.
#' @export
#' @examples
#' mvals(c(0, 0.5, 1))
mvals <- function(betas, offset = 1e-6) {
  log2((betas + offset) / (1 - betas + offset))
}

#' Detection mask derived from stored p-values
#'
#' @param detp detection p-value matrix
#' @param pthresh threshold
#' @return a logical matrix, \code{TRUE} where the probe failed detection.
#'   \code{NA} p-values count as failures.
#' @export
#' @examples
#' d <- matrix(c(0.01, 0.9, NA, 0.01), 2, 2)
#' detmask(d, 0.05)
detmask <- function(detp, pthresh = NULL) {
  pthresh <- .opt(pthresh, "detp")
  is.na(detp) | detp > pthresh
}

#' Per-sample call rate
#'
#' The fraction of probes on the array that were detected. The denominator is
#' every probe, so a probe whose p-value is missing counts against the sample.
#'
#' v2 used \code{colMeans(detP <= pthresh, na.rm = TRUE)}, which removes
#' missing probes from the denominator as well as the numerator, so a sample
#' with half its p-values missing scored a perfect 1.00.
#'
#' @param detp detection p-value matrix
#' @param pthresh threshold
#' @return a named numeric vector
#' @export
#' @examples
#' d <- matrix(c(0.01, 0.9, NA, 0.01), 2, 2, dimnames = list(NULL, c("a", "b")))
#' callrate(d, 0.05)
callrate <- function(detp, pthresh = NULL) {
  pthresh <- .opt(pthresh, "detp")
  colMeans(!is.na(detp) & detp <= pthresh)
}

#' Per-probe failure rate across samples
#'
#' @param detp detection p-value matrix
#' @param pthresh threshold
#' @return a named numeric vector; \code{NA} p-values count as failures
#' @export
#' @examples
#' d <- matrix(c(0.01, 0.9, NA, 0.01), 2, 2, dimnames = list(c("cg1", "cg2"), NULL))
#' probefail(d, 0.05)
probefail <- function(detp, pthresh = NULL) {
  pthresh <- .opt(pthresh, "detp")
  rowMeans(is.na(detp) | detp > pthresh)
}
