###############################################################################
# preprocess.R — Streaming openSesame preprocessing
#
# Processes IDAT files through SeSAMe's openSesame pipeline (quality
# masking, noob, dye-bias correction, pOOBAH detection p-values).
#
# Outputs three matrices per cohort:
#   betas — unmasked beta values (no NAs from masking)
#   mask  — logical matrix (TRUE = probe failed quality or detection)
#   detP  — detection p-values (pOOBAH) per probe per sample
#
# The user applies masks downstream via cleanmat().
###############################################################################

#' Run openSesame and extract unmasked betas, mask, and detection p-values
#'
#' @param basenames Character vector of IDAT basenames.
#' @param platform Platform string (e.g. \code{"EPIC"}).
#' @param cores Number of cores. Default: \code{cores} option (1, streaming).
#' @param collapse Collapse EPICv2 replicates. Default: \code{collapse} option.
#' @param collapsemethod Collapse method: \code{"mean"} or \code{"minPvalue"}.
#'   Default: \code{collapsemethod} option.
#' @param keepsdf If TRUE, retains the full SigDF list (memory-expensive).
#' @param logger Optional logger.
#' @return A list with:
#'   \describe{
#'     \item{betas}{Unmasked beta matrix (all probes, no NAs from masking)}
#'     \item{mvals}{Unmasked M-value matrix}
#'     \item{mask}{Logical matrix: TRUE = probe failed quality mask or
#'       detection p-value}
#'     \item{detP}{Detection p-value matrix (pOOBAH)}
#'     \item{sdfs}{List of SigDF objects (NULL if \code{keepsdf = FALSE})}
#'   }
#' @export
runsesame <- function(basenames, platform,
                      cores = NULL,
                      collapse = NULL,
                      collapsemethod = NULL,
                      keepsdf = FALSE,
                      logger = NULL) {
  cfg            <- mqcopts()
  cores          <- if (is.null(cores))          cfg$cores          else cores
  collapse       <- if (is.null(collapse))       cfg$collapse       else collapse
  collapsemethod <- if (is.null(collapsemethod)) cfg$collapsemethod else collapsemethod

  n_samples  <- length(basenames)
  sample_ids <- basename(basenames)

  if (!is.null(logger)) {
    logger$log("runsesame",
               sprintf("streaming %d samples (collapse=%s, keepsdf=%s)",
                       n_samples, collapse, keepsdf))
  }

  # --- Process first sample to determine matrix dimensions ---
  first_sdf <- sesame::openSesame(basenames[1], platform = platform,
                                  func = NULL)
  # Unmasked betas: mask=FALSE gives ALL betas regardless of QC status
  first_betas <- if (collapse) {
    sesame::getBetas(first_sdf, collapseToPfx = TRUE,
                     collapseMethod = collapsemethod, mask = FALSE)
  } else {
    sesame::getBetas(first_sdf, mask = FALSE)
  }
  n_probes  <- length(first_betas)
  probe_ids <- names(first_betas)

  # Pre-allocate all three matrices
  betas <- matrix(NA_real_, nrow = n_probes, ncol = n_samples,
                  dimnames = list(probe_ids, sample_ids))
  mask_mat <- matrix(FALSE, nrow = n_probes, ncol = n_samples,
                     dimnames = list(probe_ids, sample_ids))
  detP_mat <- matrix(NA_real_, nrow = n_probes, ncol = n_samples,
                     dimnames = list(probe_ids, sample_ids))

  # Fill first sample
  betas[, 1]    <- first_betas
  mask_mat[, 1] <- extract_mask(first_sdf, probe_ids)
  detP_mat[, 1] <- extract_detP(first_sdf, probe_ids)

  sdfs <- if (keepsdf) {
    sdf_list <- vector("list", n_samples)
    names(sdf_list) <- sample_ids
    sdf_list[[1]] <- first_sdf
    sdf_list
  } else {
    rm(first_sdf)
    NULL
  }
  rm(first_betas)

  # --- Stream remaining samples ---
  progress_step <- max(25L, as.integer(n_samples / 20))

  for (i in seq.int(2, n_samples)) {
    sdf <- sesame::openSesame(basenames[i], platform = platform,
                              func = NULL)
    b <- if (collapse) {
      sesame::getBetas(sdf, collapseToPfx = TRUE,
                       collapseMethod = collapsemethod, mask = FALSE)
    } else {
      sesame::getBetas(sdf, mask = FALSE)
    }
    betas[, i]    <- b
    mask_mat[, i] <- extract_mask(sdf, probe_ids)
    detP_mat[, i] <- extract_detP(sdf, probe_ids)

    if (keepsdf) sdfs[[i]] <- sdf else rm(sdf)
    rm(b)

    if (i %% progress_step == 0) {
      gc(verbose = FALSE)
      if (!is.null(logger)) {
        logger$log("runsesame",
                   sprintf("  processed %d / %d samples (%.0f%%)",
                           i, n_samples, 100 * i / n_samples))
      }
    }
  }

  gc(verbose = FALSE)

  # Derive M-values (also unmasked)
  mvals <- log2((betas + 1e-6) / (1 - betas + 1e-6))

  if (!is.null(logger)) {
    n_masked <- sum(mask_mat)
    n_total  <- length(mask_mat)
    logger$log("runsesame",
               sprintf("complete: %d probes x %d samples; %.2f%% of values masked",
                       nrow(betas), ncol(betas), 100 * n_masked / n_total))
  }
  list(sdfs = sdfs, betas = betas, mvals = mvals,
       mask = mask_mat, detP = detP_mat)
}

#' Extract the quality/detection mask from a SigDF
#'
#' Returns a logical vector (TRUE = masked) aligned to \code{probe_ids}.
#' The mask in the SigDF combines quality masking (design issues,
#' cross-hybridisation) and detection p-value failures (pOOBAH).
#'
#' @param sdf A SeSAMe SigDF object.
#' @param probe_ids Character vector of probe IDs to align to.
#' @return Logical vector of length(probe_ids).
#' @keywords internal
#' @noRd
extract_mask <- function(sdf, probe_ids) {
  if ("mask" %in% colnames(sdf)) {
    mask_vec <- as.logical(sdf$mask)
    names(mask_vec) <- sdf$Probe_ID
    return(mask_vec[probe_ids])
  }
  rep(FALSE, length(probe_ids))
}

#' Extract detection p-values from a SigDF
#'
#' Uses sesame::pOOBAH to compute out-of-band detection p-values.
#' Returns a numeric vector aligned to \code{probe_ids}.
#'
#' @param sdf A SeSAMe SigDF object.
#' @param probe_ids Character vector of probe IDs to align to.
#' @return Numeric vector of p-values, length(probe_ids).
#' @keywords internal
#' @noRd
extract_detP <- function(sdf, probe_ids) {
  pvals <- tryCatch({
    pv <- sesame::pOOBAH(sdf, return.pval = TRUE)
    if (is.numeric(pv)) {
      pv
    } else if (is.data.frame(pv) || methods::is(pv, "SigDF")) {
      rep(NA_real_, length(probe_ids))
    } else {
      rep(NA_real_, length(probe_ids))
    }
  }, error = function(e) {
    rep(NA_real_, length(probe_ids))
  })
  if (length(pvals) == length(probe_ids) && !is.null(names(pvals))) {
    return(pvals[probe_ids])
  }
  if (length(pvals) == length(probe_ids)) return(pvals)
  if (!is.null(names(pvals))) {
    out <- rep(NA_real_, length(probe_ids))
    m <- match(probe_ids, names(pvals))
    out[!is.na(m)] <- pvals[m[!is.na(m)]]
    return(out)
  }
  rep(NA_real_, length(probe_ids))
}

#' Null-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
