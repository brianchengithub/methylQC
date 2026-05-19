###############################################################################
# preprocess.R — Streaming openSesame preprocessing
#
# Processes IDAT files through SeSAMe's openSesame pipeline. The prep
# code is "QCDPB" by default:
#   Q — quality masking (design issues, cross-hybridization)
#   C — channel inference
#   D — dye-bias correction
#   P — pOOBAH detection p-values
#   B — noob background correction
#
# Critically, pOOBAH (P) runs BEFORE noob (B). Detection p-values are
# therefore computed on non-noob-corrected signal (noob modifies the
# out-of-band signal that pOOBAH relies on), while the output beta and
# M-value matrices ARE noob-corrected.
#
# Outputs per cohort:
#   betas  — unmasked, noob-corrected beta values (no NAs from masking)
#   mvals  — unmasked, noob-corrected M-values
#   mask   — logical matrix (TRUE = probe failed quality or detection)
#   detP   — detection p-values (pOOBAH), computed pre-noob
#   sdfs   — list of final SigDF objects (if keep_sdfs = TRUE)
###############################################################################

#' Run openSesame and extract unmasked betas, mask, and detection p-values
#'
#' @param basenames Character vector of IDAT basenames.
#' @param platform Platform string (e.g., "EPIC", "EPICv2", "HM450").
#' @param prep_code SeSAMe prep code (default from options; "QCDPB").
#' @param n_cores Number of cores (default 1 for streaming).
#' @param collapse_to_pfx Collapse EPICv2 replicates (default FALSE).
#' @param collapse_method Collapse method: "mean" or "minPvalue".
#' @param keep_sdfs If TRUE, retains full SigDF list (memory-expensive).
#' @param logger Optional logger.
#' @return A list with betas, mvals, mask, detP, and sdfs.
#' @export
run_opensesame <- function(basenames, platform,
                           prep_code = NULL,
                           n_cores = NULL,
                           collapse_to_pfx = NULL,
                           collapse_method = NULL,
                           keep_sdfs = FALSE,
                           logger = NULL) {
  cfg <- methylQC_options()
  prep_code       <- prep_code       %||% cfg$prep_code
  n_cores         <- n_cores         %||% cfg$n_cores
  collapse_to_pfx <- collapse_to_pfx %||% cfg$collapse_to_pfx
  collapse_method <- collapse_method %||% cfg$collapse_method

  # The prep code must place P (pOOBAH) before B (noob) so detection
  # p-values are computed on non-noob signal. We additionally need a
  # pre-noob prep ("everything up to and including P") to extract the
  # detection p-values themselves.
  prep_preP <- sub("B.*$", "", prep_code)   # e.g. "QCDPB" -> "QCDP"
  has_noob  <- grepl("B", prep_code)

  n_samples <- length(basenames)
  sample_ids <- basename(basenames)

  if (!is.null(logger)) {
    logger$log("opensesame",
               sprintf("streaming %d samples (prep=%s, noob=%s, keep_sdfs=%s)",
                       n_samples, prep_code, has_noob, keep_sdfs))
  }

  get_b <- function(sdf) {
    if (collapse_to_pfx) {
      sesame::getBetas(sdf, collapseToPfx = TRUE,
                       collapseMethod = collapse_method, mask = FALSE)
    } else {
      sesame::getBetas(sdf, mask = FALSE)
    }
  }

  # Process one sample: betas (noob-corrected, unmasked), mask, detP
  # (pre-noob), and the final SigDF.
  process_one <- function(bn) {
    sdf_preP <- sesame::openSesame(bn, platform = platform,
                                   prep = prep_preP, func = NULL)
    sdf_final <- if (has_noob) {
      sesame::openSesame(bn, platform = platform,
                         prep = prep_code, func = NULL)
    } else {
      sdf_preP
    }
    list(betas = get_b(sdf_final),
         mask  = extract_mask(sdf_preP, NULL),
         detP  = extract_detP(sdf_preP, NULL),
         sdf   = sdf_final)
  }

  # --- First sample sets matrix dimensions ---
  first <- process_one(basenames[1])
  probe_ids <- names(first$betas)
  n_probes <- length(probe_ids)

  betas    <- matrix(NA_real_, n_probes, n_samples,
                     dimnames = list(probe_ids, sample_ids))
  mask_mat <- matrix(FALSE, n_probes, n_samples,
                     dimnames = list(probe_ids, sample_ids))
  detP_mat <- matrix(NA_real_, n_probes, n_samples,
                     dimnames = list(probe_ids, sample_ids))

  betas[, 1]    <- first$betas
  mask_mat[, 1] <- align_to_probes(first$mask, probe_ids, fill = FALSE)
  detP_mat[, 1] <- align_to_probes(first$detP, probe_ids, fill = NA_real_)

  sdfs <- if (keep_sdfs) {
    sl <- vector("list", n_samples); names(sl) <- sample_ids
    sl[[1]] <- first$sdf; sl
  } else NULL
  rm(first); gc(verbose = FALSE)

  progress_step <- max(25L, as.integer(n_samples / 20))

  for (i in seq.int(2, n_samples)) {
    one <- process_one(basenames[i])
    betas[, i]    <- align_to_probes(one$betas, probe_ids, fill = NA_real_)
    mask_mat[, i] <- align_to_probes(one$mask,  probe_ids, fill = FALSE)
    detP_mat[, i] <- align_to_probes(one$detP,  probe_ids, fill = NA_real_)
    if (keep_sdfs) sdfs[[i]] <- one$sdf
    rm(one)

    if (i %% progress_step == 0) {
      gc(verbose = FALSE)
      if (!is.null(logger)) {
        logger$log("opensesame",
                   sprintf("  processed %d / %d samples (%.0f%%)",
                           i, n_samples, 100 * i / n_samples))
      }
    }
  }
  gc(verbose = FALSE)

  mvals <- log2((betas + 1e-6) / (1 - betas + 1e-6))

  if (!is.null(logger)) {
    logger$log("opensesame",
               sprintf("complete: %d probes x %d samples; %.2f%% of values masked",
                       nrow(betas), ncol(betas),
                       100 * sum(mask_mat) / length(mask_mat)))
  }
  list(sdfs = sdfs, betas = betas, mvals = mvals,
       mask = mask_mat, detP = detP_mat)
}

#' Align a named vector to a target probe set
#' @keywords internal
#' @noRd
align_to_probes <- function(vec, probe_ids, fill = NA_real_) {
  if (is.null(probe_ids)) return(vec)
  out <- rep(fill, length(probe_ids))
  names(out) <- probe_ids
  if (!is.null(names(vec))) {
    m <- match(probe_ids, names(vec))
    ok <- !is.na(m)
    out[ok] <- vec[m[ok]]
  } else if (length(vec) == length(probe_ids)) {
    out[] <- vec
  }
  out
}

#' Extract the quality/detection mask from a SigDF
#' @keywords internal
#' @noRd
extract_mask <- function(sdf, probe_ids = NULL) {
  if ("mask" %in% colnames(sdf)) {
    mv <- as.logical(sdf$mask)
    names(mv) <- sdf$Probe_ID
  } else {
    mv <- stats::setNames(rep(FALSE, nrow(sdf)), sdf$Probe_ID)
  }
  if (is.null(probe_ids)) return(mv)
  align_to_probes(mv, probe_ids, fill = FALSE)
}

#' Extract detection p-values (pOOBAH) from a pre-noob SigDF
#' @keywords internal
#' @noRd
extract_detP <- function(sdf, probe_ids = NULL) {
  pvals <- tryCatch({
    pv <- sesame::pOOBAH(sdf, return.pval = TRUE)
    if (is.numeric(pv)) {
      pv
    } else {
      stats::setNames(rep(NA_real_, nrow(sdf)), sdf$Probe_ID)
    }
  }, error = function(e) {
    stats::setNames(rep(NA_real_, nrow(sdf)), sdf$Probe_ID)
  })
  if (is.null(names(pvals)) && length(pvals) == nrow(sdf)) {
    names(pvals) <- sdf$Probe_ID
  }
  if (is.null(probe_ids)) return(pvals)
  align_to_probes(pvals, probe_ids, fill = NA_real_)
}

#' Null-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
