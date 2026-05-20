###############################################################################
# epidish.R — EpiDISH cell-type deconvolution
#
# Estimates cell-type proportions for blood-derived tissues using the
# EpiDISH package with the Robust Partial Correlation (RPC) method and
# the centDHSbloodDMC.m reference panel. Only triggered when the cell
# type matches the configurable blood-tissue allowlist.
#
# Cell types in the reference:
#   B cells, NK cells, CD4T, CD8T, Monocytes, Neutrophils, Eosinophils
#
# The reference panel is loaded into a local environment to avoid
# polluting the global namespace.
###############################################################################
#' Check if a cell-type label represents a blood-like tissue
#' @keywords internal
#' @noRd
is_blood_tissue <- function(cell_type, allowlist) {
  if (is.null(cell_type) || is.na(cell_type)) return(FALSE)
  normalize <- function(x) tolower(trimws(gsub("[_-]", " ", x)))
  any(normalize(cell_type) == normalize(allowlist))
}

#' Run EpiDISH deconvolution on a beta matrix
#'
#' @param betas Beta matrix (probes x samples).
#' @param ref EpiDISH reference dataset name. Default: \code{dishref} option.
#' @param method EpiDISH method. Default: \code{dishmethod} option ("RPC").
#' @param logger Optional logger.
#' @return Cell-type proportion matrix (samples x cell types).
#' @export
rundish <- function(betas, ref = NULL, method = NULL, logger = NULL) {
  cfg    <- mqcopts()
  ref    <- if (is.null(ref))    cfg$dishref    else ref
  method <- if (is.null(method)) cfg$dishmethod else method

  # Load reference into a local environment via data(), not get()
  ref_env <- new.env()
  utils::data(list = ref, package = "EpiDISH", envir = ref_env)
  ref_data <- get(ref, envir = ref_env)

  overlap <- intersect(rownames(betas), rownames(ref_data))
  if (!is.null(logger)) {
    logger$log("rundish",
               sprintf("reference '%s': %d / %d probes present (%.1f%%)",
                       ref, length(overlap), nrow(ref_data),
                       100 * length(overlap) / nrow(ref_data)))
  }

  res <- EpiDISH::epidish(beta.m = betas[overlap, , drop = FALSE],
                          ref.m  = ref_data[overlap, , drop = FALSE],
                          method = method)$estF

  rs <- rowSums(res)
  if (!is.null(logger)) {
    logger$log("rundish",
               sprintf("row-sum check: median=%.3f range=[%.3f, %.3f] (~1 expected)",
                       stats::median(rs), min(rs), max(rs)))
    for (ct in colnames(res)) {
      qs <- stats::quantile(res[, ct], c(0.25, 0.5, 0.75), na.rm = TRUE)
      logger$log("rundish",
                 sprintf("  %-12s median=%.3f IQR=[%.3f, %.3f]",
                         ct, qs[2], qs[1], qs[3]))
    }
  }
  rm(ref_env, ref_data); gc(verbose = FALSE)
  res
}
