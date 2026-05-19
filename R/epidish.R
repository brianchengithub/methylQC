###############################################################################
# epidish.R — EpiDISH cell-type deconvolution
#
# Deconvolution is run as a purity / contamination check on blood and
# isolated blood-cell populations. The default reference is the 12
# immune cell-type panel (cent12CT.m): for a cleanly sorted population
# a single cell type should dominate; for whole blood a plausible
# leukocyte mixture is expected; high immune signal in a non-blood
# tissue would indicate contamination.
#
# EpiDISH (RPC) constrains estimates to sum to 1, so the panel is only
# applied to blood-derived tissues. Running an immune-only reference on
# a non-immune tissue would force a meaningless simplex.
###############################################################################

#' Check whether a cell-type label is a blood-derived tissue
#' @keywords internal
#' @noRd
is_blood_tissue <- function(cell_type, allowlist) {
  if (is.null(cell_type) || is.na(cell_type)) return(FALSE)
  normalize <- function(x) tolower(trimws(gsub("[_-]", " ", x)))
  any(normalize(cell_type) == normalize(allowlist))
}

#' Resolve an EpiDISH reference panel by name
#'
#' Loads a named reference centroid matrix from the EpiDISH package.
#' Stops with an actionable message if the panel is not found.
#'
#' @param reference Reference object name (e.g., "cent12CT.m").
#' @return A numeric matrix (CpGs x cell types).
#' @keywords internal
#' @noRd
resolve_epidish_reference <- function(reference) {
  ref_env <- new.env()
  ok <- tryCatch({
    utils::data(list = reference, package = "EpiDISH", envir = ref_env)
    TRUE
  }, warning = function(w) FALSE, error = function(e) FALSE)
  if (!ok || !exists(reference, envir = ref_env)) {
    stop(sprintf(paste0(
      "EpiDISH reference '%s' not found in your installed EpiDISH ",
      "package.\n  Installed EpiDISH version: %s\n",
      "  Run methylQC::check_dependencies() and update EpiDISH if needed:\n",
      "    BiocManager::install(\"EpiDISH\")"),
      reference,
      as.character(utils::packageVersion("EpiDISH"))), call. = FALSE)
  }
  ref <- get(reference, envir = ref_env)
  as.matrix(ref)
}

#' Run EpiDISH deconvolution on a beta matrix
#'
#' @param betas Numeric beta matrix (probes x samples), unmasked.
#' @param reference Reference panel name, or a character vector of
#'   names to run multiple panels. Default from options
#'   (\code{cent12CT.m}).
#' @param method EpiDISH method ("RPC", "CBS", or "CP").
#' @param logger Optional logger.
#' @return If a single reference: a numeric matrix (samples x cell
#'   types). If multiple references: a named list of such matrices.
#' @export
run_epidish <- function(betas, reference = NULL, method = NULL,
                        logger = NULL) {
  cfg <- methylQC_options()
  reference <- reference %||% cfg$epidish_reference
  method    <- method    %||% cfg$epidish_method

  # Multiple references -> named list of results
  if (length(reference) > 1) {
    res_list <- lapply(reference, function(r)
      run_epidish(betas, reference = r, method = method, logger = logger))
    names(res_list) <- reference
    return(res_list)
  }

  ref_data <- resolve_epidish_reference(reference)

  overlap <- intersect(rownames(betas), rownames(ref_data))
  if (length(overlap) < 10) {
    if (!is.null(logger)) {
      logger$log("epidish",
                 sprintf("reference '%s': only %d probes overlap - skipping",
                         reference, length(overlap)))
    }
    return(NULL)
  }
  if (!is.null(logger)) {
    logger$log("epidish",
               sprintf("reference '%s': %d / %d probes present (%.1f%%)",
                       reference, length(overlap), nrow(ref_data),
                       100 * length(overlap) / nrow(ref_data)))
  }

  res <- EpiDISH::epidish(beta.m = betas[overlap, , drop = FALSE],
                          ref.m = ref_data[overlap, , drop = FALSE],
                          method = method)$estF

  rs <- rowSums(res)
  if (!is.null(logger)) {
    logger$log("epidish",
               sprintf("row-sum check: median=%.3f range=[%.3f, %.3f] (~1 expected)",
                       stats::median(rs), min(rs), max(rs)))
    for (ct in colnames(res)) {
      qs <- stats::quantile(res[, ct], c(0.25, 0.5, 0.75), na.rm = TRUE)
      logger$log("epidish",
                 sprintf("  %-12s median=%.3f IQR=[%.3f, %.3f]",
                         ct, qs[2], qs[1], qs[3]))
    }
  }
  rm(ref_data); gc(verbose = FALSE)
  res
}
