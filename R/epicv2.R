###############################################################################
# epicv2.R — EPICv2 replicate de-duplication and probe ID harmonization
#
# EPICv2 targets many CpGs with multiple replicate probes (probe IDs
# sharing a cg/ch prefix but differing by suffix). For cohort analyses
# every sample must use the SAME physical probe for a given CpG, or
# probe-design differences become a source of technical variation.
#
# methylQC resolves this in two steps:
#   1. De-duplication: for each replicated cg/ch CpG, the probe with the
#      fewest cross-sample detection failures is kept (ties -> first
#      probe by ID order). This is a COHORT-level decision, unlike
#      SeSAMe collapseToPfx(method="minPvalue") which selects per-sample.
#   2. Harmonization: surviving EPICv2 probe IDs are lifted to the
#      target platform (EPIC or HM450) via SeSAMe's mLiftOver, which
#      uses conversion mappings shipped in the sesameData cache.
#
# rs (SNP) and control probes are left untouched; only cg/ch probes are
# de-duplicated.
###############################################################################

#' De-duplicate EPICv2 replicate probes and harmonize probe IDs
#'
#' Cohort-consistent replicate resolution followed by probe ID
#' harmonization to a legacy Infinium platform.
#'
#' @param mat Numeric matrix (probes x samples): a beta-value or
#'   M-value matrix with EPICv2 probe IDs as row names. May be supplied
#'   directly or via \code{mat_path}.
#' @param detP Numeric matrix (probes x samples) of detection p-values,
#'   same row/column names as \code{mat}. May be supplied via
#'   \code{detP_path}. Required: the de-duplication criterion is the
#'   cross-sample detection failure rate.
#' @param mat_path Optional path to an .rds file holding \code{mat}.
#' @param detP_path Optional path to an .rds file holding \code{detP}.
#' @param target_platform "EPIC" (default) or "HM450".
#' @param detP_thresh Detection p-value threshold defining a "failure"
#'   (default from options, 0.05).
#' @param dedup_log_path Optional path to write the per-CpG
#'   kept/dropped log as CSV.
#' @param return_details If TRUE, return a list with the harmonized
#'   matrix plus the probe-ID correspondence needed to re-align
#'   companion matrices (used internally by the pipeline). If FALSE
#'   (default), return just the harmonized matrix.
#' @param logger Optional logger.
#' @return If \code{return_details = FALSE}, the de-duplicated,
#'   ID-harmonized matrix (probes x samples). If TRUE, a list with
#'   \code{mat}, \code{kept_ids}, \code{id_map}, and
#'   \code{target_platform}.
#' @examples
#' \dontrun{
#' # Matrices produced by methylQC Stage 1 (betas_all.rds is unmasked,
#' # detP_all.rds holds pOOBAH detection p-values; both are
#' # probes x samples with probe IDs as row names).
#' harmonized <- harmonize_epicv2(
#'   mat_path  = "results/betas_all.rds",
#'   detP_path = "results/detP_all.rds",
#'   target_platform = "EPIC",
#'   dedup_log_path  = "results/epicv2_dedup_log.csv")
#' }
#' @export
harmonize_epicv2 <- function(mat = NULL, detP = NULL,
                             mat_path = NULL, detP_path = NULL,
                             target_platform = c("EPIC", "HM450"),
                             detP_thresh = NULL,
                             dedup_log_path = NULL,
                             return_details = FALSE,
                             logger = NULL) {
  cfg <- methylQC_options()
  target_platform <- match.arg(target_platform)
  detP_thresh <- detP_thresh %||% cfg$detP_thresh

  if (is.null(mat))  mat  <- readRDS(mat_path)
  if (is.null(detP)) detP <- readRDS(detP_path)
  if (is.null(mat) || is.null(detP))
    stop("Both 'mat' and 'detP' (or their *_path arguments) are required.")
  if (!is.matrix(mat))  mat  <- as.matrix(mat)
  if (!is.matrix(detP)) detP <- as.matrix(detP)

  # Align detP to mat by probe and sample
  common_p <- intersect(rownames(mat), rownames(detP))
  common_s <- intersect(colnames(mat), colnames(detP))
  if (length(common_p) == 0 || length(common_s) == 0)
    stop("'mat' and 'detP' share no probe IDs / sample IDs.")
  mat  <- mat[common_p, common_s, drop = FALSE]
  detP <- detP[common_p, common_s, drop = FALSE]

  loginfo <- function(...) if (!is.null(logger)) logger$log("epicv2", sprintf(...))

  # --- Step 1: cohort-consistent de-duplication (cg/ch probes only) ---
  pid <- rownames(mat)
  is_cgch <- grepl("^(cg|ch)", pid)
  # The CpG-level key strips the EPICv2 replicate suffix (everything
  # from the first underscore onward).
  cpg_key <- sub("_.*$", "", pid)

  # Per-probe cross-sample detection failure rate
  fail_rate <- rowMeans(detP > detP_thresh, na.rm = TRUE)
  fail_rate[is.nan(fail_rate)] <- 1  # all-NA probe -> treat as worst

  keep <- rep(TRUE, length(pid))
  log_rows <- list()

  # Only CpGs that are actually replicated need resolving
  cgch_keys <- cpg_key[is_cgch]
  dup_keys <- unique(cgch_keys[duplicated(cgch_keys)])

  for (k in dup_keys) {
    idx <- which(is_cgch & cpg_key == k)
    fr <- fail_rate[idx]
    best_local <- which(fr == min(fr))
    chosen <- idx[best_local[1]]      # ties -> first probe by ID order
    dropped <- setdiff(idx, chosen)
    keep[dropped] <- FALSE
    log_rows[[k]] <- data.frame(
      cpg = k,
      kept_probe = pid[chosen],
      kept_fail_rate = round(fail_rate[chosen], 5),
      dropped_probes = paste(pid[dropped], collapse = ";"),
      n_replicates = length(idx),
      stringsAsFactors = FALSE)
  }

  n_dup_cpgs <- length(dup_keys)
  n_dropped  <- sum(!keep)
  loginfo("de-duplication: %d replicated CpG(s), %d probe(s) dropped",
          n_dup_cpgs, n_dropped)

  if (!is.null(dedup_log_path) && length(log_rows) > 0) {
    dedup_log <- do.call(rbind, log_rows)
    utils::write.csv(dedup_log, dedup_log_path, row.names = FALSE)
    loginfo("wrote de-duplication log: %s", dedup_log_path)
  } else if (!is.null(dedup_log_path)) {
    utils::write.csv(
      data.frame(cpg = character(0), kept_probe = character(0),
                 kept_fail_rate = numeric(0), dropped_probes = character(0),
                 n_replicates = integer(0)),
      dedup_log_path, row.names = FALSE)
  }

  mat_dedup <- mat[keep, , drop = FALSE]
  kept_ids  <- rownames(mat_dedup)

  # --- Step 2: harmonize probe IDs to the target platform ---
  # mLiftOver on a matrix lifts probe IDs to the target platform.
  # Because replicates have already been removed, this step only
  # renames IDs (there is nothing left to collapse).
  loginfo("harmonizing %d probe(s) to %s via mLiftOver",
          nrow(mat_dedup), target_platform)
  harmonized <- tryCatch({
    sesame::mLiftOver(mat_dedup, target_platform, impute = FALSE)
  }, error = function(e) {
    stop(sprintf(paste0(
      "mLiftOver failed (%s).\n  Ensure the SeSAMe data cache is ",
      "initialized: sesameData::sesameDataCacheAll()"),
      conditionMessage(e)), call. = FALSE)
  })
  if (is.null(dim(harmonized)) && is.numeric(harmonized)) {
    harmonized <- matrix(harmonized, ncol = 1,
                         dimnames = list(names(harmonized),
                                         colnames(mat_dedup)))
  }
  loginfo("harmonization complete: %d probe(s) on %s",
          nrow(harmonized), target_platform)

  if (!return_details) return(harmonized)

  # Build the EPICv2-ID -> harmonized-ID correspondence by lifting the
  # kept probe IDs alone. mLiftOver on a character vector returns the
  # ID mapping. This lets callers re-align companion matrices (mask,
  # detection p-values) with exactly the same de-dup + rename.
  id_map <- tryCatch({
    lifted <- sesame::mLiftOver(kept_ids, target_platform)
    if (!is.null(names(lifted))) {
      stats::setNames(as.character(lifted), names(lifted))
    } else {
      stats::setNames(as.character(lifted), kept_ids)
    }
  }, error = function(e) stats::setNames(rownames(harmonized), kept_ids))

  list(mat = harmonized, kept_ids = kept_ids, id_map = id_map,
       target_platform = target_platform)
}

#' Apply an existing EPICv2 de-dup + harmonization to a companion matrix
#'
#' Given the \code{kept_ids} and \code{id_map} returned by
#' \code{harmonize_epicv2(..., return_details = TRUE)}, re-align a
#' companion matrix (e.g., the quality mask or detection p-value
#' matrix) to the harmonized probe set so all matrices stay consistent.
#'
#' @param companion Numeric/logical matrix (EPICv2 probes x samples).
#' @param kept_ids Character vector of surviving EPICv2 probe IDs.
#' @param id_map Named character vector: EPICv2 ID -> harmonized ID.
#' @param fill Value for probes absent from \code{companion}.
#' @return The companion matrix de-duplicated and ID-harmonized.
#' @keywords internal
#' @noRd
apply_epicv2_map <- function(companion, kept_ids, id_map, fill = NA) {
  present <- intersect(kept_ids, rownames(companion))
  sub <- matrix(fill, length(kept_ids), ncol(companion),
                dimnames = list(kept_ids, colnames(companion)))
  if (length(present)) sub[present, ] <- companion[present, , drop = FALSE]
  rownames(sub) <- id_map[kept_ids]
  # Drop probes that did not map to the target platform
  sub[!is.na(rownames(sub)) & rownames(sub) != "", , drop = FALSE]
}
