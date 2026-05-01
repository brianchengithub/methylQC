###############################################################################
# masking.R — User-facing apply_mask() for flexible filtering + imputation
###############################################################################

#' Apply masks, filters, and optional imputation to a beta or M-value matrix
#'
#' @param mat Numeric matrix (probes x samples).
#' @param mask Logical matrix. TRUE = failed. Values set to NA. NULL to skip.
#' @param detP Numeric matrix of detection p-values. NULL to skip.
#' @param detP_thresh Detection p-value threshold (default 0.05).
#' @param exclude_probes Character vector of probe IDs to remove (rows).
#' @param exclude_samples Character vector of sample IDs to remove (columns).
#' @param exclude_sex Logical; if TRUE (default), remove sex chromosome
#'   probes using the platform manifest. Requires \code{platform} to be set.
#'   Set to FALSE to retain sex chromosome probes for sex-specific analyses.
#' @param platform Platform string (e.g., "EPIC"). Required when
#'   \code{exclude_sex = TRUE}.
#' @param probe_types Character vector of probe prefixes to KEEP.
#'   Default c("cg", "ch"). Use "cg" for CpG only. NULL keeps all.
#'   Note: sex chromosome probes are handled separately by \code{exclude_sex},
#'   not by this argument.
#' @param impute Logical; if TRUE, impute remaining NAs via chunked
#'   k-NN after all masking/filtering (default FALSE).
#' @param knn_k Number of nearest neighbors (default 10).
#' @param chunk_size Probes per imputation block (default 50000).
#' @return Filtered numeric matrix.
#' @examples
#' \dontrun{
#' betas <- readRDS("results/betas_all.rds")
#' mask  <- readRDS("results/mask_all.rds")
#' detP  <- readRDS("results/detP_all.rds")
#' excl_p <- read.csv("results/exclude_probes.csv")
#' ss     <- read.csv("results/sample_sheet.csv")
#'
#' # Standard EWAS: mask + detP + exclude flagged + autosomal CpG only
#' betas_ewas <- apply_mask(betas, mask = mask, detP = detP,
#'   exclude_probes = excl_p$probe_id,
#'   exclude_samples = ss$sample_id[ss$flagged],
#'   exclude_sex = TRUE, platform = "EPIC",
#'   probe_types = "cg")
#'
#' # Keep sex chromosomes for sex-specific analysis
#' betas_with_sex <- apply_mask(betas, mask = mask, detP = detP,
#'   exclude_probes = excl_p$probe_id[!grepl("sex_chrom", excl_p$reason)],
#'   exclude_sex = FALSE,
#'   probe_types = "cg")
#'
#' # Minimal: just quality mask, keep everything
#' betas_masked <- apply_mask(betas, mask = mask)
#'
#' # With imputation
#' betas_imputed <- apply_mask(betas, mask = mask, detP = detP,
#'   exclude_sex = TRUE, platform = "EPIC",
#'   probe_types = "cg", impute = TRUE)
#' }
#' @export
apply_mask <- function(mat,
                       mask = NULL,
                       detP = NULL,
                       detP_thresh = 0.05,
                       exclude_probes = NULL,
                       exclude_samples = NULL,
                       exclude_sex = TRUE,
                       platform = NULL,
                       probe_types = c("cg", "ch"),
                       impute = FALSE,
                       knn_k = 10L,
                       chunk_size = 50000L) {
  stopifnot(is.matrix(mat) && is.numeric(mat))

  # Step 1: Quality mask
  if (!is.null(mask)) { stopifnot(identical(dim(mask), dim(mat))); mat[mask] <- NA }

  # Step 2: Detection p-value
  if (!is.null(detP)) { stopifnot(identical(dim(detP), dim(mat))); mat[detP > detP_thresh] <- NA }

  # Step 3: Remove excluded samples
  if (!is.null(exclude_samples) && length(exclude_samples) > 0) {
    keep_s <- setdiff(colnames(mat), exclude_samples)
    if (length(keep_s) == 0) stop("No samples remaining after exclusion.")
    mat <- mat[, keep_s, drop = FALSE]
  }

  # Step 4: Remove excluded probes
  if (!is.null(exclude_probes) && length(exclude_probes) > 0) {
    keep_p <- setdiff(rownames(mat), exclude_probes)
    mat <- mat[keep_p, , drop = FALSE]
  }

  # Step 5: Remove sex chromosome probes
  if (exclude_sex) {
    if (is.null(platform)) stop("'platform' is required when exclude_sex = TRUE")
    sex_probes <- get_sex_probes(platform)
    mat <- mat[!(rownames(mat) %in% sex_probes), , drop = FALSE]
  }

  # Step 6: Filter by probe type prefix
  if (!is.null(probe_types) && length(probe_types) > 0) {
    pattern <- paste0("^(", paste(probe_types, collapse = "|"), ")")
    mat <- mat[grepl(pattern, rownames(mat)), , drop = FALSE]
  }

  # Step 7: Optional k-NN imputation
  if (impute) {
    n_na <- sum(is.na(mat))
    if (n_na > 0) mat <- impute_knn_chunked(mat, k = knn_k, chunk_size = chunk_size)
  }

  mat
}

#' Chunked k-NN imputation (internal)
#' @keywords internal
#' @noRd
impute_knn_chunked <- function(mat, k = 10L, chunk_size = 50000L) {
  n_probes <- nrow(mat); n_chunks <- ceiling(n_probes / chunk_size)
  message(sprintf("  Imputing NAs: k-NN (k=%d) in %d chunks of ~%d probes...", k, n_chunks, chunk_size))
  out <- mat
  for (i in seq_len(n_chunks)) {
    idx <- ((i - 1L) * chunk_size + 1L):min(i * chunk_size, n_probes)
    chunk <- mat[idx, , drop = FALSE]
    if (!anyNA(chunk)) next
    chunk_imp <- tryCatch({
      invisible(utils::capture.output(
        result <- impute::impute.knn(chunk, k = k, rowmax = 0.8, colmax = 0.95), file = nullfile()))
      result$data
    }, error = function(e) {
      message(sprintf("    chunk %d/%d: k-NN failed, using row means", i, n_chunks))
      rm <- rowMeans(chunk, na.rm = TRUE)
      for (j in seq_len(ncol(chunk))) { na_r <- is.na(chunk[, j]); chunk[na_r, j] <- rm[na_r] }
      chunk
    })
    out[idx, ] <- chunk_imp; rm(chunk, chunk_imp); gc(verbose = FALSE)
  }
  still_na <- which(rowSums(is.na(out)) > 0)
  if (length(still_na) > 0) {
    message(sprintf("  %d probes still NA after k-NN; filling with 0.5", length(still_na)))
    for (i in still_na) out[i, is.na(out[i, ])] <- 0.5
  }
  message(sprintf("  Imputation complete. Remaining NAs: %d", sum(is.na(out))))
  out
}
