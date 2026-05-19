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
#' @param impute Logical; if TRUE, impute remaining NAs via
#'   similarity-based k-NN after all masking/filtering (default FALSE).
#' @param knn_k Number of nearest-neighbour samples (default from
#'   options, 50).
#' @param knn_var_probes Number of most-variable probes used to compute
#'   sample-to-sample distances (default from options, 30000).
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
                       knn_k = NULL,
                       knn_var_probes = NULL) {
  stopifnot(is.matrix(mat) && is.numeric(mat))
  cfg <- methylQC_options()
  knn_k          <- knn_k          %||% cfg$knn_k
  knn_var_probes <- knn_var_probes %||% cfg$knn_var_probes

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

  # Step 7: Optional similarity-based k-NN imputation
  if (impute) {
    n_na <- sum(is.na(mat))
    if (n_na > 0) {
      mat <- impute_knn_similarity(mat, k = knn_k,
                                   n_var_probes = knn_var_probes)
    }
  }

  mat
}

#' Similarity-based k-NN imputation
#'
#' Imputes missing values from the most similar samples. Sample-to-
#' sample Euclidean distances are computed on the top
#' \code{n_var_probes} most-variable probes; restricting to variable
#' probes is what makes "nearest" biologically meaningful, since the
#' bulk of the array is near-constant and would otherwise dominate the
#' distance. Each missing value is filled with the inverse-distance-
#' weighted mean of the \code{k} nearest-neighbour samples' values at
#' that probe.
#'
#' There is no constant fallback: a value missing in a sample and in
#' all of its usable neighbours is left as NA.
#'
#' @param mat Numeric matrix (probes x samples).
#' @param k Number of nearest-neighbour samples (default 50).
#' @param n_var_probes Number of most-variable probes used for the
#'   distance computation (default 30000).
#' @return The matrix with imputable NAs filled.
#' @keywords internal
#' @noRd
impute_knn_similarity <- function(mat, k = 50L, n_var_probes = 30000L) {
  n_samples <- ncol(mat)
  if (n_samples < 2L) {
    message("  Imputation skipped: need >= 2 samples.")
    return(mat)
  }
  k <- min(as.integer(k), n_samples - 1L)

  # --- Select most-variable probes for the distance computation ---
  pvar <- matrixStats::rowVars(mat, na.rm = TRUE)
  pvar[is.na(pvar)] <- -Inf
  n_var <- min(as.integer(n_var_probes), sum(is.finite(pvar)))
  var_idx <- order(pvar, decreasing = TRUE)[seq_len(n_var)]
  anchor <- mat[var_idx, , drop = FALSE]
  message(sprintf(
    "  Imputing NAs: similarity k-NN (k=%d) on %d most-variable probes...",
    k, n_var))

  # --- Sample-to-sample Euclidean distance ---
  # Computed over pairwise-complete probes and rescaled to a full-length
  # Euclidean distance, so partial probe overlap between two samples is
  # handled gracefully.
  dmat <- matrix(NA_real_, n_samples, n_samples,
                 dimnames = list(colnames(mat), colnames(mat)))
  diag(dmat) <- 0
  for (a in seq_len(n_samples - 1L)) {
    va <- anchor[, a]
    for (b in seq.int(a + 1L, n_samples)) {
      vb <- anchor[, b]
      ok <- !is.na(va) & !is.na(vb)
      d <- if (sum(ok) >= 10L) {
        sqrt(mean((va[ok] - vb[ok])^2) * nrow(anchor))
      } else NA_real_
      dmat[a, b] <- d
      dmat[b, a] <- d
    }
  }

  # --- Impute each missing value from nearest-neighbour samples ---
  out <- mat
  eps <- 1e-6
  n_filled <- 0L
  for (s in seq_len(n_samples)) {
    na_rows <- which(is.na(mat[, s]))
    if (length(na_rows) == 0L) next

    # Rank the other samples by distance to sample s
    d_s <- dmat[, s]
    d_s[s] <- NA_real_
    ord <- order(d_s, na.last = NA)            # drops NA-distance samples
    if (length(ord) == 0L) next
    nn <- utils::head(ord, k)
    w  <- 1 / (d_s[nn] + eps)

    nn_block <- mat[na_rows, nn, drop = FALSE]  # probes x neighbours
    wmat <- matrix(w, nrow = length(na_rows), ncol = length(nn),
                   byrow = TRUE)
    wmat[is.na(nn_block)] <- 0                  # ignore NA neighbours
    wsum <- rowSums(wmat)
    num  <- rowSums(wmat * ifelse(is.na(nn_block), 0, nn_block))
    imp  <- ifelse(wsum > 0, num / wsum, NA_real_)

    fill <- !is.na(imp)
    out[na_rows[fill], s] <- imp[fill]
    n_filled <- n_filled + sum(fill)
  }

  message(sprintf("  Imputation complete. Values filled: %d; remaining NAs: %d",
                  n_filled, sum(is.na(out))))
  out
}

