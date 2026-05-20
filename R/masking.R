###############################################################################
# masking.R — User-facing applymask() for flexible filtering + imputation
#
# Probe-level masking comes ENTIRELY from the mask and detP matrices. There
# is no probe-exclusion list. Probe selection beyond masking is a category
# choice via `probes` (cg / ch / sex / snp / other); the user opts in to any
# combination. Sample exclusion is hand-assembled by the user and passed as
# `exclude` — methylQC never removes samples on its own.
###############################################################################

#' Apply masks, category filters, and optional imputation
#'
#' Probe-level masking is driven solely by the \code{mask} and \code{detP}
#' matrices. Probe SELECTION is a category choice via \code{probes}: the
#' default keeps autosomal CpG (\code{"cg"}) and non-CpG (\code{"ch"})
#' probes; sex-chromosome, SNP, and other probes are opt-in. Sample removal
#' happens only for IDs the user explicitly passes to \code{exclude}.
#'
#' @param mat Numeric matrix (probes x samples).
#' @param mask Logical matrix. TRUE = failed quality mask or detection.
#'   Matching values are set to NA. NULL to skip.
#' @param detP Numeric matrix of detection p-values. Values with
#'   \code{detP > pthresh} are set to NA. NULL to skip.
#' @param pthresh Detection p-value cutoff (default 0.05). A probe value is
#'   kept where \code{detP <= pthresh}.
#' @param exclude Character vector of sample IDs (columns) to remove. The
#'   user assembles this themselves (a vector, or read from a CSV). methylQC
#'   does not generate it. NULL keeps all samples.
#' @param probes Character vector of probe categories to KEEP. Any
#'   combination of \code{"cg"} (autosomal CpG), \code{"ch"} (non-CpG),
#'   \code{"sex"} (chrX/chrY probes), \code{"snp"} (rs probes), and
#'   \code{"other"}. Default \code{c("cg", "ch")}. \code{NULL} keeps every
#'   probe (no category filtering).
#' @param platform Platform string (e.g., \code{"EPIC"}). Required when
#'   \code{probes} contains \code{"cg"} or \code{"sex"}, because sex-
#'   chromosome probes must be identified from the manifest.
#' @param impute Logical; if TRUE, impute remaining NAs via chunked k-NN
#'   after all masking/filtering (default FALSE).
#' @param knnk Number of nearest neighbours (default 10).
#' @param chunk Probes per imputation block (default 50000).
#' @return Filtered numeric matrix.
#' @examples
#' \dontrun{
#' betas <- readRDS("results/betas_all.rds")
#' mask  <- readRDS("results/mask_all.rds")
#' detP  <- readRDS("results/detP_all.rds")
#' ss    <- read.csv("results/sample_sheet.csv")
#'
#' # Standard EWAS: mask + detP, autosomal CpG + non-CpG (the default)
#' betas_ewas <- applymask(betas, mask = mask, detP = detP,
#'                         platform = "EPIC")
#'
#' # Drop samples you chose to exclude
#' betas_ewas <- applymask(betas, mask = mask, detP = detP,
#'                         exclude = ss$sample_id[ss$flagged],
#'                         platform = "EPIC")
#'
#' # Autosomal CpG only
#' betas_cg <- applymask(betas, mask = mask, detP = detP,
#'                       probes = "cg", platform = "EPIC")
#'
#' # Keep sex-chromosome probes alongside autosomal CpG
#' betas_sex <- applymask(betas, mask = mask, detP = detP,
#'                        probes = c("cg", "sex"), platform = "EPIC")
#'
#' # SNP probes only (identity work) — no platform needed
#' betas_snp <- applymask(betas, probes = "snp")
#'
#' # Minimal: just the quality mask, keep every probe
#' betas_masked <- applymask(betas, mask = mask, probes = NULL)
#'
#' # With imputation
#' betas_imp <- applymask(betas, mask = mask, detP = detP,
#'                        platform = "EPIC", impute = TRUE)
#' }
#' @export
applymask <- function(mat,
                      mask = NULL,
                      detP = NULL,
                      pthresh = 0.05,
                      exclude = NULL,
                      probes = c("cg", "ch"),
                      platform = NULL,
                      impute = FALSE,
                      knnk = 10L,
                      chunk = 50000L) {
  stopifnot(is.matrix(mat) && is.numeric(mat))

  # Step 1: Quality + detection mask
  if (!is.null(mask)) {
    stopifnot(identical(dim(mask), dim(mat)))
    mat[mask] <- NA
  }

  # Step 2: Detection p-value
  if (!is.null(detP)) {
    stopifnot(identical(dim(detP), dim(mat)))
    mat[detP > pthresh] <- NA
  }

  # Step 3: Remove user-supplied excluded samples (columns)
  if (!is.null(exclude) && length(exclude) > 0) {
    keep_s <- setdiff(colnames(mat), exclude)
    if (length(keep_s) == 0) stop("No samples remaining after exclusion.")
    mat <- mat[, keep_s, drop = FALSE]
  }

  # Step 4: Probe category filter (rows)
  if (!is.null(probes) && length(probes) > 0) {
    valid_cats <- c("cg", "ch", "sex", "snp", "other")
    bad <- setdiff(probes, valid_cats)
    if (length(bad) > 0) {
      stop("Unknown probe category: ", paste(bad, collapse = ", "),
           ". Valid: ", paste(valid_cats, collapse = ", "))
    }
    rn     <- rownames(mat)
    is_cg  <- grepl("^cg", rn)
    is_ch  <- grepl("^ch", rn)
    is_rs  <- grepl("^rs", rn)
    needs_sex <- any(c("cg", "sex") %in% probes)
    if (needs_sex) {
      if (is.null(platform)) {
        stop("'platform' is required when 'probes' contains \"cg\" or ",
             "\"sex\" (sex-chromosome probes must be identified).")
      }
      is_sex <- rn %in% sexprobes(platform)
    } else {
      is_sex <- rep(FALSE, length(rn))
    }
    is_other <- !is_cg & !is_ch & !is_rs & !is_sex
    keep_p <- rep(FALSE, length(rn))
    if ("cg"    %in% probes) keep_p <- keep_p | (is_cg & !is_sex)
    if ("ch"    %in% probes) keep_p <- keep_p | (is_ch & !is_sex)
    if ("sex"   %in% probes) keep_p <- keep_p | is_sex
    if ("snp"   %in% probes) keep_p <- keep_p | (is_rs & !is_sex)
    if ("other" %in% probes) keep_p <- keep_p | is_other
    mat <- mat[keep_p, , drop = FALSE]
  }

  # Step 5: Optional k-NN imputation
  if (impute) {
    n_na <- sum(is.na(mat))
    if (n_na > 0) mat <- impute_knn_chunked(mat, k = knnk, chunk_size = chunk)
  }

  mat
}

#' Chunked k-NN imputation (internal)
#' @keywords internal
#' @noRd
impute_knn_chunked <- function(mat, k = 10L, chunk_size = 50000L) {
  n_probes <- nrow(mat); n_chunks <- ceiling(n_probes / chunk_size)
  message(sprintf("  Imputing NAs: k-NN (k=%d) in %d chunks of ~%d probes...",
                   k, n_chunks, chunk_size))
  out <- mat
  for (i in seq_len(n_chunks)) {
    idx <- ((i - 1L) * chunk_size + 1L):min(i * chunk_size, n_probes)
    chunk <- mat[idx, , drop = FALSE]
    if (!anyNA(chunk)) next
    chunk_imp <- tryCatch({
      invisible(utils::capture.output(
        result <- impute::impute.knn(chunk, k = k, rowmax = 0.8, colmax = 0.95),
        file = nullfile()))
      result$data
    }, error = function(e) {
      message(sprintf("    chunk %d/%d: k-NN failed, using row means",
                       i, n_chunks))
      rm <- rowMeans(chunk, na.rm = TRUE)
      for (j in seq_len(ncol(chunk))) {
        na_r <- is.na(chunk[, j]); chunk[na_r, j] <- rm[na_r]
      }
      chunk
    })
    out[idx, ] <- chunk_imp; rm(chunk, chunk_imp); gc(verbose = FALSE)
  }
  still_na <- which(rowSums(is.na(out)) > 0)
  if (length(still_na) > 0) {
    message(sprintf("  %d probes still NA after k-NN; filling with 0.5",
                     length(still_na)))
    for (i in still_na) out[i, is.na(out[i, ])] <- 0.5
  }
  message(sprintf("  Imputation complete. Remaining NAs: %d", sum(is.na(out))))
  out
}
