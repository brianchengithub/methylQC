###############################################################################
# masking.R — User-facing cleanmat() for QC application + filtering + imputation
#
# cleanmat() is the single primitive for turning the raw matrices produced by
# Stage 1 into an analysis-ready matrix. It does four things, in order:
#
#   1. Applies the per-cell quality + detection masks (NA fill).
#   2. Drops samples (columns) the user wants out.
#   3. Drops probes (rows) the user wants out.
#   4. Keeps only the requested probe categories (cg / ch / sex / snp / other).
#   5. Optionally fills the remaining NAs via chunked k-NN imputation.
#
# Sample and probe exclusions accept either an in-memory character vector or
# a path to a CSV (e.g. the failed_probes.csv and flagged_samples.csv files
# the rest of methylQC writes). The disambiguation rule is: if the argument
# is a length-1 character that satisfies file.exists(), it is treated as a
# path; otherwise it is treated as an ID vector.
#
# methylQC never assembles these lists for the user. The user passes IDs or
# file paths explicitly.
###############################################################################

#' Apply QC decisions to a matrix
#'
#' Turns a raw beta or M-value matrix from Stage 1 into an analysis-ready
#' matrix by (in order): applying the quality + detection masks; dropping
#' user-specified samples; dropping user-specified probes; keeping only
#' the requested probe categories; and optionally imputing the remaining
#' NAs by chunked k-NN.
#'
#' Probe-level masking is driven by the \code{mask} and \code{detP}
#' matrices, not by probe-ID lists. Probe SELECTION beyond masking is a
#' category choice via \code{probes}: the default keeps autosomal CpG
#' (\code{"cg"}) and non-CpG (\code{"ch"}) probes; sex-chromosome, SNP,
#' and other probes are opt-in. Drop lists (\code{dropprobes},
#' \code{dropsamples}) are user-supplied — methylQC will not pick them.
#'
#' @param mat Numeric matrix (probes x samples).
#' @param mask Logical matrix. TRUE = failed quality mask or detection.
#'   Matching values are set to NA. NULL to skip.
#' @param detP Numeric matrix of detection p-values. Values with
#'   \code{detP > pthresh} are set to NA. NULL to skip.
#' @param pthresh Detection p-value cutoff. Default: \code{mqcopts()$detp}
#'   (0.05). A probe value is kept where \code{detP <= pthresh}.
#' @param dropprobes Probes to drop (rows). One of:
#'   \itemize{
#'     \item \code{NULL}: drop no probes by ID.
#'     \item Character vector of probe IDs.
#'     \item Path to a CSV file. The function looks for a \code{probe_id}
#'           column; if absent, uses the first column. This is the form
#'           that consumes methylQC's own \code{failed_probes.csv}.
#'   }
#' @param dropsamples Samples to drop (columns). Same NULL/vector/path
#'   convention as \code{dropprobes}, with \code{sample_id} as the
#'   preferred column. Consumes \code{flagged_samples.csv}.
#' @param probes Character vector of probe categories to KEEP. Any
#'   combination of \code{"cg"} (autosomal CpG), \code{"ch"} (non-CpG),
#'   \code{"sex"} (chrX/chrY probes), \code{"snp"} (rs probes), and
#'   \code{"other"}. Default \code{c("cg", "ch")}. \code{NULL} keeps every
#'   probe (no category filtering).
#' @param platform Platform string (e.g. \code{"EPIC"}). Required when
#'   \code{probes} contains \code{"cg"} or \code{"sex"}, because sex-
#'   chromosome probes must be identified from the manifest.
#' @param impute Logical; if TRUE, impute remaining NAs via chunked k-NN
#'   after all masking/filtering. Default FALSE.
#' @param knnk Number of nearest neighbours (default 10).
#' @param chunk Probes per imputation block (default 50000).
#' @return Filtered numeric matrix.
#' @examples
#' \dontrun{
#' betas <- readRDS("results/betas_all.rds")
#' mask  <- readRDS("results/mask_all.rds")
#' detP  <- readRDS("results/detP_all.rds")
#'
#' # Standard EWAS: mask + detP, autosomal CpG + non-CpG, no exclusions
#' betas_ewas <- cleanmat(betas, mask = mask, detP = detP,
#'                        platform = "EPIC")
#'
#' # Consume the QC artifacts directly: drop probes from failed_probes.csv
#' # AND drop samples from flagged_samples.csv (paths, not vectors)
#' betas_clean <- cleanmat(
#'   betas, mask = mask, detP = detP,
#'   dropprobes  = "results/failed_probes.csv",
#'   dropsamples = "results/flagged_samples.csv",
#'   probes      = "cg", platform = "EPIC",
#'   impute      = TRUE, knnk = 10)
#'
#' # Mix forms: probe IDs from a CSV, sample IDs from an in-memory vector
#' ss <- read.csv("results/sample_sheet.csv")
#' betas_clean <- cleanmat(
#'   betas, mask = mask, detP = detP,
#'   dropprobes  = "results/failed_probes.csv",
#'   dropsamples = ss$sample_id[ss$flagged],
#'   platform = "EPIC")
#'
#' # SNP probes only (identity work) — no platform needed
#' betas_snp <- cleanmat(betas, probes = "snp")
#'
#' # Minimal: just the quality mask, keep every probe
#' betas_masked <- cleanmat(betas, mask = mask, probes = NULL)
#' }
#' @export
cleanmat <- function(mat,
                     mask = NULL,
                     detP = NULL,
                     pthresh = NULL,
                     dropprobes = NULL,
                     dropsamples = NULL,
                     probes = c("cg", "ch"),
                     platform = NULL,
                     impute = FALSE,
                     knnk = 10L,
                     chunk = 50000L) {
  stopifnot(is.matrix(mat) && is.numeric(mat))
  cfg <- mqcopts()
  if (is.null(pthresh)) pthresh <- cfg$detp

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

  # Step 3: Drop samples (columns) — vector or CSV path
  drop_s <- resolve_id_list(dropsamples, "sample_id", "Sample")
  if (length(drop_s) > 0) {
    n_before <- ncol(mat)
    keep_s   <- setdiff(colnames(mat), drop_s)
    n_drop   <- n_before - length(keep_s)
    if (length(keep_s) == 0) stop("No samples remaining after dropsamples.")
    mat <- mat[, keep_s, drop = FALSE]
    message(sprintf("  cleanmat: dropped %d sample(s); %d remain.",
                    n_drop, ncol(mat)))
  }

  # Step 4: Drop probes (rows) — vector or CSV path
  drop_p <- resolve_id_list(dropprobes, "probe_id", "Probe")
  if (length(drop_p) > 0) {
    n_before <- nrow(mat)
    keep_p   <- setdiff(rownames(mat), drop_p)
    n_drop   <- n_before - length(keep_p)
    if (length(keep_p) == 0) stop("No probes remaining after dropprobes.")
    mat <- mat[keep_p, , drop = FALSE]
    message(sprintf("  cleanmat: dropped %d probe(s); %d remain.",
                    n_drop, nrow(mat)))
  }

  # Step 5: Probe category filter (rows)
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

  # Step 6: Optional k-NN imputation
  if (impute) {
    n_na <- sum(is.na(mat))
    if (n_na > 0) mat <- impute_knn_chunked(mat, k = knnk, chunk_size = chunk)
  }

  mat
}

#' Resolve an ID exclusion argument to a character vector
#'
#' Accepts NULL, a character vector of IDs, or a length-1 file path. If
#' the argument is length 1 and \code{file.exists()} returns TRUE, the
#' file is read as CSV and the function looks for \code{id_col}; if that
#' column is absent, falls back to the first column.
#'
#' @param x NULL, character vector, or single-string file path.
#' @param id_col Preferred ID column name in the CSV.
#' @param what Label for error messages: \code{"Sample"} or \code{"Probe"}.
#' @return Character vector of unique non-empty IDs (empty if x is NULL).
#' @keywords internal
#' @noRd
resolve_id_list <- function(x, id_col, what) {
  if (is.null(x) || length(x) == 0) return(character(0))
  if (!is.character(x)) {
    stop(sprintf("%s exclusion must be NULL, character vector, or file path; got %s",
                 what, class(x)[1]), call. = FALSE)
  }
  if (length(x) == 1 && file.exists(x)) {
    df <- tryCatch(
      utils::read.csv(x, stringsAsFactors = FALSE,
                      check.names = FALSE, na.strings = c("", "NA")),
      error = function(e) {
        stop(sprintf("%s exclusion: failed to read CSV at '%s': %s",
                     what, x, conditionMessage(e)), call. = FALSE)
      })
    if (ncol(df) == 0) {
      stop(sprintf("%s exclusion CSV '%s' has no columns.", what, x),
           call. = FALSE)
    }
    ids <- if (id_col %in% colnames(df)) df[[id_col]] else df[[1]]
    ids <- as.character(ids)
    ids <- unique(ids[!is.na(ids) & nzchar(ids)])
    message(sprintf("  cleanmat: read %d %s ID(s) from '%s'.",
                    length(ids), tolower(what), x))
    return(ids)
  }
  # Length 1 with no file at that path: treat as a single ID
  unique(x[!is.na(x) & nzchar(x)])
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
