#' Build sample exclusion table from QC metrics
#'
#' Uses SeSAMe's \code{frac_dt} as the call-rate metric. Resolves the
#' sex column via flexible alias matching, and normalizes sex encodings
#' to M/F before comparing.
#'
#' @param sample_sheet_full Sample sheet merged with QC metrics.
#' @param logger Optional logger.
#' @return A \code{data.frame} with columns sample_id, frac_dt,
#'   frac_usable, reason.
#' @export
build_sample_exclusions <- function(sample_sheet_full, logger = NULL) {
  cfg <- methylQC_options()
  cols <- colnames(sample_sheet_full)

  sex_col <- resolve_column(cols, cfg$reported_sex_col,
                            cfg$reported_sex_aliases)
  has_sex_check <- !is.na(sex_col) && "inferred_sex" %in% cols
  has_intensity <- "mean_intensity" %in% cols
  has_frac_dt   <- "frac_dt" %in% cols
  has_batch     <- "batch_folder" %in% cols

  if (!has_frac_dt) {
    stop("Sample sheet does not contain `frac_dt` column. ",
         "Did you forget to merge in compute_qc_metrics() output?")
  }

  if (has_sex_check && !is.null(logger)) {
    logger$log("excl_samples",
               sprintf("resolved sex column: '%s'", sex_col))
  }

  normalized_reported_sex <- if (has_sex_check) {
    normalize_sex(sample_sheet_full[[sex_col]], col_name = sex_col)
  } else NULL
  normalized_inferred_sex <- if (has_sex_check) {
    normalize_sex(sample_sheet_full$inferred_sex, col_name = "inferred_sex")
  } else NULL

  reasons <- vapply(seq_len(nrow(sample_sheet_full)), function(i) {
    row <- sample_sheet_full[i, , drop = FALSE]
    parts <- character(0)
    if (!is.na(row$frac_dt) && row$frac_dt < cfg$sample_call_min) {
      parts <- c(parts, sprintf("low_detection(%.3f)", row$frac_dt))
    }
    if (has_intensity && !is.na(row$mean_intensity) &&
        row$mean_intensity < cfg$intensity_min) {
      parts <- c(parts, sprintf("low_intensity(%.0f)", row$mean_intensity))
    }
    if (has_sex_check) {
      rep_sex <- normalized_reported_sex[i]
      inf_sex <- normalized_inferred_sex[i]
      if (!is.na(rep_sex) && !is.na(inf_sex) && rep_sex != inf_sex) {
        parts <- c(parts, sprintf("sex_mismatch(%s/%s)",
                                  as.character(row[[sex_col]]),
                                  as.character(row$inferred_sex)))
      }
    }
    paste(parts, collapse = ";")
  }, character(1))

  base_cols <- list(
    sample_id   = sample_sheet_full$sample_id,
    frac_dt     = sample_sheet_full$frac_dt,
    frac_usable = sample_sheet_full$frac_usable,
    reason      = reasons
  )
  if (has_batch) {
    base_cols <- c(
      list(sample_id = sample_sheet_full$sample_id,
           batch_folder = sample_sheet_full$batch_folder),
      base_cols[-1]
    )
  }
  excl <- do.call(data.frame, c(base_cols, stringsAsFactors = FALSE))
  excl <- excl[nchar(excl$reason) > 0, ]

  if (!is.null(logger)) {
    logger$log("excl_samples",
               sprintf("flagged %d / %d (%.1f%%)",
                       nrow(excl), nrow(sample_sheet_full),
                       100 * nrow(excl) / nrow(sample_sheet_full)))
  }
  excl
}

#' Build probe exclusion table
#'
#' Flags probes for exclusion: low call rate, sex chromosomes (chrX/chrY),
#' SNP/rs probes, and any other non-cg/ch probes. Only probes starting
#' with "cg" or "ch" survive into downstream matrices.
#'
#' @param betas Beta matrix.
#' @param platform Platform string.
#' @param logger Optional logger.
#' @return A data.frame of excluded probes. Has attribute "all_probe_rates".
#' @export
build_probe_exclusions <- function(betas, platform, logger = NULL) {
  cfg <- methylQC_options()

  pct_failed <- rowMeans(is.na(betas)) * 100
  call_rate  <- 1 - pct_failed / 100

  sex_probes <- get_sex_probes(platform)
  is_sex <- rownames(betas) %in% sex_probes
  is_low <- call_rate < cfg$probe_call_min
  is_snp <- grepl("^rs", rownames(betas))
  is_non_cg_ch <- !grepl("^(cg|ch)", rownames(betas))

  reason <- vapply(seq_len(nrow(betas)), function(i) {
    parts <- character(0)
    if (is_low[i])      parts <- c(parts, sprintf("low_call_rate(%.3f)", call_rate[i]))
    if (is_sex[i])      parts <- c(parts, "sex_chrom")
    if (is_snp[i])      parts <- c(parts, "snp_probe")
    if (is_non_cg_ch[i] && !is_snp[i]) parts <- c(parts, "non_cg_ch_probe")
    paste(parts, collapse = ";")
  }, character(1))

  all_rates <- data.frame(
    probe_id           = rownames(betas),
    pct_samples_failed = pct_failed,
    stringsAsFactors   = FALSE
  )

  excl <- data.frame(
    probe_id  = rownames(betas),
    call_rate = call_rate,
    is_sex    = is_sex,
    is_snp    = is_snp,
    reason    = reason,
    stringsAsFactors = FALSE
  )
  attr(excl, "all_probe_rates") <- all_rates

  # Attach probe-type counts for the QC summary
  n_total <- nrow(betas)
  n_cg    <- sum(grepl("^cg", rownames(betas)))
  n_ch    <- sum(grepl("^ch", rownames(betas)))
  n_rs    <- sum(is_snp)
  n_sex   <- sum(is_sex)
  n_other <- sum(is_non_cg_ch & !is_snp)
  n_low   <- sum(is_low)
  attr(excl, "probe_counts") <- list(
    n_total = n_total, n_cg = n_cg, n_ch = n_ch,
    n_rs = n_rs, n_sex = n_sex, n_other = n_other, n_low_call = n_low
  )

  excl <- excl[nchar(excl$reason) > 0, ]

  if (!is.null(logger)) {
    logger$log("excl_probes",
               sprintf("flagged %d (low_call=%d, sex=%d, snp=%d, other_non_cg_ch=%d)",
                       nrow(excl), n_low, n_sex, n_rs, n_other))
  }
  excl
}

#' Extract the rs (SNP) probe matrix from a beta matrix
#'
#' Returns a samples-by-rs-probes numeric matrix for identity
#' verification. Rows are samples, columns are rs probes.
#'
#' @param betas Full beta matrix (probes x samples).
#' @param logger Optional logger.
#' @return Numeric matrix (samples x rs probes), or NULL if none found.
#' @export
extract_snp_betas <- function(betas, logger = NULL) {
  rs_probes <- grep("^rs", rownames(betas), value = TRUE)
  if (length(rs_probes) == 0) {
    if (!is.null(logger)) logger$log("snp_probes", "no rs probes found")
    return(NULL)
  }
  rs_mat <- t(betas[rs_probes, , drop = FALSE])
  if (!is.null(logger)) {
    logger$log("snp_probes",
               sprintf("extracted %d rs probes x %d samples",
                       ncol(rs_mat), nrow(rs_mat)))
  }
  rs_mat
}

#' Load sex-chromosome probe IDs from the platform manifest
#' @keywords internal
#' @noRd
get_sex_probes <- function(platform) {
  tryCatch({
    mft <- sesameData::sesameDataGet(paste0(platform, ".address"))
    if (!is.null(mft$hg38)) {
      gr <- mft$hg38
      chrs <- as.character(GenomicRanges::seqnames(gr))
      return(names(gr)[chrs %in% c("chrX", "chrY")])
    }
    if (!is.null(mft$ordering) && "CpG_chrm" %in% colnames(mft$ordering)) {
      return(mft$ordering$Probe_ID[mft$ordering$CpG_chrm %in% c("chrX", "chrY")])
    }
    character(0)
  }, error = function(e) character(0))
}
