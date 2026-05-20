###############################################################################
# exclusions.R — Sample flagging, SNP probe extraction, sex-probe lookup
#
# methylQC NEVER removes samples or probes automatically. This file holds:
#
#   qcflags()     — internal. Per-sample low-detection / low-intensity flags
#                   used ONLY to colour QC-report plots. Writes no file and
#                   feeds no exclusion list.
#   flagsamples() — exported. A standalone utility: given a detection p-value
#                   matrix and a call-rate threshold, writes a CSV listing
#                   samples below the threshold. Advisory only — the user
#                   decides whether to feed those IDs to cleanmat(dropsamples=).
#   snpbetas()    — exported. Extracts the rs (SNP) probe sub-matrix.
#   sexprobes()   — internal. Sex-chromosome probe IDs from the manifest.
#
# There is no probe-exclusion table and no exclude_probes.csv. Probe-level
# masking is carried entirely by the mask and detP matrices (see cleanmat()).
###############################################################################

#' Per-sample QC flags for plot colouring (internal)
#'
#' Produces low-detection and low-intensity flags from the merged sample
#' sheet. Used only to colour the detection-rate and intensity QC plots.
#' Sex mismatches and age outliers are deliberately NOT flagged here —
#' they are computed by, and belong to, their own QC-report pages.
#'
#' This function writes nothing to disk and produces no exclusion list.
#'
#' @param ss Sample sheet merged with QC metrics.
#' @param logger Optional logger.
#' @return A \code{data.frame} (one row per flagged sample) with columns
#'   sample_id, frac_dt, frac_usable, reason. Empty if nothing is flagged.
#' @keywords internal
#' @noRd
qcflags <- function(ss, logger = NULL) {
  cfg  <- mqcopts()
  cols <- colnames(ss)

  has_intensity <- "mean_intensity" %in% cols
  has_frac_dt   <- "frac_dt" %in% cols
  has_batch     <- "batch_folder" %in% cols

  if (!has_frac_dt) {
    stop("Sample sheet does not contain `frac_dt` column. ",
         "Did you forget to merge in qcmetrics() output?")
  }

  reasons <- vapply(seq_len(nrow(ss)), function(i) {
    row <- ss[i, , drop = FALSE]
    parts <- character(0)
    if (!is.na(row$frac_dt) && row$frac_dt < cfg$samplemin) {
      parts <- c(parts, sprintf("low_detection(%.3f)", row$frac_dt))
    }
    if (has_intensity && !is.na(row$mean_intensity) &&
        row$mean_intensity < cfg$intmin) {
      parts <- c(parts, sprintf("low_intensity(%.0f)", row$mean_intensity))
    }
    paste(parts, collapse = ";")
  }, character(1))

  base_cols <- list(
    sample_id   = ss$sample_id,
    frac_dt     = ss$frac_dt,
    frac_usable = ss$frac_usable,
    reason      = reasons
  )
  if (has_batch) {
    base_cols <- c(
      list(sample_id = ss$sample_id, batch_folder = ss$batch_folder),
      base_cols[-1]
    )
  }
  flags <- do.call(data.frame, c(base_cols, stringsAsFactors = FALSE))
  flags <- flags[nchar(flags$reason) > 0, ]

  if (!is.null(logger)) {
    logger$log("qcflags",
               sprintf("flagged %d / %d (%.1f%%) for low detection/intensity",
                       nrow(flags), nrow(ss),
                       100 * nrow(flags) / max(nrow(ss), 1)))
  }
  flags
}

#' Flag samples by detection call rate
#'
#' A standalone, opt-in utility. Computes each sample's call rate from a
#' detection p-value matrix (fraction of probes with \code{detP <= pthresh})
#' and writes a CSV of every sample whose call rate falls below
#' \code{callrate}.
#'
#' This function performs NO automatic exclusion. It only produces a list.
#' To act on it, the user passes the IDs themselves to
#' \code{\link{cleanmat}(dropsamples = ...)} — either the returned vector or a
#' CSV they assemble. \code{flagsamples()} does not accept hand-picked IDs;
#' it flags purely on the call-rate threshold.
#'
#' @param detP Numeric detection p-value matrix (probes x samples), e.g.
#'   \code{readRDS("detP_all.rds")}.
#' @param callrate Minimum acceptable call rate. Samples below this are
#'   flagged. Default: the \code{samplemin} option (0.95).
#' @param pthresh Detection p-value cutoff; a probe is counted as detected
#'   when \code{detP <= pthresh}. Default: the \code{detp} option (0.05).
#' @param csv Output CSV path. Default \code{"flagged_samples.csv"}.
#' @param logger Optional logger.
#' @return Invisibly, a data.frame of ALL samples with columns sample_id,
#'   call_rate, and flagged (logical). The CSV contains only flagged rows.
#' @examples
#' \dontrun{
#' detP <- readRDS("results/detP_all.rds")
#' flagged <- flagsamples(detP, callrate = 0.95,
#'                        csv = "results/flagged_samples.csv")
#' bad_ids <- flagged$sample_id[flagged$flagged]
#'
#' betas <- readRDS("results/betas_all.rds")
#' mask  <- readRDS("results/mask_all.rds")
#' betas_clean <- cleanmat(betas, mask = mask, detP = detP,
#'                          dropsamples = bad_ids, platform = "EPIC")
#' }
#' @export
flagsamples <- function(detP, callrate = NULL, pthresh = NULL,
                        csv = "flagged_samples.csv", logger = NULL) {
  stopifnot(is.matrix(detP) && is.numeric(detP))
  cfg      <- mqcopts()
  callrate <- if (is.null(callrate)) cfg$samplemin else callrate
  pthresh  <- if (is.null(pthresh))  cfg$detp      else pthresh

  call_rate <- colMeans(detP <= pthresh, na.rm = TRUE)
  sample_id <- colnames(detP)
  if (is.null(sample_id)) sample_id <- as.character(seq_along(call_rate))

  out <- data.frame(sample_id = sample_id,
                    call_rate = call_rate,
                    flagged   = call_rate < callrate,
                    stringsAsFactors = FALSE)
  flagged <- out[out$flagged, c("sample_id", "call_rate"), drop = FALSE]

  if (!is.null(csv)) {
    utils::write.csv(flagged, csv, row.names = FALSE)
  }
  if (!is.null(logger)) {
    logger$log("flagsamples",
               sprintf("call rate < %.3f: %d / %d sample(s); wrote %s",
                       callrate, nrow(flagged), nrow(out),
                       if (is.null(csv)) "(no file)" else csv))
  }
  message(sprintf(
    "flagsamples(): %d / %d sample(s) below call rate %.3f. No samples removed.",
    nrow(flagged), nrow(out), callrate))
  invisible(out)
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
snpbetas <- function(betas, logger = NULL) {
  rs_probes <- grep("^rs", rownames(betas), value = TRUE)
  if (length(rs_probes) == 0) {
    if (!is.null(logger)) logger$log("snpbetas", "no rs probes found")
    return(NULL)
  }
  rs_mat <- t(betas[rs_probes, , drop = FALSE])
  if (!is.null(logger)) {
    logger$log("snpbetas",
               sprintf("extracted %d rs probes x %d samples",
                       ncol(rs_mat), nrow(rs_mat)))
  }
  rs_mat
}

#' Load sex-chromosome probe IDs from the platform manifest
#' @keywords internal
#' @noRd
sexprobes <- function(platform) {
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
