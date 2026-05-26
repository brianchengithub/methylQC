###############################################################################
# pipeline.R - Main pipeline orchestration
#
# Stage 1 (qc):
#   - discover() IDATs
#   - runsesame() to produce betas/mvals/mask/detP
#   - extract rs probe betas via snpbetas()
#   - qcmetrics() / qcstream() for per-sample QC
#   - checkmeta() for metadata cross-checks
#   - writes betas_all.rds, mvals_all.rds, mask_all.rds, detP_all.rds,
#     snp_betas.rds, sample_sheet.csv, metadata.rds
#
# Stage 2 (prep, prepcell):
#   - qcflags() for per-sample low-detection / low-intensity flags
#     (USED ONLY FOR PLOT COLOURING; nothing is excluded automatically)
#   - qcreport() writes the multi-page QC PDF + pc_scores.csv +
#     failed_probes.csv
#   - rundish() for blood tissues
#   - writes a per-cell-type sample_sheet.csv
#
# There is no exclude_samples.csv, exclude_probes.csv, or
# probe_call_rates.csv. The user applies their own exclusion criteria
# via cleanmat() and flagsamples().
###############################################################################

#' Stage 1: Preprocessing and QC metric computation
#'
#' @param indir Input directory containing IDATs (recursive).
#' @param outdir Output directory (created if needed).
#' @param platform Optional platform string; auto-detected if NULL.
#' @param ... Option overrides forwarded to \code{\link{mqcset}}.
#' @return Invisibly, a list with betas, mvals, mask, detP, sample_sheet,
#'   platform.
#' @export
qc <- function(indir, outdir, platform = NULL, ...) {
  old_opts <- mqcset(...)
  if (length(old_opts) > 0) on.exit(do.call(options, old_opts), add = TRUE)
  cfg <- mqcopts()
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  logger <- makelog(outdir)
  logger$log("start", "=== Stage 1 (qc): preprocessing ===")

  ss <- discover(indir, platform = platform, logger = logger)
  detected <- unique(stats::na.omit(ss$detected_platform))
  if (is.null(platform)) {
    if (length(detected) == 0) stop("Could not auto-detect platform.")
    if (length(detected) > 1)
      stop("Inconsistent platforms: ", paste(detected, collapse = ", "))
    platform <- detected[1]
    logger$log("start", sprintf("auto-detected platform: %s", platform))
  }

  saveRDS(list(platform        = platform,
               indir           = normalizePath(indir, mustWork = FALSE),
               n_samples       = nrow(ss),
               n_batches       = length(unique(ss$batch_folder)),
               stage1_time     = as.character(Sys.time()),
               package_version = utils::packageVersion("methylQC")),
          file.path(outdir, "metadata.rds"))

  res <- runsesame(ss$Basename, platform,
                   keepsdf = cfg$savesdf, logger = logger)

  saveRDS(res$betas, file.path(outdir, "betas_all.rds"))
  saveRDS(res$mvals, file.path(outdir, "mvals_all.rds"))
  if (!is.null(res$mask)) saveRDS(res$mask, file.path(outdir, "mask_all.rds"))
  if (!is.null(res$detP)) saveRDS(res$detP, file.path(outdir, "detP_all.rds"))
  if (cfg$savesdf && !is.null(res$sdfs))
    saveRDS(res$sdfs, file.path(outdir, "sdfs_all.rds"))

  rs_mat <- snpbetas(res$betas, logger = logger)
  if (!is.null(rs_mat)) {
    saveRDS(rs_mat, file.path(outdir, "snp_betas.rds"))
    logger$log("start",
               sprintf("wrote snp_betas.rds (%d x %d)",
                       nrow(rs_mat), ncol(rs_mat)))
  }
  rm(rs_mat); gc(verbose = FALSE)

  qcdf <- if (cfg$savesdf && !is.null(res$sdfs)) {
    qcmetrics(res$sdfs, platform = platform, logger = logger)
  } else {
    qcstream(ss$Basename, platform, logger = logger)
  }
  if (!is.null(res$sdfs)) { res$sdfs <- NULL; gc(verbose = FALSE) }

  sample_sheet <- dplyr::left_join(ss, qcdf, by = "sample_id")
  checkmeta(sample_sheet, qcdf, logger = logger)
  utils::write.csv(sample_sheet,
                   file.path(outdir, "sample_sheet.csv"),
                   row.names = FALSE)

  logger$log("done", "Stage 1 complete.")
  invisible(list(betas        = res$betas,
                 mvals        = res$mvals,
                 mask         = res$mask,
                 detP         = res$detP,
                 sample_sheet = sample_sheet,
                 platform     = platform))
}

#' Stage 2 per cell type: flagging, QC report, EpiDISH, sample sheet
#'
#' Produces QC plots, optional EpiDISH cell-type proportions, and a
#' consolidated sample sheet. Performs NO automatic sample or probe
#' exclusion: only flag columns are written. The user applies exclusion
#' criteria downstream via \code{\link{cleanmat}}.
#'
#' @param celltype Cell-type label for this slice.
#' @param betas Beta matrix (probes x samples) for this slice.
#' @param mvals M-value matrix for this slice (currently unused inside
#'   prepcell; retained because callers commonly hold both matrices).
#' @param ss Sample sheet for this slice (already merged with QC metrics).
#' @param platform Platform string.
#' @param outdir Output directory for this slice.
#' @param mask Optional logical mask matrix (probes x samples).
#' @param detP Optional detection p-value matrix.
#' @param ... Option overrides forwarded to \code{\link{mqcset}}.
#' @return Invisibly, the number of samples in the slice.
#' @export
prepcell <- function(celltype, betas, mvals, ss, platform, outdir,
                     mask = NULL, detP = NULL, dish = TRUE, ...) {
  old_opts <- mqcset(...)
  if (length(old_opts) > 0) on.exit(do.call(options, old_opts), add = TRUE)
  cfg <- mqcopts()
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  logger <- makelog(outdir)
  logger$log("start",
             sprintf("=== prep for cell type: %s (n=%d) ===",
                     celltype, ncol(betas)))
  n_total <- ncol(betas)

  # --- Flag samples (plot colouring only; no CSV, no exclusion) ---
  flags <- qcflags(ss, logger = logger)
  if (nrow(flags) / n_total > 0.25) {
    logger$log("qcflags",
               sprintf("WARNING: %.0f%% flagged (low detection/intensity)",
                       100 * nrow(flags) / n_total))
  }
  n_low_det   <- sum(grepl("low_detection", flags$reason))
  n_low_int   <- sum(grepl("low_intensity", flags$reason))
  low_int_ids <- flags$sample_id[grepl("low_intensity", flags$reason)]

  # --- Masked beta matrix for the scree/PCA panels ---
  # Apply quality mask + detection p-value mask, restrict to cg/ch probes.
  # No imputation (PCA panel handles complete.cases internally).
  betasok <- cleanmat(betas, mask = mask, detP = detP,
                       pthresh = cfg$detp,
                       probes = c("cg", "ch"),
                       platform = platform,
                       impute = FALSE)

  # --- QC diagnostic report PDF ---
  qc_flags <- qcreport(ss       = ss,
                       betas    = betas,
                       betasok  = betasok,
                       flagged  = flags,
                       mask     = mask,
                       detP     = detP,
                       platform = platform,
                       pdf      = file.path(outdir, "qc_plots.pdf"),
                       pccsv    = file.path(outdir, "pc_scores.csv"),
                       tailcsv  = file.path(outdir, "failed_probes.csv"),
                       logger   = logger)
  rm(betasok); gc(verbose = FALSE)
  n_sex_mm  <- length(qc_flags$sex_mismatch_ids)
  n_sex_unc <- length(qc_flags$sex_unclear_ids)
  n_age_out <- length(qc_flags$age_outlier_ids)

  # --- EpiDISH on raw (unmasked) betas ---
  # Skipped if dish = FALSE. There is no tissue allowlist: rundish()
  # decides on its own (it errors if the reference's CpGs don't
  # overlap the matrix). Caller controls which cell types to run.
  cell_props <- NULL
  if (!isTRUE(dish)) {
    logger$log("rundish",
               sprintf("skipping EpiDISH (dish = FALSE) for '%s'", celltype))
  } else {
    cell_props <- tryCatch(
      rundish(betas, platform = platform, logger = logger),
      error = function(e) {
        logger$log("rundish",
                   sprintf("EpiDISH skipped for '%s': %s",
                           celltype, conditionMessage(e)))
        NULL
      })
    if (!is.null(cell_props)) {
      utils::write.csv(data.frame(sample_id = colnames(betas), cell_props,
                                  check.names = FALSE),
                       file.path(outdir, "cell_proportions.csv"),
                       row.names = FALSE)
    }
  }

  # --- Consolidated sample sheet ---
  sample_sheet <- build_consolidated_sample_sheet(ss, flags, qc_flags,
                                                  cell_props, cfg)
  utils::write.csv(sample_sheet, file.path(outdir, "sample_sheet.csv"),
                   row.names = FALSE)

  # --- QC summary printed to console ---
  message(""); message("============ QC SUMMARY ============")
  message(sprintf("Total samples:                  %d", n_total))
  message("--- Sample flagging (call rate + intensity) ---")
  message(sprintf("  Low detection fraction:       %d / %d (%.1f%%)",
                  n_low_det, n_total, 100 * n_low_det / n_total))
  message(sprintf("  Low mean intensity:           %d / %d (%.1f%%)",
                  n_low_int, n_total, 100 * n_low_int / n_total))
  if (n_low_int > 0)
    message(sprintf("    IDs: %s", paste(low_int_ids, collapse = ", ")))
  message("--- Sex check ---")
  message(sprintf("  Mismatches:                   %d", n_sex_mm))
  if (n_sex_mm > 0)
    message(sprintf("    IDs: %s",
                    paste(qc_flags$sex_mismatch_ids, collapse = ", ")))
  message(sprintf("  Predicted sex unclear:        %d", n_sex_unc))
  if (n_sex_unc > 0)
    message(sprintf("    IDs: %s",
                    paste(qc_flags$sex_unclear_ids, collapse = ", ")))
  message("--- Age check (reported vs Horvath) ---")
  message(sprintf("  Age outliers (>3 SD):         %d", n_age_out))
  if (n_age_out > 0)
    message(sprintf("    IDs: %s",
                    paste(qc_flags$age_outlier_ids, collapse = ", ")))

  # --- Probe breakdown computed inline (no CSV) ---
  probe_ids  <- rownames(betas)
  is_cg      <- grepl("^cg", probe_ids)
  is_ch      <- grepl("^ch", probe_ids)
  is_rs      <- grepl("^rs", probe_ids)
  is_other   <- !is_cg & !is_ch & !is_rs
  sex_pids   <- sexprobes(platform)
  is_sex     <- probe_ids %in% sex_pids
  is_cg_auto <- is_cg & !is_sex
  is_cg_sex  <- is_cg & is_sex

  det_mask  <- if (!is.null(detP))
                 detP > cfg$detp
               else
                 matrix(FALSE, nrow = nrow(betas), ncol = ncol(betas))
  qual_mask <- if (!is.null(mask))
                 mask & !det_mask
               else
                 matrix(FALSE, nrow = nrow(betas), ncol = ncol(betas))
  combined_mask <- det_mask | qual_mask

  probe_pass_rate <- 1 - rowMeans(combined_mask)
  probe_passes    <- probe_pass_rate >= cfg$probemin

  n_total_probes <- nrow(betas)
  n_rs    <- sum(is_rs)
  n_other <- sum(is_other)

  fmt <- function(pass, total) sprintf("%s / %s (%.1f%%)",
    format(pass, big.mark = ","), format(total, big.mark = ","),
    if (total > 0) 100 * pass / total else 0)

  message("--- Probe breakdown ---")
  message(sprintf("  Total probes on array:        %s",
                  format(n_total_probes, big.mark = ",")))
  message(sprintf("    Autosomal CpG (cg):         %s",
                  fmt(sum(is_cg_auto & probe_passes), sum(is_cg_auto))))
  message(sprintf("    Non-CpG (ch):               %s",
                  fmt(sum(is_ch & probe_passes), sum(is_ch))))
  message(sprintf("    Sex chromosome (cg):        %s",
                  fmt(sum(is_cg_sex & probe_passes), sum(is_cg_sex))))
  message(sprintf("    SNP (rs):                   %d", n_rs))
  message(sprintf("    Other non-cg/ch:            %d", n_other))

  n_cells  <- length(combined_mask)
  n_det    <- sum(det_mask)
  n_qual   <- sum(qual_mask)
  pct_det  <- 100 * n_det  / n_cells
  pct_qual <- 100 * n_qual / n_cells
  message(sprintf("  Masking (detection p > %.2f): %.2f%% of probe-sample values",
                  cfg$detp, pct_det))
  message(sprintf("  Masking (quality/design):     %.2f%% of probe-sample values",
                  pct_qual))
  message("====================================")
  message("")

  summary_line <- sprintf(
    "[%s] summary :: total=%d low_det=%d low_int=%d sex_mm=%d sex_unc=%d age_out=%d",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    n_total, n_low_det, n_low_int, n_sex_mm, n_sex_unc, n_age_out)
  cat(summary_line, "\n", file = logger$path, append = TRUE)

  rm(sample_sheet, combined_mask, det_mask, qual_mask); gc(verbose = FALSE)
  done_line <- sprintf("[%s] done :: prep for %s complete.",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"), celltype)
  cat(done_line, "\n", file = logger$path, append = TRUE)
  invisible(n_total)
}

#' Build consolidated sample sheet
#' @keywords internal
#' @noRd
build_consolidated_sample_sheet <- function(ss, flags, qc_flags,
                                            cell_props, cfg) {
  out <- ss
  out$flagged <- out$sample_id %in% flags$sample_id
  m_excl <- match(out$sample_id, flags$sample_id)
  out$flag_reason <- ifelse(is.na(m_excl), "", flags$reason[m_excl])

  sex_col <- resolve_column(colnames(out), cfg$sexcol, cfg$sexaliases)
  if (!is.na(sex_col))
    out$reported_sex_normalized <- normalize_sex(out[[sex_col]],
                                                 col_name = sex_col)
  if (!is.null(qc_flags$sex_detail) && nrow(qc_flags$sex_detail) > 0) {
    m_sex <- match(out$sample_id, qc_flags$sex_detail$sample_id)
    out$inferred_sex_intensity <-
      qc_flags$sex_detail$inferred_sex_intensity[m_sex]
  } else {
    out$inferred_sex_intensity <- NA_character_
  }
  out$sex_mismatch <- out$sample_id %in% qc_flags$sex_mismatch_ids
  out$sex_unclear  <- out$sample_id %in% qc_flags$sex_unclear_ids
  out$age_outlier  <- out$sample_id %in% qc_flags$age_outlier_ids

  age_col <- resolve_column(colnames(out), cfg$agecol, cfg$agealiases)
  if (!is.na(age_col))
    out$reported_age_parsed <- parse_age_robust(out[[age_col]])

  if (!is.null(cell_props)) {
    props_df <- data.frame(sample_id = rownames(cell_props), cell_props,
                           stringsAsFactors = FALSE, check.names = FALSE)
    out <- dplyr::left_join(out, props_df, by = "sample_id")
  }
  out
}

#' Stage 2 wrapper: optionally split by cell type
#'
#' @param dir Stage 1 output directory.
#' @param platform Optional platform; defaults to value in metadata.rds.
#' @param bycell If TRUE, split processing by detected cell type.
#' @param dish If TRUE (default), run EpiDISH deconvolution within each
#'   prepcell() call. There is no built-in tissue allowlist; rundish()
#'   errors gracefully (and is caught) if the configured reference's
#'   CpGs don't sufficiently overlap the matrix. Set FALSE to skip
#'   entirely and call \code{\link{rundish}} on the output directory
#'   afterwards.
#' @param ... Option overrides forwarded to \code{\link{mqcset}}.
#' @export
prep <- function(dir, platform = NULL, bycell = TRUE, dish = TRUE, ...) {
  cfg <- mqcopts()
  metadata <- readRDS(file.path(dir, "metadata.rds"))
  if (is.null(platform)) platform <- metadata$platform

  ss <- utils::read.csv(file.path(dir, "sample_sheet.csv"),
                        stringsAsFactors = FALSE)
  mask_path <- file.path(dir, "mask_all.rds")
  detP_path <- file.path(dir, "detP_all.rds")
  mask_all <- if (file.exists(mask_path)) readRDS(mask_path) else NULL
  detP_all <- if (file.exists(detP_path)) readRDS(detP_path) else NULL

  cell_col <- resolve_column(colnames(ss), cfg$cellcol, cfg$cellaliases)

  if (!bycell || is.na(cell_col) ||
      length(unique(stats::na.omit(ss[[cell_col]]))) <= 1) {
    message("Running single-cohort prep.")
    betas <- readRDS(file.path(dir, "betas_all.rds"))
    mvals <- readRDS(file.path(dir, "mvals_all.rds"))
    ids <- intersect(ss$sample_id, colnames(betas))
    mask_s <- if (!is.null(mask_all)) mask_all[, ids, drop = FALSE] else NULL
    detP_s <- if (!is.null(detP_all)) detP_all[, ids, drop = FALSE] else NULL
    prepcell("all",
             betas[, ids, drop = FALSE], mvals[, ids, drop = FALSE],
             ss[match(ids, ss$sample_id), ],
             platform, dir, mask = mask_s, detP = detP_s, dish = dish, ...)
    rm(betas, mvals, mask_s, detP_s); gc(verbose = FALSE)
    return(invisible(NULL))
  }

  cts <- sort(unique(stats::na.omit(ss[[cell_col]])))
  message(sprintf("Splitting by '%s': %d cell types", cell_col, length(cts)))
  message("Loading matrices...")
  ba <- readRDS(file.path(dir, "betas_all.rds"))
  ma <- readRDS(file.path(dir, "mvals_all.rds"))
  gc(verbose = FALSE)
  message(sprintf("Loaded: %d probes x %d samples", nrow(ba), ncol(ba)))

  sdf <- data.frame(cell_type = character(0),
                    n_input   = integer(0),
                    n_final   = integer(0),
                    status    = character(0),
                    stringsAsFactors = FALSE)

  for (i in seq_along(cts)) {
    ct <- cts[i]
    ids <- intersect(ss$sample_id[!is.na(ss[[cell_col]]) &
                                  ss[[cell_col]] == ct],
                     colnames(ba))
    if (length(ids) < 2) {
      message(sprintf("[%d/%d] skipping %s", i, length(cts), ct))
      sdf <- rbind(sdf,
                   data.frame(cell_type = ct, n_input = length(ids),
                              n_final = length(ids), status = "skipped",
                              stringsAsFactors = FALSE))
      next
    }
    message(sprintf("[%d/%d] processing %s (n=%d)...",
                    i, length(cts), ct, length(ids)))
    mask_s <- if (!is.null(mask_all)) mask_all[, ids, drop = FALSE] else NULL
    detP_s <- if (!is.null(detP_all)) detP_all[, ids, drop = FALSE] else NULL
    nf <- tryCatch(prepcell(ct,
                            ba[, ids, drop = FALSE], ma[, ids, drop = FALSE],
                            ss[match(ids, ss$sample_id), ],
                            platform, file.path(dir, ct),
                            mask = mask_s, detP = detP_s, dish = dish, ...),
                   error = function(e) {
                     message("  ERROR: ", e$message); NA_integer_
                   })
    sdf <- rbind(sdf,
                 data.frame(cell_type = ct, n_input = length(ids),
                            n_final = if (is.null(nf) || is.na(nf))
                                       NA_integer_ else as.integer(nf),
                            status = if (is.null(nf) || is.na(nf))
                                       "failed" else "done",
                            stringsAsFactors = FALSE))
    rm(mask_s, detP_s); gc(verbose = FALSE)
  }

  rm(ba, ma, mask_all, detP_all); gc(verbose = FALSE)
  utils::write.csv(sdf, file.path(dir, "cell_type_summary.csv"),
                   row.names = FALSE)
  message(sprintf("Wrote cell_type_summary.csv (%d rows)", nrow(sdf)))
  invisible(sdf)
}

#' End-to-end pipeline
#'
#' @param indir Input directory (passed to \code{\link{qc}}).
#' @param outdir Output directory (passed to \code{\link{qc}} and used as
#'   the Stage 1 directory for \code{\link{prep}}).
#' @param platform Optional platform string.
#' @param bycell If TRUE, Stage 2 splits by cell type.
#' @param dish If TRUE (default), run EpiDISH deconvolution per cell
#'   type during Stage 2. Set to FALSE to skip; call
#'   \code{\link{rundish}} on the output directory afterwards.
#' @param ... Option overrides forwarded to \code{\link{mqcset}}.
#' @export
pipeline <- function(indir, outdir, platform = NULL, bycell = TRUE,
                     dish = TRUE, ...) {
  message("=== Running end-to-end pipeline ===")
  qc(indir = indir, outdir = outdir, platform = platform, ...)
  prep(dir = outdir, platform = platform, bycell = bycell, dish = dish, ...)
}

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
