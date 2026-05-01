###############################################################################
# pipeline.R — Main pipeline orchestration
###############################################################################

#' Stage 1: Preprocessing and QC metric computation
#' @export
qc <- function(in_dir, out_dir, platform = NULL, ...) {
  old_opts <- methylQC_set(...); if (length(old_opts) > 0) on.exit(do.call(options, old_opts), add = TRUE)
  cfg <- methylQC_options(); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  logger <- make_logger(out_dir); logger$log("start", "=== Stage 1 (qc): preprocessing ===")
  ss <- discover_idats(in_dir, expected_platform = platform, logger = logger)
  detected <- unique(stats::na.omit(ss$detected_platform))
  if (is.null(platform)) { if (length(detected) == 0) stop("Could not auto-detect platform.")
    if (length(detected) > 1) stop("Inconsistent platforms: ", paste(detected, collapse = ", "))
    platform <- detected[1]; logger$log("start", sprintf("auto-detected platform: %s", platform)) }
  saveRDS(list(platform = platform, in_dir = normalizePath(in_dir, mustWork = FALSE),
    n_samples = nrow(ss), n_batches = length(unique(ss$batch_folder)),
    stage1_time = as.character(Sys.time()), package_version = utils::packageVersion("methylQC")),
    file.path(out_dir, "metadata.rds"))
  res <- run_opensesame(ss$Basename, platform, keep_sdfs = cfg$save_sdfs, logger = logger)
  saveRDS(res$betas, file.path(out_dir, "betas_all.rds"))
  saveRDS(res$mvals, file.path(out_dir, "mvals_all.rds"))
  if (!is.null(res$mask)) saveRDS(res$mask, file.path(out_dir, "mask_all.rds"))
  if (!is.null(res$detP)) saveRDS(res$detP, file.path(out_dir, "detP_all.rds"))
  if (cfg$save_sdfs && !is.null(res$sdfs)) saveRDS(res$sdfs, file.path(out_dir, "sdfs_all.rds"))
  rs_mat <- extract_snp_betas(res$betas, logger = logger)
  if (!is.null(rs_mat)) { saveRDS(rs_mat, file.path(out_dir, "snp_betas.rds"))
    logger$log("start", sprintf("wrote snp_betas.rds (%d x %d)", nrow(rs_mat), ncol(rs_mat))) }
  rm(rs_mat); gc(verbose = FALSE)
  qc_df <- if (cfg$save_sdfs && !is.null(res$sdfs)) compute_qc_metrics(res$sdfs, platform = platform, logger = logger)
    else compute_qc_metrics_streaming(ss$Basename, platform, logger = logger)
  if (!is.null(res$sdfs)) { res$sdfs <- NULL; gc(verbose = FALSE) }
  sample_sheet <- dplyr::left_join(ss, qc_df, by = "sample_id")
  check_sample_metadata(sample_sheet, qc_df, logger = logger)
  utils::write.csv(sample_sheet, file.path(out_dir, "sample_sheet.csv"), row.names = FALSE)
  logger$log("done", "Stage 1 complete.")
  invisible(list(betas = res$betas, mvals = res$mvals, mask = res$mask, detP = res$detP,
                 sample_sheet = sample_sheet, platform = platform))
}

#' Stage 2: Flagging, QC report, EpiDISH, consolidated sample sheet
#' @param mask_sub Logical matrix (probes x samples). NULL if not available.
#' @param detP_sub Numeric matrix of detection p-values. NULL if not available.
#' @export
prep_single <- function(cell_type, betas_sub, mvals_sub, ss_sub, platform, out_subdir,
                        mask_sub = NULL, detP_sub = NULL, ...) {
  old_opts <- methylQC_set(...); if (length(old_opts) > 0) on.exit(do.call(options, old_opts), add = TRUE)
  cfg <- methylQC_options(); dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  logger <- make_logger(out_subdir)
  logger$log("start", sprintf("=== prep for cell type: %s (n=%d) ===", cell_type, ncol(betas_sub)))
  n_total <- ncol(betas_sub)

  # --- Flag samples and probes ---
  excl_s <- build_sample_exclusions(ss_sub, logger = logger)
  excl_p <- build_probe_exclusions(betas_sub, platform, logger = logger)
  if (nrow(excl_s) / n_total > 0.25)
    logger$log("excl_samples", sprintf("WARNING: %.0f%% flagged", 100 * nrow(excl_s) / n_total))
  apr <- attr(excl_p, "all_probe_rates")
  if (!is.null(apr)) utils::write.csv(apr, file.path(out_subdir, "probe_call_rates.csv"), row.names = FALSE)
  utils::write.csv(excl_s, file.path(out_subdir, "exclude_samples.csv"), row.names = FALSE)
  utils::write.csv(excl_p, file.path(out_subdir, "exclude_probes.csv"), row.names = FALSE)
  n_low_det <- sum(grepl("low_detection", excl_s$reason))
  n_low_int <- sum(grepl("low_intensity", excl_s$reason))
  low_int_ids <- excl_s$sample_id[grepl("low_intensity", excl_s$reason)]

  # --- Temporary filtered matrix for QC plots (not saved to disk) ---
  keep_s <- setdiff(colnames(betas_sub), excl_s$sample_id)
  keep_p <- setdiff(rownames(betas_sub), excl_p$probe_id)
  betas_qc <- betas_sub[keep_p, keep_s, drop = FALSE]
  n_final <- ncol(betas_qc)
  logger$log("apply_excl", sprintf("retained %d probes x %d samples (for QC plots)",
             nrow(betas_qc), n_final))

  # --- QC diagnostic report PDF ---
  qc_flags <- write_qc_report(ss_all = ss_sub, betas_all_sub = betas_sub, betas_qc = betas_qc,
    excluded_ids = excl_s$sample_id, platform = platform,
    pdf_path = file.path(out_subdir, "qc_plots.pdf"),
    pc_csv_path = file.path(out_subdir, "pc_scores.csv"),
    probe_csv_path = file.path(out_subdir, "probe_call_rates.csv"), logger = logger)
  n_sex_mm <- length(qc_flags$sex_mismatch_ids); n_sex_unc <- length(qc_flags$sex_unclear_ids)
  n_age_out <- length(qc_flags$age_outlier_ids)

  rm(betas_qc); gc(verbose = FALSE)

  # --- EpiDISH on raw (unmasked) betas ---
  resolved_ct <- infer_cell_type_label(cell_type, ss_sub, cfg, logger)
  cell_props <- NULL
  if (is_blood_tissue(resolved_ct, cfg$epidish_cell_types)) {
    cell_props <- run_epidish(betas_sub, logger = logger)
    utils::write.csv(data.frame(sample_id = colnames(betas_sub), cell_props,
                                 check.names = FALSE),
      file.path(out_subdir, "cell_proportions.csv"), row.names = FALSE)
  } else {
    logger$log("epidish", sprintf("skipping EpiDISH for '%s'", resolved_ct))
  }

  # --- Consolidated sample sheet ---
  sample_sheet <- build_consolidated_sample_sheet(ss_sub, excl_s, qc_flags, cell_props, cfg)
  utils::write.csv(sample_sheet, file.path(out_subdir, "sample_sheet.csv"), row.names = FALSE)

  # --- QC Summary ---
  message(""); message("============ QC SUMMARY ============")
  message(sprintf("Total samples:                  %d", n_total))
  message("--- Sample flagging (call rate + intensity) ---")
  message(sprintf("  Low detection fraction:       %d / %d (%.1f%%)", n_low_det, n_total, 100 * n_low_det / n_total))
  message(sprintf("  Low mean intensity:           %d / %d (%.1f%%)", n_low_int, n_total, 100 * n_low_int / n_total))
  if (n_low_int > 0) message(sprintf("    IDs: %s", paste(low_int_ids, collapse = ", ")))
  message("--- Sex check ---")
  message(sprintf("  Mismatches:                   %d", n_sex_mm))
  if (n_sex_mm > 0) message(sprintf("    IDs: %s", paste(qc_flags$sex_mismatch_ids, collapse = ", ")))
  message(sprintf("  Predicted sex unclear:        %d", n_sex_unc))
  if (n_sex_unc > 0) message(sprintf("    IDs: %s", paste(qc_flags$sex_unclear_ids, collapse = ", ")))
  message("--- Age check (reported vs Horvath) ---")
  message(sprintf("  Age outliers (>3 SD):         %d", n_age_out))
  if (n_age_out > 0) message(sprintf("    IDs: %s", paste(qc_flags$age_outlier_ids, collapse = ", ")))

  # --- Probe breakdown with masking stats ---
  pc <- attr(excl_p, "probe_counts")
  if (!is.null(pc)) {
    message("--- Probe breakdown ---")
    probe_ids <- rownames(betas_sub)
    is_cg    <- grepl("^cg", probe_ids)
    is_ch    <- grepl("^ch", probe_ids)
    is_rs    <- grepl("^rs", probe_ids)
    is_other <- !is_cg & !is_ch & !is_rs
    sex_probes <- get_sex_probes(platform)
    is_sex   <- probe_ids %in% sex_probes
    # Autosomal cg = cg probes NOT on sex chromosomes
    is_cg_auto <- is_cg & !is_sex
    # Sex cg = cg probes ON sex chromosomes
    is_cg_sex  <- is_cg & is_sex

    # Build detection mask and quality mask separately
    det_mask <- matrix(FALSE, nrow = nrow(betas_sub), ncol = ncol(betas_sub))
    qual_mask <- matrix(FALSE, nrow = nrow(betas_sub), ncol = ncol(betas_sub))
    if (!is.null(detP_sub)) det_mask <- detP_sub > cfg$detP_thresh
    if (!is.null(mask_sub)) {
      # The SeSAMe mask combines quality + detection; quality-only is
      # the portion of the mask NOT explained by detection failure
      qual_mask <- mask_sub & !det_mask
    }
    combined_mask <- det_mask | qual_mask

    # Per-probe pass rate (fraction of samples where the probe passed)
    probe_pass_rate <- 1 - rowMeans(combined_mask)
    probe_passes <- probe_pass_rate >= cfg$probe_call_min

    fmt <- function(pass, total) sprintf("%s / %s (%.1f%%)",
      format(pass, big.mark = ","), format(total, big.mark = ","),
      if (total > 0) 100 * pass / total else 0)

    message(sprintf("  Total probes on array:        %s", format(pc$n_total, big.mark = ",")))
    message(sprintf("    Autosomal CpG (cg):         %s", fmt(sum(is_cg_auto & probe_passes), sum(is_cg_auto))))
    message(sprintf("    Non-CpG (ch):               %s", fmt(sum(is_ch & probe_passes), sum(is_ch))))
    message(sprintf("    Sex chromosome (cg):        %s", fmt(sum(is_cg_sex & probe_passes), sum(is_cg_sex))))
    message(sprintf("    SNP (rs):                   %d", pc$n_rs))
    message(sprintf("    Other non-cg/ch:            %d", pc$n_other))

    # Masking breakdown: detection vs quality
    n_cells <- length(combined_mask)
    n_det   <- sum(det_mask)
    n_qual  <- sum(qual_mask)
    pct_det  <- 100 * n_det / n_cells
    pct_qual <- 100 * n_qual / n_cells
    message(sprintf("  Masking (detection p > %.2f): %.2f%% of probe-sample values", cfg$detP_thresh, pct_det))
    message(sprintf("  Masking (quality/design):     %.2f%% of probe-sample values", pct_qual))
  }
  message("===================================="); message("")

  # Write summary to log file only
  summary_line <- sprintf("[%s] summary :: total=%d low_det=%d low_int=%d sex_mm=%d sex_unc=%d age_out=%d",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n_total, n_low_det, n_low_int, n_sex_mm, n_sex_unc, n_age_out)
  cat(summary_line, "\n", file = logger$path, append = TRUE)
  rm(sample_sheet); gc(verbose = FALSE)
  done_line <- sprintf("[%s] done :: prep for %s complete.",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), cell_type)
  cat(done_line, "\n", file = logger$path, append = TRUE)
  invisible(n_final)
}

#' Build consolidated sample sheet
#' @keywords internal
#' @noRd
build_consolidated_sample_sheet <- function(ss_sub, excl_s, qc_flags, cell_props, cfg) {
  ss <- ss_sub
  ss$flagged <- ss$sample_id %in% excl_s$sample_id
  m_excl <- match(ss$sample_id, excl_s$sample_id)
  ss$flag_reason <- ifelse(is.na(m_excl), "", excl_s$reason[m_excl])
  sex_col <- resolve_column(colnames(ss), cfg$reported_sex_col, cfg$reported_sex_aliases)
  if (!is.na(sex_col)) ss$reported_sex_normalized <- normalize_sex(ss[[sex_col]], col_name = sex_col)
  if (!is.null(qc_flags$sex_detail) && nrow(qc_flags$sex_detail) > 0) {
    m_sex <- match(ss$sample_id, qc_flags$sex_detail$sample_id)
    ss$inferred_sex_intensity <- qc_flags$sex_detail$inferred_sex_intensity[m_sex]
  } else ss$inferred_sex_intensity <- NA_character_
  ss$sex_mismatch <- ss$sample_id %in% qc_flags$sex_mismatch_ids
  ss$sex_unclear  <- ss$sample_id %in% qc_flags$sex_unclear_ids
  ss$age_outlier  <- ss$sample_id %in% qc_flags$age_outlier_ids
  age_col <- resolve_column(colnames(ss), cfg$reported_age_col, cfg$reported_age_aliases)
  if (!is.na(age_col)) ss$reported_age_parsed <- parse_age_robust(ss[[age_col]])
  if (!is.null(cell_props)) {
    props_df <- data.frame(sample_id = rownames(cell_props), cell_props,
                            stringsAsFactors = FALSE, check.names = FALSE)
    ss <- dplyr::left_join(ss, props_df, by = "sample_id")
  }
  ss
}

#' Infer cell-type label from sample sheet
#' @keywords internal
#' @noRd
infer_cell_type_label <- function(cell_type, ss, cfg, logger = NULL) {
  if (!is.null(cell_type) && cell_type != "all") return(cell_type)
  cell_col <- resolve_column(colnames(ss), cfg$cell_col, cfg$cell_aliases)
  if (!is.na(cell_col)) {
    vals <- stats::na.omit(ss[[cell_col]])
    if (length(vals) > 0) {
      dom <- names(sort(table(vals), decreasing = TRUE))[1]
      if (!is.null(logger)) logger$log("epidish", sprintf("inferred '%s' from '%s'", dom, cell_col))
      return(dom)
    }
  }
  if (!is.null(logger)) logger$log("epidish", "could not infer cell type"); cell_type
}

#' Stage 2 wrapper: optionally split by cell type
#' @export
prep <- function(stage1_dir, platform = NULL, split_by_cell = TRUE, ...) {
  cfg <- methylQC_options(); metadata <- readRDS(file.path(stage1_dir, "metadata.rds"))
  if (is.null(platform)) platform <- metadata$platform
  ss <- utils::read.csv(file.path(stage1_dir, "sample_sheet.csv"), stringsAsFactors = FALSE)
  mask_path <- file.path(stage1_dir, "mask_all.rds")
  detP_path <- file.path(stage1_dir, "detP_all.rds")
  mask_all <- if (file.exists(mask_path)) readRDS(mask_path) else NULL
  detP_all <- if (file.exists(detP_path)) readRDS(detP_path) else NULL
  cell_col <- resolve_column(colnames(ss), cfg$cell_col, cfg$cell_aliases)
  if (!split_by_cell || is.na(cell_col) || length(unique(stats::na.omit(ss[[cell_col]]))) <= 1) {
    message("Running single-cohort prep."); betas <- readRDS(file.path(stage1_dir, "betas_all.rds"))
    mvals <- readRDS(file.path(stage1_dir, "mvals_all.rds")); ids <- intersect(ss$sample_id, colnames(betas))
    mask_s <- if (!is.null(mask_all)) mask_all[, ids, drop = FALSE] else NULL
    detP_s <- if (!is.null(detP_all)) detP_all[, ids, drop = FALSE] else NULL
    prep_single("all", betas[,ids,drop=FALSE], mvals[,ids,drop=FALSE], ss[match(ids, ss$sample_id),],
                platform, stage1_dir, mask_sub = mask_s, detP_sub = detP_s, ...)
    rm(betas, mvals, mask_s, detP_s); gc(verbose = FALSE); return(invisible(NULL)) }
  cts <- sort(unique(stats::na.omit(ss[[cell_col]]))); message(sprintf("Splitting by '%s': %d cell types", cell_col, length(cts)))
  message("Loading matrices..."); ba <- readRDS(file.path(stage1_dir, "betas_all.rds"))
  ma <- readRDS(file.path(stage1_dir, "mvals_all.rds")); gc(verbose = FALSE)
  message(sprintf("Loaded: %d probes x %d samples", nrow(ba), ncol(ba)))
  sdf <- data.frame(cell_type=character(0), n_input=integer(0), n_final=integer(0), status=character(0), stringsAsFactors=FALSE)
  for (i in seq_along(cts)) { ct <- cts[i]
    ids <- intersect(ss$sample_id[!is.na(ss[[cell_col]]) & ss[[cell_col]]==ct], colnames(ba))
    if (length(ids) < 2) { message(sprintf("[%d/%d] skipping %s", i, length(cts), ct))
      sdf <- rbind(sdf, data.frame(cell_type=ct, n_input=length(ids), n_final=length(ids), status="skipped", stringsAsFactors=FALSE)); next }
    message(sprintf("[%d/%d] processing %s (n=%d)...", i, length(cts), ct, length(ids)))
    mask_s <- if (!is.null(mask_all)) mask_all[, ids, drop = FALSE] else NULL
    detP_s <- if (!is.null(detP_all)) detP_all[, ids, drop = FALSE] else NULL
    nf <- tryCatch(prep_single(ct, ba[,ids,drop=FALSE], ma[,ids,drop=FALSE], ss[match(ids,ss$sample_id),],
      platform, file.path(stage1_dir, ct), mask_sub = mask_s, detP_sub = detP_s, ...),
      error = function(e) { message("  ERROR: ", e$message); NA_integer_ })
    sdf <- rbind(sdf, data.frame(cell_type=ct, n_input=length(ids),
      n_final=if(is.null(nf)||is.na(nf)) NA_integer_ else as.integer(nf),
      status=if(is.null(nf)||is.na(nf)) "failed" else "done", stringsAsFactors=FALSE))
    rm(mask_s, detP_s); gc(verbose = FALSE) }
  rm(ba, ma, mask_all, detP_all); gc(verbose = FALSE)
  utils::write.csv(sdf, file.path(stage1_dir, "cell_type_summary.csv"), row.names = FALSE)
  message(sprintf("Wrote cell_type_summary.csv (%d rows)", nrow(sdf))); invisible(sdf)
}

#' End-to-end pipeline
#' @export
pipeline <- function(in_dir, out_dir, platform = NULL, split_by_cell = TRUE, ...) {
  message("=== Running end-to-end pipeline ===")
  qc(in_dir = in_dir, out_dir = out_dir, platform = platform, ...)
  prep(stage1_dir = out_dir, platform = platform, split_by_cell = split_by_cell, ...)
}

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Edit exclusion lists post-hoc
#' @export
edit_exclusions <- function(stage1_dir, add_samples = character(0), remove_samples = character(0),
  add_probes = character(0), remove_probes = character(0), add_reason = "manual") {
  sp <- file.path(stage1_dir, "exclude_samples.csv"); pp <- file.path(stage1_dir, "exclude_probes.csv")
  if (!file.exists(sp) || !file.exists(pp)) stop("Exclusion CSVs not found")
  es <- utils::read.csv(sp, stringsAsFactors = FALSE); ep <- utils::read.csv(pp, stringsAsFactors = FALSE)
  logger <- make_logger(stage1_dir)
  es <- es[!es$sample_id %in% remove_samples, , drop = FALSE]; ep <- ep[!ep$probe_id %in% remove_probes, , drop = FALSE]
  as2 <- setdiff(add_samples, es$sample_id)
  if (length(as2)) { nr <- data.frame(sample_id = as2, stringsAsFactors = FALSE)
    for (col in setdiff(colnames(es), "sample_id")) nr[[col]] <- if (col == "reason") add_reason else NA
    es <- rbind(es, nr[, colnames(es), drop = FALSE]) }
  ap <- setdiff(add_probes, ep$probe_id)
  if (length(ap)) { nr <- data.frame(probe_id = ap, stringsAsFactors = FALSE)
    for (col in setdiff(colnames(ep), "probe_id")) nr[[col]] <- if (col == "reason") add_reason else NA
    ep <- rbind(ep, nr[, colnames(ep), drop = FALSE]) }
  utils::write.csv(es, sp, row.names = FALSE); utils::write.csv(ep, pp, row.names = FALSE)
  logger$log("edit_excl", sprintf("updated: %d samples, %d probes", nrow(es), nrow(ep)))
  invisible(list(exclude_samples = es, exclude_probes = ep))
}
