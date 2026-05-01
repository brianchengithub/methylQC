#' Write all QC diagnostic plots to a single PDF
#' @export
write_qc_report <- function(ss_all, betas_all_sub, betas_qc, excluded_ids, platform,
                             pdf_path, pc_csv_path, probe_csv_path,
                             n_top = 10000, logger = NULL) {
  cfg <- methylQC_options()
  theme_qc <- ggplot2::theme_minimal() + ggplot2::theme(
    axis.title = ggplot2::element_text(size = 14), axis.text = ggplot2::element_text(size = 12),
    plot.title = ggplot2::element_text(size = 15, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 13))
  grDevices::pdf(pdf_path, width = 10, height = 7); on.exit(grDevices::dev.off())
  ss_all$excluded <- ss_all$sample_id %in% excluded_ids
  ss_all$status <- ifelse(ss_all$excluded, "Flagged", "OK")
  sex_probes <- get_sex_probes(platform)
  donor_col <- resolve_column(colnames(ss_all), cfg$donor_col, cfg$donor_aliases)
  # Page 1: Detection rate
  theme_det <- theme_qc + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
  print(ggplot2::ggplot(ss_all, ggplot2::aes(x = stats::reorder(sample_id, frac_dt), y = frac_dt, fill = status)) +
    ggplot2::geom_col() + ggplot2::scale_fill_manual(values = c("OK" = "steelblue", "Flagged" = "red3"), name = NULL) +
    ggplot2::geom_hline(yintercept = cfg$sample_call_min, linetype = "dashed", color = "black", linewidth = 0.7) +
    ggplot2::coord_flip() + theme_det +
    ggplot2::labs(title = "Detection rate per sample (frac_dt)", subtitle = sprintf("All %d samples; red = flagged; dashed = %.2f", nrow(ss_all), cfg$sample_call_min), x = NULL, y = "Fraction detected"))
  # Page 2: MDS
  mds_outlier_ids <- plot_mds_all(betas_all_sub, ss_all, donor_col, sex_probes, n_top, theme_qc, logger)
  # Page 3: Intensity
  if ("mean_intensity" %in% colnames(ss_all)) {
    ss_all$mds_status <- ifelse(ss_all$sample_id %in% mds_outlier_ids, "MDS outlier", "OK")
    print(ggplot2::ggplot(ss_all, ggplot2::aes(mean_intensity, fill = mds_status)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.8, position = "identity") +
      ggplot2::scale_fill_manual(values = c("OK" = "steelblue", "MDS outlier" = "red3"), name = NULL) +
      ggplot2::geom_vline(xintercept = cfg$intensity_min, linetype = "dashed", color = "black", linewidth = 0.7) +
      theme_qc + ggplot2::labs(title = "Mean intensity per sample", subtitle = sprintf("All %d samples; colored by MDS outlier status", nrow(ss_all)), x = "Mean intensity", y = "Count"))
  }
  # Page 4: Probe failure
  pct_failed <- rowMeans(is.na(betas_all_sub)) * 100
  n_high_fail <- sum(pct_failed > (1 - cfg$probe_call_min) * 100)
  print(ggplot2::ggplot(data.frame(pf = pct_failed), ggplot2::aes(x = pf)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    ggplot2::geom_vline(xintercept = (1 - cfg$probe_call_min) * 100, linetype = "dashed", color = "red", linewidth = 0.7) +
    theme_qc + ggplot2::labs(title = "Per-probe sample failure rate (all probes, all samples)",
      subtitle = sprintf("%d probes total; %d (%.1f%%) with >%.0f%% NA", nrow(betas_all_sub), n_high_fail, 100 * n_high_fail / nrow(betas_all_sub), (1 - cfg$probe_call_min) * 100),
      x = "% samples with NA per probe", y = "Number of probes"))
  utils::write.csv(data.frame(probe_id = rownames(betas_qc), pct_samples_failed = rowMeans(is.na(betas_qc)) * 100, stringsAsFactors = FALSE), probe_csv_path, row.names = FALSE)
  rm(pct_failed); gc(verbose = FALSE)
  # Page 5: Sex check
  sex_result <- plot_sex_check_optimal(ss_all, theme_qc, logger)
  # Page 6: Age check
  age_result <- plot_age_check(ss_all, sex_mismatch_ids = sex_result$mismatch_ids, theme_qc, logger)
  # Pages 7+: PCA
  plot_pca_retained(betas_qc, ss_all, excluded_ids, sex_probes, n_top, pc_csv_path, theme_qc, logger)
  rm(ss_all); gc(verbose = FALSE)
  invisible(list(sex_mismatch_ids = sex_result$mismatch_ids, sex_unclear_ids = sex_result$unclear_ids, age_outlier_ids = age_result$outlier_ids, sex_detail = sex_result$detail, age_detail = age_result$detail))
}

#' Sex check: optimal threshold + colored regression bands
#' @keywords internal
#' @noRd
plot_sex_check_optimal <- function(ss_all, theme_qc, logger = NULL) {
  cfg <- methylQC_options()
  empty <- list(mismatch_ids = character(0), unclear_ids = character(0),
                detail = data.frame(sample_id = character(0), reported_sex = character(0),
                                     inferred_sex_intensity = character(0), stringsAsFactors = FALSE))
  has_x <- "sex_chrX_intensity" %in% colnames(ss_all)
  has_y <- "sex_chrY_intensity" %in% colnames(ss_all)
  if (!has_x || !has_y) return(empty)
  sex_df <- data.frame(sample_id = ss_all$sample_id, chrX = ss_all$sex_chrX_intensity,
                        chrY = ss_all$sex_chrY_intensity, stringsAsFactors = FALSE)
  sex_df <- sex_df[!is.na(sex_df$chrX) & !is.na(sex_df$chrY), ]
  if (nrow(sex_df) < 6) return(empty)
  m <- match(sex_df$sample_id, ss_all$sample_id)
  sex_col <- resolve_column(colnames(ss_all), cfg$reported_sex_col, cfg$reported_sex_aliases)
  sex_df$reported_sex <- if (!is.na(sex_col))
    normalize_sex(ss_all[[sex_col]][m], col_name = sex_col) else NA_character_
  n_missing_sex <- sum(is.na(sex_df$reported_sex))
  # Optimize chrY threshold
  candidates <- stats::quantile(sex_df$chrY, probs = seq(0.15, 0.85, by = 0.01), na.rm = TRUE)
  best_thresh <- stats::median(sex_df$chrY); best_cost <- Inf
  for (thr in candidates) {
    fi <- which(sex_df$chrY <= thr); mi <- which(sex_df$chrY > thr)
    if (length(fi) < 3 || length(mi) < 3) next
    ff <- tryCatch(stats::lm(chrY ~ chrX, data = sex_df[fi, ]), error = function(e) NULL)
    fm <- tryCatch(stats::lm(chrY ~ chrX, data = sex_df[mi, ]), error = function(e) NULL)
    if (is.null(ff) || is.null(fm)) next
    cost <- sum(abs(stats::residuals(ff))) + sum(abs(stats::residuals(fm)))
    if (cost < best_cost) { best_cost <- cost; best_thresh <- thr }
  }
  f_idx <- which(sex_df$chrY <= best_thresh); m_idx <- which(sex_df$chrY > best_thresh)
  fit_f <- if (length(f_idx) >= 3) tryCatch(stats::lm(chrY ~ chrX, data = sex_df[f_idx, ]), error = function(e) NULL) else NULL
  fit_m <- if (length(m_idx) >= 3) tryCatch(stats::lm(chrY ~ chrX, data = sex_df[m_idx, ]), error = function(e) NULL) else NULL
  # Orthogonal distance
  ortho_dist <- function(x, y, fit) {
    if (is.null(fit)) return(rep(Inf, length(x)))
    a <- stats::coef(fit)[2]; c <- stats::coef(fit)[1]
    abs(a * x - y + c) / sqrt(a^2 + 1)
  }
  d_f <- ortho_dist(sex_df$chrX, sex_df$chrY, fit_f)
  d_m <- ortho_dist(sex_df$chrX, sex_df$chrY, fit_m)
  sd_f <- if (length(f_idx) >= 3) stats::sd(d_f[f_idx]) else Inf
  sd_m <- if (length(m_idx) >= 3) stats::sd(d_m[m_idx]) else Inf
  band_sd <- 5.0  # wider bands
  thresh_f <- band_sd * sd_f; thresh_m <- band_sd * sd_m
  # Classification: purely DNA methylation based
  sex_df$inferred_sex_intensity <- ifelse(
    d_f <= thresh_f & d_m > thresh_m, "F",
    ifelse(d_m <= thresh_m & d_f > thresh_f, "M",
    ifelse(d_f <= thresh_f & d_m <= thresh_m,
           ifelse(sex_df$chrY <= best_thresh, "F", "M"), "Unclear")))
  # Mismatch: has reported sex, confident call, disagrees. Never for missing sex.
  sex_df$mismatch <- !is.na(sex_df$reported_sex) &
    sex_df$inferred_sex_intensity != "Unclear" &
    sex_df$reported_sex != sex_df$inferred_sex_intensity
  sex_df$unclear <- sex_df$inferred_sex_intensity == "Unclear"
  mismatch_ids <- sex_df$sample_id[sex_df$mismatch]
  unclear_ids  <- sex_df$sample_id[sex_df$unclear]
  if (!is.null(logger)) {
    logger$log("qc_plots", sprintf("sex clustering: threshold=%.0f, F=%d, M=%d, mismatch=%d, unclear=%d, missing_sex=%d",
      best_thresh, length(f_idx), length(m_idx), length(mismatch_ids), length(unclear_ids), n_missing_sex))
    if (length(mismatch_ids) > 0) logger$log("qc_plots", sprintf("sex mismatch IDs: %s", paste(mismatch_ids, collapse = ", ")))
    if (length(unclear_ids) > 0) logger$log("qc_plots", sprintf("sex unclear IDs: %s", paste(unclear_ids, collapse = ", ")))
  }
  detail <- sex_df[, c("sample_id", "reported_sex", "inferred_sex_intensity"), drop = FALSE]
  # Plot: F, M, n/a in legend
  sex_df$color_var <- ifelse(is.na(sex_df$reported_sex), "n/a", sex_df$reported_sex)
  sex_df$color_var <- factor(sex_df$color_var, levels = c("F", "M", "n/a"))
  p <- ggplot2::ggplot(sex_df, ggplot2::aes(chrX, chrY, color = color_var, shape = color_var)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::scale_color_manual(values = c("F" = "#E41A1C", "M" = "#377EB8", "n/a" = "grey50"), name = "Reported sex") +
    ggplot2::scale_shape_manual(values = c("F" = 16, "M" = 17, "n/a" = 1), name = "Reported sex")
  # Colored confidence bands
  x_range <- range(sex_df$chrX)
  x_seq <- seq(x_range[1] - diff(x_range) * 0.05, x_range[2] + diff(x_range) * 0.05, length.out = 200)
  for (cl in list(list(fit = fit_f, thresh = thresh_f, col = "#E41A1C"),
                  list(fit = fit_m, thresh = thresh_m, col = "#377EB8"))) {
    if (is.null(cl$fit) || !is.finite(cl$thresh)) next
    a <- stats::coef(cl$fit)[2]; b_int <- stats::coef(cl$fit)[1]
    y_pred <- b_int + a * x_seq; y_off <- cl$thresh / cos(atan(a))
    band_df <- data.frame(x = c(x_seq, rev(x_seq)), y = c(y_pred + y_off, rev(y_pred - y_off)))
    p <- p + ggplot2::geom_polygon(data = band_df, ggplot2::aes(x, y),
      inherit.aes = FALSE, fill = cl$col, alpha = 0.20)
  }
  subtitle_txt <- sprintf("All %d samples; %d mismatch(es), %d predicted sex unclear",
    nrow(sex_df), length(mismatch_ids), length(unclear_ids))
  if (n_missing_sex > 0) subtitle_txt <- paste0(subtitle_txt, sprintf("; %d missing reported sex", n_missing_sex))
  p <- p + theme_qc + ggplot2::labs(title = "Sex check: chrX vs chrY median intensity (curated probes)",
    subtitle = subtitle_txt, x = "Median total intensity (chrX.xlinked)", y = "Median total intensity (chrY.clean)")
  print(p)
  list(mismatch_ids = mismatch_ids, unclear_ids = unclear_ids, detail = detail)
}

#' Age check: regression band, colored by sex mismatch only
#' @keywords internal
#' @noRd
plot_age_check <- function(ss_all, sex_mismatch_ids, theme_qc, logger = NULL) {
  cfg <- methylQC_options()
  empty <- list(outlier_ids = character(0), detail = data.frame(sample_id = character(0),
    reported_age = numeric(0), horvath_age = numeric(0), ortho_dist = numeric(0), stringsAsFactors = FALSE))
  age_col <- resolve_column(colnames(ss_all), cfg$reported_age_col, cfg$reported_age_aliases)
  if (is.na(age_col)) { if (!is.null(logger)) logger$log("qc_plots", "age check skipped: no age column"); return(empty) }
  if (!"horvath_age" %in% colnames(ss_all)) { if (!is.null(logger)) logger$log("qc_plots", "age check skipped: no horvath_age"); return(empty) }
  age_df <- data.frame(sample_id = ss_all$sample_id, reported = parse_age_robust(ss_all[[age_col]]),
                        horvath = ss_all$horvath_age, stringsAsFactors = FALSE)
  if (!is.null(logger)) logger$log("qc_plots", sprintf("age check: %d reported, %d horvath",
    sum(!is.na(age_df$reported)), sum(!is.na(age_df$horvath))))
  age_df <- age_df[!is.na(age_df$reported) & !is.na(age_df$horvath), ]
  if (nrow(age_df) < 3) { if (!is.null(logger)) logger$log("qc_plots", sprintf("age check skipped: %d valid pairs", nrow(age_df))); return(empty) }
  fit <- stats::lm(horvath ~ reported, data = age_df)
  a <- stats::coef(fit)[2]; b <- stats::coef(fit)[1]
  age_df$ortho_dist <- abs(a * age_df$reported - age_df$horvath + b) / sqrt(a^2 + 1)
  dist_sd <- stats::sd(age_df$ortho_dist)
  age_df$age_outlier <- age_df$ortho_dist > 3 * dist_sd
  outlier_ids <- age_df$sample_id[age_df$age_outlier]
  detail <- age_df[age_df$age_outlier, c("sample_id", "reported", "horvath", "ortho_dist"), drop = FALSE]
  colnames(detail) <- c("sample_id", "reported_age", "horvath_age", "ortho_dist")
  r <- suppressWarnings(stats::cor(age_df$reported, age_df$horvath))
  mae <- mean(abs(age_df$reported - age_df$horvath))
  # Color ONLY by sex mismatch (no unclear)
  age_df$color_var <- ifelse(age_df$sample_id %in% sex_mismatch_ids, "Sex mismatch", "OK")
  age_df$color_var <- factor(age_df$color_var, levels = c("OK", "Sex mismatch"))
  n_mm_in_age <- sum(age_df$sample_id %in% sex_mismatch_ids)
  if (!is.null(logger)) logger$log("qc_plots", sprintf("age plot: %d/%d sex mismatches have valid ages",
    n_mm_in_age, length(sex_mismatch_ids)))
  # Identical axes
  all_vals <- c(age_df$reported, age_df$horvath)
  ax_lim <- range(all_vals) + c(-1, 1) * diff(range(all_vals)) * 0.05
  # Blue confidence band
  x_seq <- seq(ax_lim[1], ax_lim[2], length.out = 200)
  y_pred <- b + a * x_seq; y_off <- 3 * dist_sd / cos(atan(a))
  band_df <- data.frame(x = c(x_seq, rev(x_seq)), y = c(y_pred + y_off, rev(y_pred - y_off)))
  print(ggplot2::ggplot(age_df, ggplot2::aes(reported, horvath, color = color_var)) +
    ggplot2::geom_polygon(data = band_df, ggplot2::aes(x, y), inherit.aes = FALSE, fill = "steelblue", alpha = 0.15) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +
    ggplot2::geom_point(size = 2.5, alpha = 0.8) +
    ggplot2::scale_color_manual(values = c("OK" = "steelblue", "Sex mismatch" = "red3"), name = NULL) +
    ggplot2::coord_fixed(xlim = ax_lim, ylim = ax_lim) + theme_qc +
    ggplot2::labs(title = "Reported vs predicted (Horvath) age",
      subtitle = sprintf("r = %.3f, MAE = %.1f yrs, %d age outlier(s) (n = %d)", r, mae, length(outlier_ids), nrow(age_df)),
      x = "Reported age", y = "Predicted age (Horvath)"))
  if (!is.null(logger)) logger$log("qc_plots", sprintf("age check: r=%.3f, MAE=%.1f, %d outlier(s)", r, mae, length(outlier_ids)))
  list(outlier_ids = outlier_ids, detail = detail)
}

#' MDS with Donor labels
#' @keywords internal
#' @noRd
plot_mds_all <- function(betas_all_sub, ss_all, donor_col, sex_probes, n_top, theme_qc, logger) {
  keep <- grepl("^(cg|ch)", rownames(betas_all_sub))
  if (length(sex_probes) > 0) keep <- keep & !(rownames(betas_all_sub) %in% sex_probes)
  ba <- betas_all_sub[keep, , drop = FALSE]; ba <- ba[stats::complete.cases(ba), , drop = FALSE]
  rm(keep); gc(verbose = FALSE)
  if (nrow(ba) < 100 || ncol(ba) < 3) { if (!is.null(logger)) logger$log("qc_plots", "too few for MDS"); return(character(0)) }
  vars <- matrixStats::rowVars(ba); top <- order(vars, decreasing = TRUE)[seq_len(min(n_top, length(vars)))]
  bt <- ba[top, , drop = FALSE]; rm(vars, ba); gc(verbose = FALSE)
  d <- stats::dist(t(bt)); mds <- stats::cmdscale(d, k = 2); rm(d, bt); gc(verbose = FALSE)
  colnames(mds) <- c("MDS1", "MDS2")
  mds_df <- data.frame(MDS1 = mds[, 1], MDS2 = mds[, 2], sample_id = rownames(mds), stringsAsFactors = FALSE)
  if (!is.na(donor_col)) {
    m <- match(mds_df$sample_id, ss_all$sample_id)
    mds_df$label <- ss_all[[donor_col]][m]; mds_df$label[is.na(mds_df$label)] <- mds_df$sample_id[is.na(mds_df$label)]
  } else mds_df$label <- mds_df$sample_id
  medoid <- geometric_median(as.matrix(mds_df[, c("MDS1", "MDS2")]))
  all_dists <- sqrt((mds_df$MDS1 - medoid[1])^2 + (mds_df$MDS2 - medoid[2])^2)
  r4 <- 4 * stats::sd(all_dists); mds_df$outlier <- all_dists > r4
  outlier_ids <- mds_df$sample_id[mds_df$outlier]
  p <- ggplot2::ggplot(mds_df, ggplot2::aes(MDS1, MDS2, label = label, color = outlier)) +
    ggplot2::geom_text(size = 5, alpha = 0.8, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red3")) +
    ggplot2::annotate("point", x = medoid[1], y = medoid[2], color = "black", shape = 4, size = 4, stroke = 1.2)
  if (!is.na(r4) && r4 > 0) {
    th <- seq(0, 2 * pi, length.out = 200)
    circ <- data.frame(MDS1 = medoid[1] + r4 * cos(th), MDS2 = medoid[2] + r4 * sin(th))
    p <- p + ggplot2::geom_path(data = circ, ggplot2::aes(MDS1, MDS2), inherit.aes = FALSE, linetype = "dashed", color = "grey40")
  }
  p <- p + theme_qc + ggplot2::labs(title = "MDS (all samples, autosomal cg/ch probes)",
    subtitle = sprintf("%d samples; %d outlier(s) beyond 4 SD (red)", nrow(mds_df), length(outlier_ids)))
  print(p)
  if (!is.null(logger)) logger$log("qc_plots", sprintf("MDS: %d samples, %d outlier(s)", nrow(mds_df), length(outlier_ids)))
  rm(mds, mds_df); gc(verbose = FALSE); outlier_ids
}

#' PCA on OK samples
#' @keywords internal
#' @noRd
plot_pca_retained <- function(betas_qc, ss_all, excluded_ids, sex_probes, n_top, pc_csv_path, theme_qc, logger) {
  keep <- grepl("^(cg|ch)", rownames(betas_qc))
  if (length(sex_probes) > 0) keep <- keep & !(rownames(betas_qc) %in% sex_probes)
  bc <- betas_qc[keep, , drop = FALSE]; bc <- bc[stats::complete.cases(bc), , drop = FALSE]; rm(keep); gc(verbose = FALSE)
  if (nrow(bc) < 100 || ncol(bc) < 3) { if (!is.null(logger)) logger$log("qc_plots", "too few for PCA"); return(invisible(NULL)) }
  vars <- matrixStats::rowVars(bc); top <- order(vars, decreasing = TRUE)[seq_len(min(n_top, length(vars)))]
  bt <- bc[top, , drop = FALSE]; rm(vars, bc); gc(verbose = FALSE)
  ns <- ncol(bt); if (ns < 3) return(invisible(NULL))
  pca <- stats::prcomp(t(bt), scale. = TRUE); rm(bt); gc(verbose = FALSE)
  npc <- ncol(pca$x); pcs <- as.data.frame(pca$x); pcs$sample_id <- rownames(pca$x)
  utils::write.csv(pcs, pc_csv_path, row.names = FALSE)
  if (!is.null(logger)) logger$log("qc_plots", sprintf("wrote %d PC scores to %s", npc, pc_csv_path))
  ve <- (pca$sdev^2 / sum(pca$sdev^2)) * 100; nsc <- min(20, length(ve))
  print(ggplot2::ggplot(data.frame(PC = factor(paste0("PC", 1:nsc), levels = paste0("PC", 1:nsc)), vp = ve[1:nsc]),
    ggplot2::aes(PC, vp)) + ggplot2::geom_col(fill = "steelblue") + theme_qc +
    ggplot2::labs(title = "Scree plot (OK samples, autosomal cg/ch)", subtitle = sprintf("Top %d of %d PCs (%d samples)", nsc, npc, ns),
      y = "% variance explained", x = NULL))
  npp <- min(6L, npc); m <- match(pcs$sample_id, ss_all$sample_id); ss <- ss_all[m, ]
  qp <- "^frac_|^n_|^num_|^mean_|^na_|^med[RG]$|^top[RG]$|^InfI_|^RG|^sex_chr"
  skip <- c("sample_id","Basename","batch_folder","sheet_path","horvath_age","detected_platform",
    "Fname","Grn","Red","excluded","status","mds_status", grep(qp, colnames(ss), value = TRUE))
  cands <- setdiff(colnames(ss), skip); is_cont <- vapply(cands, function(c) is.numeric(ss[[c]]), logical(1))

  # Build all PC scatter plots, then arrange 3 per page
  pc_plots <- list()
  for (pi in seq_len(npp)) {
    pv <- pca$x[, pi]; partner <- if (pi == 1) 2L else 1L; px <- min(pi, partner); py <- max(pi, partner)
    bv <- NA_character_; bs <- Inf; bl <- NA_character_
    for (col in cands) {
      v <- ss[[col]]; if (all(is.na(v))) next
      if (is_cont[col]) { r <- suppressWarnings(stats::cor(pv, v, use = "complete.obs"))
        if (!is.na(r) && -abs(r) < bs) { bs <- -abs(r); bv <- col; bl <- sprintf("%s (|r|=%.3f)", col, abs(r)) }
      } else { f <- factor(v); nl <- nlevels(droplevels(f[!is.na(f)])); if (nl < 2 || nl > 20) next
        fit <- tryCatch(stats::lm(pv ~ f), error = function(e) NULL); if (is.null(fit)) next
        fs <- summary(fit)$fstatistic; if (is.null(fs)) next
        pval <- stats::pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
        if (!is.na(pval) && log10(pval) < bs && pval < 0.05) { bs <- log10(pval); bv <- col; bl <- sprintf("%s (ANOVA p=%.2e)", col, pval) } }
    }
    df <- data.frame(x = pca$x[, px], y = pca$x[, py])
    theme_pc <- ggplot2::theme_minimal() + ggplot2::theme(
      axis.title = ggplot2::element_text(size = 10), axis.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 8), legend.title = ggplot2::element_text(size = 9),
      legend.key.size = ggplot2::unit(0.4, "cm"))
    if (!is.na(bv)) { cv <- ss[[bv]]; df$color <- if (is.numeric(cv)) cv else factor(cv)
      p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = color)) + ggplot2::geom_point(size = 1.5, alpha = 0.7) + theme_pc +
        ggplot2::labs(title = sprintf("PC%d vs PC%d - %s", px, py, bl), x = sprintf("PC%d (%.1f%%)", px, ve[px]),
          y = sprintf("PC%d (%.1f%%)", py, ve[py]), color = bv)
      if (is.numeric(cv)) p <- p + ggplot2::scale_color_viridis_c()
    } else { p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_point(size = 1.5, alpha = 0.7, color = "steelblue") +
      theme_pc + ggplot2::labs(title = sprintf("PC%d vs PC%d (no significant variable)", px, py),
        x = sprintf("PC%d (%.1f%%)", px, ve[px]), y = sprintf("PC%d (%.1f%%)", py, ve[py])) }
    pc_plots[[pi]] <- p
    if (!is.null(logger)) logger$log("qc_plots", sprintf("PC%d -> %s", pi, if (is.na(bv)) "none" else bl))
  }

  # Arrange 3 plots per page using gridExtra
  for (page_start in seq(1, length(pc_plots), by = 3)) {
    page_end <- min(page_start + 2, length(pc_plots))
    page_plots <- pc_plots[page_start:page_end]
    # Pad with nullGrob if fewer than 3 plots on last page
    while (length(page_plots) < 3) page_plots[[length(page_plots) + 1]] <- ggplot2::ggplot() + ggplot2::theme_void()
    gridExtra::grid.arrange(grobs = page_plots, ncol = 1, nrow = 3)
  }
  rm(pca, pcs, ss, pc_plots); gc(verbose = FALSE); invisible(NULL)
}

#' Geometric median
#' @keywords internal
#' @noRd
geometric_median <- function(X, eps = 1e-5, max_iter = 200) {
  y <- colMeans(X)
  for (i in seq_len(max_iter)) { d <- sqrt(rowSums((X - matrix(y, nrow(X), ncol(X), byrow = TRUE))^2))
    d[d < eps] <- eps; w <- 1/d; y_new <- colSums(X * w) / sum(w)
    if (sqrt(sum((y_new - y)^2)) < eps) break; y <- y_new }; y
}
