###############################################################################
# plots.R - QC report PDF
#
# Single entry point: qcreport(). Writes a 9-page diagnostic PDF and two
# accompanying CSVs (PC scores, tail-probe failure rates).
#
# Page order (REORDERED in v2.0.0):
#   1. Detection rate per sample           (col by low-detection flag)
#   2. Per-probe failure histogram          (full 0-1 range; Y axis capped)
#   3. MDS on ALL non-sex cg/ch probes     (col by MDS outlier)
#   4. Mean intensity per sample           (MDS outlier > low-intensity > OK)
#   5. Beta density per sample             (NEW; col by low-intensity flag)
#   6. Sex check (chrX vs chrY intensity)  (col by reported sex)
#   7. Reported vs Horvath age             (col by sex mismatch)
#   8. Scree plot
#   9. PC vs associated variable (2x3)     (REPLACES PC-vs-PC scatters)
#
# Colour-coding rules (per design v2.0.0):
#   detection rate  - low-detection flag only
#   probe failure   - single colour; Y axis capped so the tail is visible.
#                     Dashed vertical line at failmin marks the CSV cutoff
#   MDS             - sample outside geometric-median 4-SD ring
#   intensity       - MDS outlier wins; else low-intensity flag; else OK
#   beta density    - low-intensity flag only
#   sex check       - reported sex (F / M / n/a); regression band per cluster
#   age check       - sex-mismatch flag (no "unclear")
#   scree           - none
#   PC vs variable  - by the variable's levels/values
###############################################################################

#' Write the QC report PDF and accompanying CSVs
#'
#' @param ss Sample sheet with QC metrics merged (rows = samples, must
#'   include \code{sample_id}, \code{frac_dt}, and -- where available --
#'   \code{mean_intensity}, \code{sex_chrX_intensity}, \code{sex_chrY_intensity},
#'   \code{horvath_age}).
#' @param betas Unmasked beta matrix (probes x samples) used for the
#'   beta-density panel and for MDS.
#' @param betasok Quality-masked / detection-masked beta matrix (probes
#'   x samples) used for the scree and PCA panels.
#' @param flagged Output of internal \code{qcflags()} -- data.frame of
#'   per-sample low-detection / low-intensity flags. Used ONLY for
#'   colouring the detection-rate, intensity, and beta-density panels.
#' @param mask Quality + detection mask matrix (probes x samples,
#'   logical). May be NULL.
#' @param detP Detection p-value matrix (probes x samples). May be NULL;
#'   if so, the probe-failure histogram falls back to \code{mask} or
#'   to \code{is.na(betas)}.
#' @param platform Platform string (e.g. "EPIC"), used to look up sex probes.
#' @param pdf Output PDF path.
#' @param pccsv Output CSV path for PC scores.
#' @param tailcsv Output CSV path for tail-probe failure rates.
#' @param logger Optional logger.
#' @return Invisibly, a list with sex_mismatch_ids, sex_unclear_ids,
#'   age_outlier_ids, sex_detail, age_detail.
#' @export
qcreport <- function(ss, betas, betasok, flagged,
                     mask, detP, platform,
                     pdf, pccsv, tailcsv,
                     logger = NULL) {
  cfg <- mqcopts()

  theme_qc <- ggplot2::theme_minimal() + ggplot2::theme(
    axis.title   = ggplot2::element_text(size = 14),
    axis.text    = ggplot2::element_text(size = 12),
    plot.title   = ggplot2::element_text(size = 15, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 12),
    legend.text  = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 13))

  grDevices::pdf(pdf, width = 10, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Derive flag-ID vectors from qcflags() output. Sex mismatch is NOT here.
  low_det_ids <- character(0)
  low_int_ids <- character(0)
  if (!is.null(flagged) && is.data.frame(flagged) && nrow(flagged) > 0) {
    low_det_ids <- flagged$sample_id[grepl("low_detection", flagged$reason)]
    low_int_ids <- flagged$sample_id[grepl("low_intensity", flagged$reason)]
  }
  ss$low_det <- ss$sample_id %in% low_det_ids
  ss$low_int <- ss$sample_id %in% low_int_ids

  sex_probes <- sexprobes(platform)
  donor_col  <- resolve_column(colnames(ss), cfg$donorcol, cfg$donoraliases)

  # ---- Page 1: Detection rate per sample (col by low-detection flag) ----
  ss$status_det <- ifelse(ss$low_det, "Low detection", "OK")
  theme_det <- theme_qc + ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 7))
  print(ggplot2::ggplot(ss, ggplot2::aes(
        x = stats::reorder(sample_id, frac_dt),
        y = frac_dt, fill = status_det)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("OK" = "steelblue", "Low detection" = "red3"),
      name = NULL) +
    ggplot2::geom_hline(yintercept = cfg$samplemin, linetype = "dashed",
                        color = "black", linewidth = 0.7) +
    ggplot2::coord_flip() + theme_det +
    ggplot2::labs(
      title = "Detection rate per sample (frac_dt)",
      subtitle = sprintf(
        "All %d samples; red = below %.2f threshold; dashed line = threshold",
        nrow(ss), cfg$samplemin),
      x = NULL, y = "Fraction detected"))

  # ---- Page 2: Per-probe failure histogram (full range, Y capped) ----
  mds_outlier_ids <- character(0)  # filled in by plot_mds()
  plot_probe_failure(betas, mask, detP, cfg, theme_qc, tailcsv, logger)

  # ---- Page 3: MDS on ALL non-sex cg/ch probes ----
  mds_outlier_ids <- plot_mds(betas, ss, donor_col, sex_probes, theme_qc, logger)

  # ---- Page 4: Mean intensity (MDS outlier > low-intensity > OK) ----
  if ("mean_intensity" %in% colnames(ss)) {
    ss$intensity_color <- ifelse(
      ss$sample_id %in% mds_outlier_ids, "MDS outlier",
      ifelse(ss$low_int, "Low intensity", "OK"))
    ss$intensity_color <- factor(
      ss$intensity_color, levels = c("OK", "Low intensity", "MDS outlier"))
    print(ggplot2::ggplot(ss, ggplot2::aes(mean_intensity,
                                           fill = intensity_color)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.85,
                              position = "identity") +
      ggplot2::scale_fill_manual(
        values = c("OK" = "steelblue",
                   "Low intensity" = "goldenrod3",
                   "MDS outlier" = "red3"),
        name = NULL, drop = FALSE) +
      ggplot2::geom_vline(xintercept = cfg$intmin, linetype = "dashed",
                          color = "black", linewidth = 0.7) +
      theme_qc + ggplot2::labs(
        title = "Mean intensity per sample",
        subtitle = sprintf(
          paste0("All %d samples; MDS-outlier colour takes precedence ",
                 "over low-intensity"),
          nrow(ss)),
        x = "Mean intensity", y = "Count"))
  }

  # ---- Page 5: Sample beta density (col by low-intensity flag) ----
  plot_beta_density(betas, ss, theme_qc, logger)

  # ---- Page 6: Sex check ----
  sex_result <- plot_sex_check(ss, theme_qc, logger)

  # ---- Page 7: Age check ----
  age_result <- plot_age_check(ss, sex_mismatch_ids = sex_result$mismatch_ids,
                               theme_qc, logger)

  # ---- Page 8 + 9: Scree plot + PC-vs-variable 2x3 panel ----
  plot_pca_retained(betasok, ss, sex_probes, pccsv, theme_qc, logger)

  invisible(list(sex_mismatch_ids = sex_result$mismatch_ids,
                 sex_unclear_ids  = sex_result$unclear_ids,
                 age_outlier_ids  = age_result$outlier_ids,
                 sex_detail       = sex_result$detail,
                 age_detail       = age_result$detail))
}

###############################################################################
# Per-probe sample-failure histogram (Page 2)
###############################################################################
#' Per-probe sample-failure histogram (full range, Y axis capped)
#'
#' Plots a histogram of per-probe sample-failure rate across the full
#' 0-1 range. The Y axis is capped so the long tail is visible, at the
#' cost of clipping the leftmost spike (probes that fail in near-0% of
#' samples). A dashed vertical line marks \code{failmin}, the threshold
#' above which probes are written to \code{tailcsv}.
#'
#' Failure rate is derived from detP by default
#' (\code{rowMeans(detP > pthresh)}), then mask, then \code{is.na(betas)}.
#'
#' Quality-masked probes are excluded by default; pass \code{inclqual = TRUE}
#' via \code{mqcopts()} to keep them.
#' @keywords internal
#' @noRd
plot_probe_failure <- function(betas, mask, detP, cfg, theme_qc,
                               tailcsv, logger = NULL) {
  pthresh  <- cfg$detp
  failmin  <- cfg$failmin
  inclqual <- isTRUE(cfg$inclqual)

  # Per-probe failure rate -------------------------------------------------
  source <- "detP"
  if (!is.null(detP) && is.matrix(detP)) {
    fail_rate <- rowMeans(detP > pthresh, na.rm = TRUE)
  } else if (!is.null(mask) && is.matrix(mask)) {
    fail_rate <- rowMeans(mask, na.rm = TRUE); source <- "mask"
  } else {
    fail_rate <- rowMeans(is.na(betas));        source <- "is.na(betas)"
  }
  probe_ids <- if (!is.null(detP) && !is.null(rownames(detP))) rownames(detP) else
               if (!is.null(mask) && !is.null(rownames(mask))) rownames(mask) else
                                                                rownames(betas)

  # Identify quality-masked probes ----------------------------------------
  is_qual_masked <- rep(FALSE, length(fail_rate))
  if (!inclqual && !is.null(mask) && is.matrix(mask)) {
    if (!is.null(detP) && is.matrix(detP)) {
      # Quality mask catches probes flagged TRUE in mask even when detP <= pthresh
      qual_frac <- rowMeans(mask & (detP <= pthresh), na.rm = TRUE)
    } else {
      qual_frac <- rowMeans(mask, na.rm = TRUE)
    }
    is_qual_masked <- qual_frac > 0.5
    if (!is.null(logger))
      logger$log("qc_plots",
                 sprintf("probe-failure: excluding %d quality-masked probe(s)",
                         sum(is_qual_masked)))
  }

  keep <- !is_qual_masked
  fail_rate <- fail_rate[keep]
  probe_ids <- probe_ids[keep]

  n_total  <- length(fail_rate)
  tail_idx <- which(fail_rate >= failmin)
  n_tail   <- length(tail_idx)

  # Write the tail probes to CSV (probes with fail_rate >= failmin) -------
  tail_df <- data.frame(
    probe_id  = probe_ids[tail_idx],
    fail_rate = fail_rate[tail_idx],
    stringsAsFactors = FALSE)
  tail_df <- tail_df[order(-tail_df$fail_rate), ]
  utils::write.csv(tail_df, tailcsv, row.names = FALSE)
  if (!is.null(logger))
    logger$log("qc_plots",
               sprintf("probe-failure: %d probes with fail_rate >= %.2f written to %s (source=%s)",
                       n_tail, failmin, tailcsv, source))

  # Bin counts manually so the Y cap can be derived from the data ---------
  nbins  <- 40L
  breaks <- seq(0, 1, length.out = nbins + 1L)
  hh     <- graphics::hist(fail_rate, breaks = breaks, plot = FALSE)
  counts <- hh$counts
  mids   <- hh$mids

  # Y cap: 1.5x the tallest bin centred at fail_rate >= 0.05. The leftmost
  # spike of probes passing in ~all samples is deliberately ignored from
  # the cap; that bar still extends above the visible area. The 0.05
  # cap-exclusion cutoff is independent of failmin: it's about clipping
  # the leftmost spike, not the CSV threshold.
  rest   <- counts[mids >= 0.05]
  y_cap  <- if (length(rest) > 0 && max(rest) > 0) 1.5 * max(rest)
            else                                    max(counts)
  if (!is.finite(y_cap) || y_cap <= 0) y_cap <- max(counts, 1)

  spike_bin <- counts[1]
  spike_clipped <- spike_bin > y_cap

  qmask_note <- if (inclqual) " (incl. quality-masked)" else ""
  subtitle_txt <- sprintf(
    "%s probes from %s%s. %s with fail >= %.2f written to %s.%s",
    format(n_total, big.mark = ","),
    source,
    qmask_note,
    format(n_tail,  big.mark = ","),
    failmin,
    basename(tailcsv),
    if (spike_clipped)
      sprintf("\nLeftmost bar clipped: %s probes at fail_rate ~ 0 (y capped at %s).",
              format(spike_bin, big.mark = ","),
              format(round(y_cap), big.mark = ","))
    else "")

  print(ggplot2::ggplot(
        data.frame(fr = fail_rate),
        ggplot2::aes(x = fr)) +
    ggplot2::geom_histogram(breaks = breaks,
                            fill = "steelblue", color = "white") +
    ggplot2::geom_vline(xintercept = failmin, linetype = "dashed",
                        color = "red3", linewidth = 0.7) +
    ggplot2::annotate("text", x = failmin, y = y_cap,
                      label = sprintf(" failmin = %.2f", failmin),
                      hjust = 0, vjust = 1.1,
                      color = "red3", size = 4) +
    ggplot2::scale_x_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, by = 0.1)) +
    ggplot2::coord_cartesian(ylim = c(0, y_cap)) +
    theme_qc + ggplot2::labs(
      title = "Per-probe sample-failure rate",
      subtitle = subtitle_txt,
      x = sprintf("Fraction of samples where probe failed (detP > %.2f)",
                  pthresh),
      y = "Number of probes"))
  invisible(NULL)
}

###############################################################################
# MDS (Page 3) - uses ALL non-sex cg/ch probes
###############################################################################
#' MDS panel, using all complete-case non-sex cg/ch probes
#' @keywords internal
#' @noRd
plot_mds <- function(betas, ss, donor_col, sex_probes, theme_qc, logger) {
  keep <- grepl("^(cg|ch)", rownames(betas))
  if (length(sex_probes) > 0)
    keep <- keep & !(rownames(betas) %in% sex_probes)
  ba <- betas[keep, , drop = FALSE]
  ba <- ba[stats::complete.cases(ba), , drop = FALSE]
  rm(keep); gc(verbose = FALSE)

  if (nrow(ba) < 100 || ncol(ba) < 3) {
    if (!is.null(logger))
      logger$log("qc_plots", "too few probes/samples for MDS")
    return(character(0))
  }

  if (!is.null(logger))
    logger$log("qc_plots",
               sprintf("MDS using ALL %s complete-case cg/ch probes",
                       format(nrow(ba), big.mark = ",")))

  d <- stats::dist(t(ba))
  mds <- stats::cmdscale(d, k = 2)
  rm(d, ba); gc(verbose = FALSE)
  colnames(mds) <- c("MDS1", "MDS2")

  mds_df <- data.frame(MDS1 = mds[, 1], MDS2 = mds[, 2],
                       sample_id = rownames(mds),
                       stringsAsFactors = FALSE)

  if (!is.na(donor_col)) {
    m <- match(mds_df$sample_id, ss$sample_id)
    mds_df$label <- ss[[donor_col]][m]
    mds_df$label[is.na(mds_df$label)] <- mds_df$sample_id[is.na(mds_df$label)]
  } else {
    mds_df$label <- mds_df$sample_id
  }

  medoid <- geometric_median(as.matrix(mds_df[, c("MDS1", "MDS2")]))
  all_dists <- sqrt((mds_df$MDS1 - medoid[1])^2 +
                    (mds_df$MDS2 - medoid[2])^2)
  r4 <- 4 * stats::sd(all_dists)
  mds_df$outlier <- all_dists > r4
  outlier_ids <- mds_df$sample_id[mds_df$outlier]

  p <- ggplot2::ggplot(mds_df, ggplot2::aes(MDS1, MDS2, label = label,
                                            color = outlier)) +
    ggplot2::geom_text(size = 5, alpha = 0.8, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red3")) +
    ggplot2::annotate("point", x = medoid[1], y = medoid[2],
                      color = "black", shape = 4, size = 4, stroke = 1.2)

  if (!is.na(r4) && r4 > 0) {
    th <- seq(0, 2 * pi, length.out = 200)
    circ <- data.frame(MDS1 = medoid[1] + r4 * cos(th),
                       MDS2 = medoid[2] + r4 * sin(th))
    p <- p + ggplot2::geom_path(data = circ, ggplot2::aes(MDS1, MDS2),
                                inherit.aes = FALSE,
                                linetype = "dashed", color = "grey40")
  }

  p <- p + theme_qc + ggplot2::labs(
    title = "MDS (all complete-case cg/ch probes; sex chromosomes excluded)",
    subtitle = sprintf("%d samples; %d outlier(s) beyond 4 SD ring (red)",
                       nrow(mds_df), length(outlier_ids)))
  print(p)

  if (!is.null(logger))
    logger$log("qc_plots",
               sprintf("MDS: %d samples, %d outlier(s)",
                       nrow(mds_df), length(outlier_ids)))
  rm(mds, mds_df); gc(verbose = FALSE)
  outlier_ids
}

###############################################################################
# Beta density (Page 5, NEW)
###############################################################################
#' Per-sample beta-value density, overlaid on one panel
#'
#' One density curve per sample, computed from noob-corrected betas on
#' autosomal cg/ch probes. Lines are coloured by the low-intensity
#' flag (red3) vs OK (steelblue). No legend, per design.
#' @keywords internal
#' @noRd
plot_beta_density <- function(betas, ss, theme_qc, logger = NULL) {
  cg_ch <- grepl("^(cg|ch)", rownames(betas))
  bs <- betas[cg_ch, , drop = FALSE]
  if (nrow(bs) < 100 || ncol(bs) < 1) {
    if (!is.null(logger))
      logger$log("qc_plots", "beta-density skipped: too few cg/ch probes")
    return(invisible(NULL))
  }

  n_grid <- 512L
  sample_ids <- colnames(bs)
  parts <- lapply(seq_along(sample_ids), function(j) {
    v <- bs[, j]
    v <- v[!is.na(v) & v >= 0 & v <= 1]
    if (length(v) < 50) return(NULL)
    d <- stats::density(v, from = 0, to = 1, n = n_grid,
                        na.rm = TRUE)
    data.frame(beta = d$x, density = d$y,
               sample_id = sample_ids[j],
               stringsAsFactors = FALSE)
  })
  dens_df <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
  if (is.null(dens_df) || nrow(dens_df) == 0) {
    if (!is.null(logger))
      logger$log("qc_plots", "beta-density skipped: no valid samples")
    return(invisible(NULL))
  }

  dens_df$low_int <- dens_df$sample_id %in%
                     ss$sample_id[isTRUE_vec(ss$low_int)]

  n_low <- length(unique(dens_df$sample_id[dens_df$low_int]))
  n_ok  <- length(unique(dens_df$sample_id)) - n_low

  print(ggplot2::ggplot(dens_df, ggplot2::aes(x = beta, y = density,
                                              group = sample_id,
                                              color = low_int)) +
    ggplot2::geom_line(alpha = 0.4, linewidth = 0.3) +
    ggplot2::scale_color_manual(values = c("FALSE" = "steelblue",
                                           "TRUE"  = "red3")) +
    ggplot2::scale_x_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, 0.2)) +
    theme_qc + ggplot2::theme(legend.position = "none") +
    ggplot2::labs(
      title = "Sample beta-value density (noob-corrected)",
      subtitle = sprintf(
        "One curve per sample; %d OK (blue), %d low-intensity (red); cg/ch probes",
        n_ok, n_low),
      x = "Beta value", y = "Density"))

  if (!is.null(logger))
    logger$log("qc_plots",
               sprintf("beta-density: %d sample curve(s) drawn (%d low-intensity)",
                       length(unique(dens_df$sample_id)), n_low))
  rm(bs, dens_df, parts); gc(verbose = FALSE)
  invisible(NULL)
}

# Coerces NA -> FALSE in a logical vector
isTRUE_vec <- function(x) !is.na(x) & as.logical(x)

###############################################################################
# Sex check (Page 6)
###############################################################################
#' Sex check via chrX vs chrY median intensity, with regression bands
#' @keywords internal
#' @noRd
plot_sex_check <- function(ss, theme_qc, logger = NULL) {
  cfg <- mqcopts()
  empty <- list(mismatch_ids = character(0), unclear_ids = character(0),
                detail = data.frame(sample_id = character(0),
                                     reported_sex = character(0),
                                     inferred_sex_intensity = character(0),
                                     stringsAsFactors = FALSE))
  has_x <- "sex_chrX_intensity" %in% colnames(ss)
  has_y <- "sex_chrY_intensity" %in% colnames(ss)
  if (!has_x || !has_y) return(empty)

  sex_df <- data.frame(sample_id = ss$sample_id,
                       chrX = ss$sex_chrX_intensity,
                       chrY = ss$sex_chrY_intensity,
                       stringsAsFactors = FALSE)
  sex_df <- sex_df[!is.na(sex_df$chrX) & !is.na(sex_df$chrY), ]
  if (nrow(sex_df) < 6) return(empty)

  m <- match(sex_df$sample_id, ss$sample_id)
  sex_col <- resolve_column(colnames(ss), cfg$sexcol, cfg$sexaliases)
  sex_df$reported_sex <- if (!is.na(sex_col))
    normalize_sex(ss[[sex_col]][m], col_name = sex_col) else NA_character_
  n_missing_sex <- sum(is.na(sex_df$reported_sex))

  # Optimise chrY threshold to minimise within-cluster regression residuals
  candidates <- stats::quantile(sex_df$chrY,
                                probs = seq(0.15, 0.85, by = 0.01),
                                na.rm = TRUE)
  best_thresh <- stats::median(sex_df$chrY); best_cost <- Inf
  for (thr in candidates) {
    fi <- which(sex_df$chrY <= thr); mi <- which(sex_df$chrY > thr)
    if (length(fi) < 3 || length(mi) < 3) next
    ff <- tryCatch(stats::lm(chrY ~ chrX, data = sex_df[fi, ]),
                   error = function(e) NULL)
    fm <- tryCatch(stats::lm(chrY ~ chrX, data = sex_df[mi, ]),
                   error = function(e) NULL)
    if (is.null(ff) || is.null(fm)) next
    cost <- sum(abs(stats::residuals(ff))) + sum(abs(stats::residuals(fm)))
    if (cost < best_cost) { best_cost <- cost; best_thresh <- thr }
  }
  f_idx <- which(sex_df$chrY <= best_thresh)
  m_idx <- which(sex_df$chrY >  best_thresh)
  fit_f <- if (length(f_idx) >= 3)
    tryCatch(stats::lm(chrY ~ chrX, data = sex_df[f_idx, ]),
             error = function(e) NULL) else NULL
  fit_m <- if (length(m_idx) >= 3)
    tryCatch(stats::lm(chrY ~ chrX, data = sex_df[m_idx, ]),
             error = function(e) NULL) else NULL

  ortho_dist <- function(x, y, fit) {
    if (is.null(fit)) return(rep(Inf, length(x)))
    a <- stats::coef(fit)[2]; c <- stats::coef(fit)[1]
    abs(a * x - y + c) / sqrt(a^2 + 1)
  }
  d_f <- ortho_dist(sex_df$chrX, sex_df$chrY, fit_f)
  d_m <- ortho_dist(sex_df$chrX, sex_df$chrY, fit_m)
  sd_f <- if (length(f_idx) >= 3) stats::sd(d_f[f_idx]) else Inf
  sd_m <- if (length(m_idx) >= 3) stats::sd(d_m[m_idx]) else Inf
  band_sd <- 5.0
  thresh_f <- band_sd * sd_f; thresh_m <- band_sd * sd_m

  sex_df$inferred_sex_intensity <- ifelse(
    d_f <= thresh_f & d_m >  thresh_m, "F",
    ifelse(d_m <= thresh_m & d_f >  thresh_f, "M",
    ifelse(d_f <= thresh_f & d_m <= thresh_m,
           ifelse(sex_df$chrY <= best_thresh, "F", "M"), "Unclear")))

  sex_df$mismatch <- !is.na(sex_df$reported_sex) &
                     sex_df$inferred_sex_intensity != "Unclear" &
                     sex_df$reported_sex != sex_df$inferred_sex_intensity
  sex_df$unclear <- sex_df$inferred_sex_intensity == "Unclear"
  mismatch_ids <- sex_df$sample_id[sex_df$mismatch]
  unclear_ids  <- sex_df$sample_id[sex_df$unclear]

  if (!is.null(logger)) {
    logger$log("qc_plots", sprintf(
      "sex clustering: threshold=%.0f, F=%d, M=%d, mismatch=%d, unclear=%d, missing_sex=%d",
      best_thresh, length(f_idx), length(m_idx),
      length(mismatch_ids), length(unclear_ids), n_missing_sex))
    if (length(mismatch_ids) > 0)
      logger$log("qc_plots", sprintf("sex mismatch IDs: %s",
                                     paste(mismatch_ids, collapse = ", ")))
    if (length(unclear_ids) > 0)
      logger$log("qc_plots", sprintf("sex unclear IDs: %s",
                                     paste(unclear_ids, collapse = ", ")))
  }
  detail <- sex_df[, c("sample_id", "reported_sex",
                       "inferred_sex_intensity"), drop = FALSE]

  sex_df$color_var <- ifelse(is.na(sex_df$reported_sex), "n/a",
                             sex_df$reported_sex)
  sex_df$color_var <- factor(sex_df$color_var, levels = c("F", "M", "n/a"))
  p <- ggplot2::ggplot(sex_df, ggplot2::aes(chrX, chrY, color = color_var,
                                            shape = color_var)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::scale_color_manual(
      values = c("F" = "#E41A1C", "M" = "#377EB8", "n/a" = "grey50"),
      name = "Reported sex") +
    ggplot2::scale_shape_manual(
      values = c("F" = 16, "M" = 17, "n/a" = 1),
      name = "Reported sex")

  x_range <- range(sex_df$chrX)
  x_seq <- seq(x_range[1] - diff(x_range) * 0.05,
               x_range[2] + diff(x_range) * 0.05, length.out = 200)
  for (cl in list(list(fit = fit_f, thresh = thresh_f, col = "#E41A1C"),
                  list(fit = fit_m, thresh = thresh_m, col = "#377EB8"))) {
    if (is.null(cl$fit) || !is.finite(cl$thresh)) next
    a <- stats::coef(cl$fit)[2]; b_int <- stats::coef(cl$fit)[1]
    y_pred <- b_int + a * x_seq
    y_off  <- cl$thresh / cos(atan(a))
    band_df <- data.frame(x = c(x_seq, rev(x_seq)),
                          y = c(y_pred + y_off, rev(y_pred - y_off)))
    p <- p + ggplot2::geom_polygon(data = band_df, ggplot2::aes(x, y),
                                   inherit.aes = FALSE, fill = cl$col,
                                   alpha = 0.20)
  }
  subtitle_txt <- sprintf(
    "All %d samples; %d mismatch(es), %d predicted sex unclear",
    nrow(sex_df), length(mismatch_ids), length(unclear_ids))
  if (n_missing_sex > 0)
    subtitle_txt <- paste0(subtitle_txt,
                           sprintf("; %d missing reported sex", n_missing_sex))
  p <- p + theme_qc + ggplot2::labs(
    title = "Sex check: chrX vs chrY median intensity (curated probes)",
    subtitle = subtitle_txt,
    x = "Median total intensity (chrX.xlinked)",
    y = "Median total intensity (chrY.clean)")
  print(p)

  list(mismatch_ids = mismatch_ids, unclear_ids = unclear_ids,
       detail = detail)
}

###############################################################################
# Age check (Page 7)
###############################################################################
#' Reported vs Horvath age, coloured by sex mismatch
#' @keywords internal
#' @noRd
plot_age_check <- function(ss, sex_mismatch_ids, theme_qc, logger = NULL) {
  cfg <- mqcopts()
  empty <- list(outlier_ids = character(0),
                detail = data.frame(sample_id = character(0),
                                    reported_age = numeric(0),
                                    horvath_age = numeric(0),
                                    ortho_dist = numeric(0),
                                    stringsAsFactors = FALSE))
  age_col <- resolve_column(colnames(ss), cfg$agecol, cfg$agealiases)
  if (is.na(age_col)) {
    if (!is.null(logger))
      logger$log("qc_plots", "age check skipped: no age column")
    return(empty)
  }
  if (!"horvath_age" %in% colnames(ss)) {
    if (!is.null(logger))
      logger$log("qc_plots", "age check skipped: no horvath_age")
    return(empty)
  }
  age_df <- data.frame(sample_id = ss$sample_id,
                       reported = parse_age_robust(ss[[age_col]]),
                       horvath  = ss$horvath_age,
                       stringsAsFactors = FALSE)
  age_df <- age_df[!is.na(age_df$reported) & !is.na(age_df$horvath), ]
  if (nrow(age_df) < 3) {
    if (!is.null(logger))
      logger$log("qc_plots",
                 sprintf("age check skipped: %d valid pairs", nrow(age_df)))
    return(empty)
  }
  fit <- stats::lm(horvath ~ reported, data = age_df)
  a <- stats::coef(fit)[2]; b <- stats::coef(fit)[1]
  age_df$ortho_dist <- abs(a * age_df$reported - age_df$horvath + b) /
                       sqrt(a^2 + 1)
  dist_sd <- stats::sd(age_df$ortho_dist)
  age_df$age_outlier <- age_df$ortho_dist > 3 * dist_sd
  outlier_ids <- age_df$sample_id[age_df$age_outlier]
  detail <- age_df[age_df$age_outlier,
                   c("sample_id", "reported", "horvath", "ortho_dist"),
                   drop = FALSE]
  colnames(detail) <- c("sample_id", "reported_age",
                        "horvath_age", "ortho_dist")
  r   <- suppressWarnings(stats::cor(age_df$reported, age_df$horvath))
  mae <- mean(abs(age_df$reported - age_df$horvath))

  age_df$color_var <- ifelse(age_df$sample_id %in% sex_mismatch_ids,
                             "Sex mismatch", "OK")
  age_df$color_var <- factor(age_df$color_var,
                             levels = c("OK", "Sex mismatch"))
  n_mm_in_age <- sum(age_df$sample_id %in% sex_mismatch_ids)

  all_vals <- c(age_df$reported, age_df$horvath)
  ax_lim   <- range(all_vals) + c(-1, 1) * diff(range(all_vals)) * 0.05
  x_seq    <- seq(ax_lim[1], ax_lim[2], length.out = 200)
  y_pred   <- b + a * x_seq
  y_off    <- 3 * dist_sd / cos(atan(a))
  band_df  <- data.frame(x = c(x_seq, rev(x_seq)),
                         y = c(y_pred + y_off, rev(y_pred - y_off)))

  print(ggplot2::ggplot(age_df, ggplot2::aes(reported, horvath,
                                             color = color_var)) +
    ggplot2::geom_polygon(data = band_df, ggplot2::aes(x, y),
                          inherit.aes = FALSE, fill = "steelblue",
                          alpha = 0.15) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey70") +
    ggplot2::geom_point(size = 2.5, alpha = 0.8) +
    ggplot2::scale_color_manual(
      values = c("OK" = "steelblue", "Sex mismatch" = "red3"),
      name = NULL) +
    ggplot2::coord_fixed(xlim = ax_lim, ylim = ax_lim) + theme_qc +
    ggplot2::labs(
      title = "Reported vs predicted (Horvath) age",
      subtitle = sprintf(
        "r = %.3f, MAE = %.1f yrs, %d age outlier(s) (n = %d); %d sex-mismatch sample(s) shown",
        r, mae, length(outlier_ids), nrow(age_df), n_mm_in_age),
      x = "Reported age", y = "Predicted age (Horvath)"))

  if (!is.null(logger))
    logger$log("qc_plots",
               sprintf("age check: r=%.3f, MAE=%.1f, %d outlier(s)",
                       r, mae, length(outlier_ids)))
  list(outlier_ids = outlier_ids, detail = detail)
}

###############################################################################
# Scree + PCA panel (Pages 8 + 9)
###############################################################################
#' Scree plot + a single 2x3 page of PC-vs-associated-variable panels
#' @keywords internal
#' @noRd
plot_pca_retained <- function(betasok, ss, sex_probes,
                              pccsv, theme_qc, logger) {
  cfg <- mqcopts()
  ntop <- cfg$ntop

  keep <- grepl("^(cg|ch)", rownames(betasok))
  if (length(sex_probes) > 0)
    keep <- keep & !(rownames(betasok) %in% sex_probes)
  bc <- betasok[keep, , drop = FALSE]
  bc <- bc[stats::complete.cases(bc), , drop = FALSE]
  rm(keep); gc(verbose = FALSE)

  if (nrow(bc) < 100 || ncol(bc) < 3) {
    if (!is.null(logger))
      logger$log("qc_plots", "too few probes/samples for PCA")
    return(invisible(NULL))
  }

  vars <- matrixStats::rowVars(bc)
  top  <- order(vars, decreasing = TRUE)[seq_len(min(ntop, length(vars)))]
  bt   <- bc[top, , drop = FALSE]
  rm(vars, bc); gc(verbose = FALSE)

  ns <- ncol(bt)
  if (ns < 3) return(invisible(NULL))
  pca <- stats::prcomp(t(bt), scale. = TRUE)
  rm(bt); gc(verbose = FALSE)
  npc <- ncol(pca$x)
  pcs <- as.data.frame(pca$x); pcs$sample_id <- rownames(pca$x)
  utils::write.csv(pcs, pccsv, row.names = FALSE)

  if (!is.null(logger))
    logger$log("qc_plots",
               sprintf("PCA: %d probes (top by variance), %d samples; wrote %s",
                       nrow(pca$rotation), ns, pccsv))

  ve <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  nsc <- min(20, length(ve))

  # ---- Page 8: scree (variance labels KEPT here per design) ----
  print(ggplot2::ggplot(
        data.frame(PC = factor(paste0("PC", seq_len(nsc)),
                               levels = paste0("PC", seq_len(nsc))),
                   vp = ve[seq_len(nsc)]),
        ggplot2::aes(PC, vp)) +
    ggplot2::geom_col(fill = "steelblue") + theme_qc +
    ggplot2::theme(
      axis.line  = ggplot2::element_line(color = "black", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.4)) +
    ggplot2::labs(
      title = "Scree plot (cg/ch probes, sex excluded)",
      subtitle = sprintf("Top %d of %d PCs (%d samples); ntop = %s probes",
                         nsc, npc, ns, format(nrow(pca$rotation),
                                              big.mark = ",")),
      y = "% variance explained", x = NULL))

  # ---- Page 9: PC vs associated variable, 2x3 panel on one page ----
  npp <- min(6L, npc)
  m   <- match(pcs$sample_id, ss$sample_id)
  ss2 <- ss[m, ]

  qp <- "^frac_|^n_|^num_|^mean_|^na_|^med[RG]$|^top[RG]$|^InfI_|^RG|^sex_chr"
  skip <- c("sample_id", "Basename", "batch_folder", "sheet_path",
            "horvath_age", "detected_platform",
            "Fname", "Grn", "Red", "low_det", "low_int", "status_det",
            "intensity_color",
            grep(qp, colnames(ss2), value = TRUE))
  cands  <- setdiff(colnames(ss2), skip)
  is_cont <- vapply(cands, function(c) is.numeric(ss2[[c]]), logical(1))

  panel_theme <- ggplot2::theme_minimal() + ggplot2::theme(
    axis.title    = ggplot2::element_text(size = 10),
    axis.text     = ggplot2::element_text(size = 8),
    axis.line     = ggplot2::element_line(color = "black", linewidth = 0.4),
    axis.ticks    = ggplot2::element_line(color = "black", linewidth = 0.4),
    plot.title    = ggplot2::element_text(size = 11, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 9),
    legend.text   = ggplot2::element_text(size = 8),
    legend.title  = ggplot2::element_text(size = 9),
    legend.key.size = ggplot2::unit(0.4, "cm"))

  pc_panels <- vector("list", npp)
  for (pi in seq_len(npp)) {
    pv <- pca$x[, pi]
    bv <- NA_character_; bs <- Inf; bl <- NA_character_; b_kind <- NA_character_
    for (col in cands) {
      v <- ss2[[col]]
      if (all(is.na(v))) next
      if (is_cont[col]) {
        r <- suppressWarnings(stats::cor(pv, v, use = "complete.obs"))
        if (!is.na(r) && -abs(r) < bs) {
          bs <- -abs(r); bv <- col
          bl <- sprintf("|r|=%.3f", abs(r)); b_kind <- "continuous"
        }
      } else {
        f <- factor(v)
        nl <- nlevels(droplevels(f[!is.na(f)]))
        if (nl < 2 || nl > 20) next
        fit <- tryCatch(stats::lm(pv ~ f), error = function(e) NULL)
        if (is.null(fit)) next
        fs <- summary(fit)$fstatistic
        if (is.null(fs)) next
        pval <- stats::pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
        if (!is.na(pval) && log10(pval) < bs && pval < 0.05) {
          bs <- log10(pval); bv <- col
          bl <- sprintf("ANOVA p=%.2e", pval); b_kind <- "categorical"
        }
      }
    }

    if (is.na(bv)) {
      # No significant associated variable: PC vs sample index, no variance label
      df <- data.frame(x = seq_along(pv), y = pv)
      p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
        ggplot2::geom_point(size = 1.5, alpha = 0.7, color = "steelblue") +
        panel_theme + ggplot2::labs(
          title = sprintf("PC%d (no significant variable)", pi),
          x = "Sample index", y = sprintf("PC%d", pi))
    } else if (b_kind == "continuous") {
      df <- data.frame(x = ss2[[bv]], y = pv)
      p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
        ggplot2::geom_point(size = 1.5, alpha = 0.7, color = "steelblue") +
        panel_theme + ggplot2::labs(
          title    = sprintf("PC%d vs %s", pi, bv),
          subtitle = bl, x = bv, y = sprintf("PC%d", pi))
    } else {
      df <- data.frame(x = factor(ss2[[bv]]), y = pv)
      df <- df[!is.na(df$x), ]
      p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = x)) +
        ggplot2::geom_boxplot(outlier.shape = NA, color = "grey40",
                              fill = NA) +
        ggplot2::geom_jitter(width = 0.15, size = 1.3, alpha = 0.75) +
        ggplot2::guides(color = "none") +
        panel_theme + ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
        ggplot2::labs(title    = sprintf("PC%d by %s", pi, bv),
                      subtitle = bl, x = bv, y = sprintf("PC%d", pi))
    }
    pc_panels[[pi]] <- p

    if (!is.null(logger))
      logger$log("qc_plots",
                 sprintf("PC%d -> %s", pi,
                         if (is.na(bv)) "none" else sprintf("%s (%s)", bv, bl)))
  }

  # Pad to 6, arrange 2 rows x 3 cols, single page
  while (length(pc_panels) < 6)
    pc_panels[[length(pc_panels) + 1]] <-
      ggplot2::ggplot() + ggplot2::theme_void()
  gridExtra::grid.arrange(grobs = pc_panels, nrow = 2, ncol = 3)

  rm(pca, pcs, ss2, pc_panels); gc(verbose = FALSE)
  invisible(NULL)
}

#' Geometric median (Weiszfeld)
#' @keywords internal
#' @noRd
geometric_median <- function(X, eps = 1e-5, max_iter = 200) {
  y <- colMeans(X)
  for (i in seq_len(max_iter)) {
    d <- sqrt(rowSums((X - matrix(y, nrow(X), ncol(X), byrow = TRUE))^2))
    d[d < eps] <- eps; w <- 1 / d
    y_new <- colSums(X * w) / sum(w)
    if (sqrt(sum((y_new - y)^2)) < eps) break
    y <- y_new
  }
  y
}
