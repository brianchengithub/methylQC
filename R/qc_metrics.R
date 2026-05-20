###############################################################################
# qc_metrics.R — Per-sample QC metric computation
#
# Three exported entry points:
#   qcmetrics()  — from an in-memory list of SigDF objects
#   qcstream()   — streaming variant (reads IDATs one at a time)
#   checkmeta()  — cross-checks reported metadata vs computed metrics
#
# Both qcmetrics() and qcstream() return the SAME columns; downstream
# code in pipeline.R / plots.R must not branch on which path was used.
###############################################################################

#' Compute per-sample QC metrics from a list of SigDFs
#'
#' @param sdfs Named list of SeSAMe SigDF objects.
#' @param platform Platform string (unused; kept for symmetry with
#'   \code{qcstream}).
#' @param logger Optional logger.
#' @return data.frame with one row per sample.
#' @export
qcmetrics <- function(sdfs, platform, logger = NULL) {
  if (!is.null(logger))
    logger$log("qcmetrics", "computing per-sample QC stats...")

  qc_list <- lapply(sdfs, function(sdf) {
    q <- sesame::sesameQC_calcStats(
      sdf, funs = c("detection", "intensity", "numProbes", "channel", "dyeBias"))
    stat <- q@stat
    scalars <- stat[vapply(stat, function(x) length(x) == 1, logical(1))]
    as.data.frame(scalars, stringsAsFactors = FALSE)
  })
  all_cols <- unique(unlist(lapply(qc_list, names)))
  qc_list <- lapply(qc_list, function(d) {
    for (m in setdiff(all_cols, names(d))) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  qc_df <- do.call(rbind, qc_list)
  qc_df$sample_id <- names(sdfs)

  decomp <- lapply(sdfs, function(sdf) {
    tryCatch({
      n_total <- nrow(sdf)
      apriori_mask <- tryCatch({
        sc <- sdf; if ("mask" %in% colnames(sc)) sc$mask <- FALSE
        m <- sesame::qualityMask(sc)
        if (methods::is(m, "SigDF") || is.data.frame(m)) as.logical(m$mask)
        else if (is.logical(m)) m else rep(NA, n_total)
      }, error = function(e) rep(NA, n_total))
      final_mask <- if ("mask" %in% colnames(sdf)) as.logical(sdf$mask)
                    else rep(FALSE, n_total)
      if (all(is.na(apriori_mask)))
        return(data.frame(n_probes_total = n_total, n_masked_apriori = NA_integer_,
          n_failed_detection = NA_integer_, n_detected = NA_integer_,
          n_usable = sum(!final_mask, na.rm = TRUE),
          frac_apriori_mask = NA_real_,
          frac_detection_among_unmasked = NA_real_,
          frac_usable = sum(!final_mask, na.rm = TRUE) / n_total,
          stringsAsFactors = FALSE))
      apriori_mask[is.na(apriori_mask)] <- FALSE
      final_mask[is.na(final_mask)] <- FALSE
      nma <- sum(apriori_mask); nfd <- sum(final_mask & !apriori_mask)
      nu <- sum(!final_mask)
      data.frame(n_probes_total = n_total, n_masked_apriori = nma,
        n_failed_detection = nfd, n_detected = n_total - nma - nfd,
        n_usable = nu,
        frac_apriori_mask = nma / n_total,
        frac_detection_among_unmasked =
          if (n_total - nma > 0) (n_total - nma - nfd) / (n_total - nma)
          else NA_real_,
        frac_usable = nu / n_total, stringsAsFactors = FALSE)
    }, error = function(e) data.frame(n_probes_total = nrow(sdf),
      n_masked_apriori = NA_integer_, n_failed_detection = NA_integer_,
      n_detected = NA_integer_, n_usable = NA_integer_,
      frac_apriori_mask = NA_real_,
      frac_detection_among_unmasked = NA_real_, frac_usable = NA_real_,
      stringsAsFactors = FALSE))
  })
  decomp_df <- do.call(rbind, decomp)
  decomp_df$sample_id <- names(sdfs)
  qc_df <- merge(qc_df, decomp_df, by = "sample_id", sort = FALSE)

  qc_df$horvath_age <- vapply(sdfs, function(s)
    predict_horvath_age(sesame::getBetas(s)), numeric(1))

  sex_int <- compute_sex_intensities_from_sdfs(sdfs)
  if (!is.null(sex_int))
    qc_df <- merge(qc_df, sex_int, by = "sample_id", sort = FALSE)

  if (!is.null(logger))
    logger$log("qcmetrics", sprintf("frac_dt: median=%.3f min=%.3f",
               stats::median(qc_df$frac_dt, na.rm = TRUE),
               min(qc_df$frac_dt, na.rm = TRUE)))
  qc_df
}

#' Streaming per-sample QC metrics (reads IDATs one at a time)
#'
#' Identical output schema to \code{\link{qcmetrics}}.
#'
#' @param basenames Character vector of IDAT basenames.
#' @param platform Platform string.
#' @param logger Optional logger.
#' @return data.frame with one row per sample.
#' @export
qcstream <- function(basenames, platform, logger = NULL) {
  if (!is.null(logger))
    logger$log("qcstream",
               sprintf("streaming QC stats for %d samples...", length(basenames)))

  qc_rows <- lapply(basenames, function(bn) {
    sdf <- sesame::openSesame(bn, platform = platform, func = NULL)
    q <- sesame::sesameQC_calcStats(sdf,
      funs = c("detection", "intensity", "numProbes", "channel", "dyeBias"))
    stat <- q@stat
    scalars <- stat[vapply(stat, function(x) length(x) == 1, logical(1))]
    qr <- as.data.frame(scalars, stringsAsFactors = FALSE)

    n_total <- nrow(sdf)
    sc <- sdf; if ("mask" %in% colnames(sc)) sc$mask <- FALSE
    ap <- tryCatch({
      m <- sesame::qualityMask(sc)
      if (methods::is(m, "SigDF") || is.data.frame(m)) as.logical(m$mask)
      else if (is.logical(m)) m else rep(NA, n_total)
    }, error = function(e) rep(NA, n_total))
    fm <- if ("mask" %in% colnames(sdf)) as.logical(sdf$mask)
          else rep(FALSE, n_total)

    if (all(is.na(ap))) {
      qr$n_probes_total <- n_total; qr$n_masked_apriori <- NA_integer_
      qr$n_failed_detection <- NA_integer_; qr$n_detected <- NA_integer_
      qr$n_usable <- sum(!fm, na.rm = TRUE)
      qr$frac_apriori_mask <- NA_real_
      qr$frac_detection_among_unmasked <- NA_real_
      qr$frac_usable <- qr$n_usable / n_total
    } else {
      ap[is.na(ap)] <- FALSE; fm[is.na(fm)] <- FALSE
      nma <- sum(ap); nfd <- sum(fm & !ap); nu <- sum(!fm)
      qr$n_probes_total <- n_total; qr$n_masked_apriori <- nma
      qr$n_failed_detection <- nfd; qr$n_detected <- n_total - nma - nfd
      qr$n_usable <- nu; qr$frac_apriori_mask <- nma / n_total
      qr$frac_detection_among_unmasked <-
        if (n_total - nma > 0) qr$n_detected / (n_total - nma) else NA_real_
      qr$frac_usable <- nu / n_total
    }
    qr$horvath_age <- predict_horvath_age(sesame::getBetas(sdf))
    si <- compute_sex_intensities_single(sdf)
    qr$sex_chrX_intensity <- si$chrX
    qr$sex_chrY_intensity <- si$chrY
    rm(sdf); gc(verbose = FALSE)
    qr
  })

  all_cols <- unique(unlist(lapply(qc_rows, names)))
  qc_rows <- lapply(qc_rows, function(d) {
    for (m in setdiff(all_cols, names(d))) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  qc_df <- do.call(rbind, qc_rows)
  qc_df$sample_id <- basename(basenames)
  if (!is.null(logger))
    logger$log("qcstream", sprintf("frac_dt: median=%.3f min=%.3f",
               stats::median(qc_df$frac_dt, na.rm = TRUE),
               min(qc_df$frac_dt, na.rm = TRUE)))
  qc_df
}

#' Predict Horvath (2013) epigenetic age
#' @keywords internal
#' @noRd
predict_horvath_age <- function(betas) {
  tryCatch({
    if (is.null(names(betas))) return(NA_real_)
    present <- intersect(.horvath_probes, names(betas))
    if (length(present) < 200) return(NA_real_)
    idx <- match(present, .horvath_probes)
    bv <- betas[present]; bv[is.na(bv)] <- 0.5
    x <- .horvath_intercept + sum(.horvath_coefficients[idx] * bv)
    if (x < 0) 21 * exp(x) - 1 else 21 * x + 20
  }, error = function(e) NA_real_)
}

#' Sex intensities from single SigDF
#' @keywords internal
#' @noRd
compute_sex_intensities_single <- function(sdf) {
  pid <- sdf$Probe_ID
  ti <- rowSums(cbind(ifelse(is.na(sdf$MG), 0, sdf$MG),
                       ifelse(is.na(sdf$MR), 0, sdf$MR),
                       ifelse(is.na(sdf$UG), 0, sdf$UG),
                       ifelse(is.na(sdf$UR), 0, sdf$UR)),
                na.rm = TRUE)
  xi <- which(pid %in% .sesame_chrX_xlinked)
  yi <- which(pid %in% .sesame_chrY_clean)
  list(chrX = if (length(xi) >= 10) stats::median(ti[xi], na.rm = TRUE) else NA_real_,
       chrY = if (length(yi) >= 10) stats::median(ti[yi], na.rm = TRUE) else NA_real_)
}

#' Sex intensities from list of SigDFs
#' @keywords internal
#' @noRd
compute_sex_intensities_from_sdfs <- function(sdfs) {
  rows <- lapply(seq_along(sdfs), function(i) {
    si <- compute_sex_intensities_single(sdfs[[i]])
    data.frame(sample_id = names(sdfs)[i],
               sex_chrX_intensity = si$chrX,
               sex_chrY_intensity = si$chrY,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Cross-check reported metadata against computed QC metrics (logging only)
#'
#' @param ss Sample sheet.
#' @param qcdf QC metrics data.frame from \code{\link{qcmetrics}} or
#'   \code{\link{qcstream}}.
#' @param logger Optional logger.
#' @return Invisibly NULL.
#' @export
checkmeta <- function(ss, qcdf, logger = NULL) {
  cfg <- mqcopts()
  m <- match(qcdf$sample_id, ss$sample_id)
  age_col <- resolve_column(colnames(ss), cfg$agecol, cfg$agealiases)
  if (!is.na(age_col) && !is.null(logger)) {
    age <- parse_age_robust(ss[[age_col]][m])
    nv <- sum(!is.na(age) & !is.na(qcdf$horvath_age))
    logger$log("checkmeta",
               sprintf("horvath age: %d / %d valid predictions",
                       nv, nrow(qcdf)))
    if (nv >= 3) {
      r <- suppressWarnings(stats::cor(qcdf$horvath_age, age,
                                       use = "complete.obs"))
      logger$log("checkmeta",
                 sprintf("horvath vs reported: r=%.3f", r))
    }
  }
  invisible(NULL)
}

#' Robust age parsing
#' @keywords internal
#' @noRd
parse_age_robust <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x[tolower(trimws(x)) %in% c("na", "n/a", "")] <- NA
  cleaned <- gsub("[^0-9.\\-]", "", x)
  cleaned[cleaned == "" | cleaned == "." | cleaned == "-"] <- NA
  suppressWarnings(as.numeric(cleaned))
}
