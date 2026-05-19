#' Compute per-sample QC metrics from a list of SigDFs
#' @export
compute_qc_metrics <- function(sdfs, platform, logger = NULL) {
  if (!is.null(logger))
    logger$log("qc_metrics", "computing per-sample QC stats...")

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
        return(data.frame(n_probes_total=n_total, n_masked_apriori=NA_integer_,
          n_failed_detection=NA_integer_, n_detected=NA_integer_,
          n_usable=sum(!final_mask,na.rm=TRUE), frac_apriori_mask=NA_real_,
          frac_detection_among_unmasked=NA_real_,
          frac_usable=sum(!final_mask,na.rm=TRUE)/n_total, stringsAsFactors=FALSE))
      apriori_mask[is.na(apriori_mask)] <- FALSE
      final_mask[is.na(final_mask)] <- FALSE
      nma <- sum(apriori_mask); nfd <- sum(final_mask & !apriori_mask)
      nu <- sum(!final_mask)
      data.frame(n_probes_total=n_total, n_masked_apriori=nma,
        n_failed_detection=nfd, n_detected=n_total-nma-nfd, n_usable=nu,
        frac_apriori_mask=nma/n_total,
        frac_detection_among_unmasked=if(n_total-nma>0)(n_total-nma-nfd)/(n_total-nma) else NA_real_,
        frac_usable=nu/n_total, stringsAsFactors=FALSE)
    }, error = function(e) data.frame(n_probes_total=nrow(sdf),
      n_masked_apriori=NA_integer_, n_failed_detection=NA_integer_,
      n_detected=NA_integer_, n_usable=NA_integer_, frac_apriori_mask=NA_real_,
      frac_detection_among_unmasked=NA_real_, frac_usable=NA_real_,
      stringsAsFactors=FALSE))
  })
  decomp_df <- do.call(rbind, decomp)
  decomp_df$sample_id <- names(sdfs)
  qc_df <- merge(qc_df, decomp_df, by = "sample_id", sort = FALSE)

  # Horvath age: unmasked, noob-corrected betas (sdfs are final QCDPB
  # SigDFs when keep_sdfs = TRUE). mask = FALSE keeps every probe.
  qc_df$horvath_age <- vapply(sdfs, function(s)
    predict_horvath_age(sesame::getBetas(s, mask = FALSE)), numeric(1))

  sex_int <- compute_sex_intensities_from_sdfs(sdfs)
  if (!is.null(sex_int)) qc_df <- merge(qc_df, sex_int, by="sample_id", sort=FALSE)

  if (!is.null(logger))
    logger$log("qc_metrics", sprintf("frac_dt: median=%.3f min=%.3f",
               stats::median(qc_df$frac_dt,na.rm=TRUE), min(qc_df$frac_dt,na.rm=TRUE)))
  qc_df
}

#' Streaming QC metrics
#'
#' For each sample two SigDFs are derived: a pre-noob SigDF (prep code
#' up to and including P) used for the QC statistics, sex chromosome
#' intensities, and mean intensity; and a noob-corrected SigDF (full
#' prep code) used only for the Horvath clock betas. Sex and intensity
#' metrics intentionally use the non-noob signal (noob only subtracts
#' background and cannot recover weak signal, so it offers no benefit
#' for chromosome-presence calls).
#' @export
compute_qc_metrics_streaming <- function(basenames, platform, logger = NULL) {
  cfg <- methylQC_options()
  prep_code <- cfg$prep_code
  prep_preP <- sub("B.*$", "", prep_code)
  has_noob  <- grepl("B", prep_code)
  if (!is.null(logger))
    logger$log("qc_metrics", sprintf("streaming QC stats for %d samples...", length(basenames)))

  qc_rows <- lapply(basenames, function(bn) {
    # Pre-noob SigDF: QC stats, sex intensities, mean intensity
    sdf <- sesame::openSesame(bn, platform=platform, prep=prep_preP, func=NULL)
    q <- sesame::sesameQC_calcStats(sdf, funs=c("detection","intensity","numProbes","channel","dyeBias"))
    stat <- q@stat
    scalars <- stat[vapply(stat, function(x) length(x)==1, logical(1))]
    qr <- as.data.frame(scalars, stringsAsFactors=FALSE)

    n_total <- nrow(sdf)
    sc <- sdf; if ("mask" %in% colnames(sc)) sc$mask <- FALSE
    ap <- tryCatch({
      m <- sesame::qualityMask(sc)
      if (methods::is(m,"SigDF")||is.data.frame(m)) as.logical(m$mask)
      else if (is.logical(m)) m else rep(NA,n_total)
    }, error=function(e) rep(NA,n_total))
    fm <- if ("mask" %in% colnames(sdf)) as.logical(sdf$mask) else rep(FALSE,n_total)

    if (all(is.na(ap))) {
      qr$n_probes_total <- n_total; qr$n_masked_apriori <- NA_integer_
      qr$n_failed_detection <- NA_integer_; qr$n_detected <- NA_integer_
      qr$n_usable <- sum(!fm,na.rm=TRUE); qr$frac_apriori_mask <- NA_real_
      qr$frac_detection_among_unmasked <- NA_real_
      qr$frac_usable <- qr$n_usable/n_total
    } else {
      ap[is.na(ap)] <- FALSE; fm[is.na(fm)] <- FALSE
      nma <- sum(ap); nfd <- sum(fm & !ap); nu <- sum(!fm)
      qr$n_probes_total <- n_total; qr$n_masked_apriori <- nma
      qr$n_failed_detection <- nfd; qr$n_detected <- n_total-nma-nfd
      qr$n_usable <- nu; qr$frac_apriori_mask <- nma/n_total
      qr$frac_detection_among_unmasked <- if(n_total-nma>0) qr$n_detected/(n_total-nma) else NA_real_
      qr$frac_usable <- nu/n_total
    }
    # Sex chromosome intensities: non-noob signal
    si <- compute_sex_intensities_single(sdf)
    qr$sex_chrX_intensity <- si$chrX; qr$sex_chrY_intensity <- si$chrY

    # Horvath clock: noob-corrected, unmasked betas
    sdf_noob <- if (has_noob) {
      sesame::openSesame(bn, platform=platform, prep=prep_code, func=NULL)
    } else sdf
    qr$horvath_age <- predict_horvath_age(
      sesame::getBetas(sdf_noob, mask = FALSE))

    rm(sdf, sdf_noob); gc(verbose=FALSE)
    qr
  })

  all_cols <- unique(unlist(lapply(qc_rows, names)))
  qc_rows <- lapply(qc_rows, function(d) {
    for (m in setdiff(all_cols, names(d))) d[[m]] <- NA
    d[, all_cols, drop=FALSE]
  })
  qc_df <- do.call(rbind, qc_rows)
  qc_df$sample_id <- basename(basenames)
  if (!is.null(logger))
    logger$log("qc_metrics", sprintf("frac_dt: median=%.3f min=%.3f",
               stats::median(qc_df$frac_dt,na.rm=TRUE), min(qc_df$frac_dt,na.rm=TRUE)))
  qc_df
}

#' Predict Horvath (2013) epigenetic age
#'
#' Per-CpG value precedence for the 353 clock probes:
#' \enumerate{
#'   \item Probe present on the platform: use the actual (noob-corrected,
#'     unmasked) beta value. A present-but-failed probe still contributes
#'     its real value; the clock is a coarse QC screen and is not
#'     re-imputed here.
#'   \item Probe absent from the platform manifest: use the hard-coded
#'     zero-shot reference median (see \code{.horvath_zeroshot_betas}),
#'     derived from HM450 blood samples.
#' }
#' There is no 0.5 fallback. A clock CpG that is neither present nor in
#' the zero-shot reference contributes nothing (its coefficient is
#' dropped), and the prediction is returned as NA if fewer than 200 of
#' the 353 CpGs can be resolved.
#'
#' @param betas Named numeric vector of beta values.
#' @return Predicted epigenetic age in years, or NA.
#' @keywords internal
#' @noRd
predict_horvath_age <- function(betas) {
  tryCatch({
    if (is.null(names(betas))) return(NA_real_)
    bv <- rep(NA_real_, length(.horvath_probes))
    names(bv) <- .horvath_probes

    # (1) Probe present on the platform -> actual beta
    present <- intersect(.horvath_probes, names(betas))
    if (length(present)) bv[present] <- betas[present]

    # (2) Probe absent from the platform -> zero-shot reference median
    absent <- .horvath_probes[is.na(bv)]
    zs <- intersect(absent, names(.horvath_zeroshot_betas))
    if (length(zs)) bv[zs] <- .horvath_zeroshot_betas[zs]

    # A present-but-failed probe may still be NA in `betas`; keep its
    # actual (possibly NA) value rather than imputing. Drop only the
    # CpGs that remain unresolvable.
    usable <- which(!is.na(bv))
    if (length(usable) < 200) return(NA_real_)

    x <- .horvath_intercept +
         sum(.horvath_coefficients[usable] * bv[usable])
    if (x < 0) 21 * exp(x) - 1 else 21 * x + 20
  }, error = function(e) NA_real_)
}

#' Sex intensities from single SigDF
#' @keywords internal
#' @noRd
compute_sex_intensities_single <- function(sdf) {
  pid <- sdf$Probe_ID
  ti <- rowSums(cbind(ifelse(is.na(sdf$MG),0,sdf$MG), ifelse(is.na(sdf$MR),0,sdf$MR),
                       ifelse(is.na(sdf$UG),0,sdf$UG), ifelse(is.na(sdf$UR),0,sdf$UR)),
                na.rm=TRUE)
  xi <- which(pid %in% .sesame_chrX_xlinked); yi <- which(pid %in% .sesame_chrY_clean)
  list(chrX = if(length(xi)>=10) stats::median(ti[xi],na.rm=TRUE) else NA_real_,
       chrY = if(length(yi)>=10) stats::median(ti[yi],na.rm=TRUE) else NA_real_)
}

#' Sex intensities from list of SigDFs
#' @keywords internal
#' @noRd
compute_sex_intensities_from_sdfs <- function(sdfs) {
  rows <- lapply(seq_along(sdfs), function(i) {
    si <- compute_sex_intensities_single(sdfs[[i]])
    data.frame(sample_id=names(sdfs)[i], sex_chrX_intensity=si$chrX,
               sex_chrY_intensity=si$chrY, stringsAsFactors=FALSE)
  })
  do.call(rbind, rows)
}

#' Cross-check reported metadata against computed QC metrics (logging only)
#' @export
check_sample_metadata <- function(sample_sheet, qc_df, logger = NULL) {
  cfg <- methylQC_options()
  m <- match(qc_df$sample_id, sample_sheet$sample_id)
  age_col <- resolve_column(colnames(sample_sheet), cfg$reported_age_col, cfg$reported_age_aliases)
  if (!is.na(age_col) && !is.null(logger)) {
    age <- parse_age_robust(sample_sheet[[age_col]][m])
    nv <- sum(!is.na(age) & !is.na(qc_df$horvath_age))
    logger$log("metadata_check", sprintf("horvath age: %d / %d valid predictions", nv, nrow(qc_df)))
    if (nv >= 3) {
      r <- suppressWarnings(stats::cor(qc_df$horvath_age, age, use="complete.obs"))
      logger$log("metadata_check", sprintf("horvath vs reported: r=%.3f", r))
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
