###############################################################################
# snp_identity.R — Sample identity verification via SNP probe MDS
#
# Uses the 59 rs (SNP genotyping) probes to verify sample identity.
# These probes target common SNPs and produce genotype-clustered beta
# values (~0 for AA, ~0.5 for AB, ~1 for BB), creating a genetic
# fingerprint per individual.
#
# Algorithm:
#   1. Impute NAs with probe medians (rs probes have clear genotype clusters)
#   2. Compute Euclidean distance matrix across all samples
#   3. Classical MDS (cmdscale) to 2 dimensions
#   4. If donor IDs are available:
#      a. Compute per-donor centroids in MDS space
#      b. For each sample, find the nearest centroid
#      c. Flag samples whose nearest centroid belongs to a different donor
#         than their reported ID (potential sample swap)
#   5. Output: 2-page PDF (all samples + zoomed flagged donors), MDS
#      coordinates CSV, and flags CSV
###############################################################################
#' Sample identity check via MDS on SNP (rs) probes
#'
#' Runs multidimensional scaling on the rs probe beta matrix, merges with
#' sample sheet metadata to identify donor IDs, plots each sample as its
#' donor label, and flags samples that cluster near a different donor's
#' group.
#'
#' @param rsbetas Numeric matrix with samples as rows, rs probes as
#'   columns (as returned by \code{\link{snpbetas}}). May also be the
#'   path to an RDS file containing such a matrix.
#' @param ss Optional data.frame or path to a CSV sample sheet. If
#'   provided, merged with \code{rsbetas} row names to look up donor IDs.
#'   If NULL, only the MDS coordinates CSV is produced.
#' @param donorcol Optional column name in the sample sheet containing
#'   donor/participant IDs. If NULL, the function tries to auto-detect
#'   from aliases.
#' @param outdir Output directory (created if needed).
#' @param width Width of output PDFs in inches (default 12).
#' @param height Height of output PDFs in inches (default 10).
#' @param logger Optional logger.
#' @return Invisibly, a list with \code{mds_coords} (data.frame),
#'   \code{flagged} (data.frame of mismatched samples), and
#'   \code{donorcol} (resolved column name).
#' @export
snpcheck <- function(rsbetas,
                     ss = NULL,
                     donorcol = NULL,
                     outdir = ".",
                     width = 12,
                     height = 10,
                     logger = NULL) {
  # Load rs betas from file if a path was given
  if (is.character(rsbetas) && length(rsbetas) == 1) {
    if (!file.exists(rsbetas)) stop("rsbetas file not found: ", rsbetas)
    rsbetas <- readRDS(rsbetas)
  }
  stopifnot(is.matrix(rsbetas) && is.numeric(rsbetas))

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # --- Step 1: Impute NAs with probe medians ---
  # rs probes cluster at 0/0.5/1 (genotypes), so median is sensible
  na_count <- sum(is.na(rsbetas))
  if (na_count > 0) {
    col_medians <- apply(rsbetas, 2, stats::median, na.rm = TRUE)
    for (j in seq_len(ncol(rsbetas))) {
      na_rows <- is.na(rsbetas[, j])
      rsbetas[na_rows, j] <- col_medians[j]
    }
    if (!is.null(logger)) {
      logger$log("snpcheck",
                 sprintf("imputed %d NAs with probe medians", na_count))
    }
  }

  # --- Step 2: MDS on Euclidean distances ---
  d <- stats::dist(rsbetas, method = "euclidean")
  mds <- stats::cmdscale(d, k = 2)
  colnames(mds) <- c("MDS1", "MDS2")
  rm(d); gc(verbose = FALSE)

  mds_df <- data.frame(
    sample_id = rownames(mds),
    MDS1 = mds[, 1], MDS2 = mds[, 2],
    stringsAsFactors = FALSE)
  rm(mds); gc(verbose = FALSE)

  # Load and merge sample sheet
  if (!is.null(ss)) {
    if (is.character(ss) && length(ss) == 1)
      ss <- utils::read.csv(ss, stringsAsFactors = FALSE)
    stopifnot(is.data.frame(ss))

    id_col <- find_matching_id_column(ss, mds_df$sample_id)
    if (!is.null(id_col)) {
      mds_df <- merge(mds_df, ss,
                      by.x = "sample_id", by.y = id_col,
                      all.x = TRUE, sort = FALSE)
      if (!is.null(logger)) {
        logger$log("snpcheck",
                   sprintf("merged sample sheet via column '%s'", id_col))
      }
    } else if (!is.null(logger)) {
      logger$log("snpcheck",
                 "WARNING: could not find matching ID column in sample sheet")
    }
  }

  # Resolve donor column
  cfg <- mqcopts()
  if (!is.null(donorcol)) {
    if (!donorcol %in% colnames(mds_df))
      stop("donorcol '", donorcol,
           "' not found. Available: ",
           paste(colnames(mds_df), collapse = ", "))
    has_donor <- TRUE
  } else {
    aliases <- cfg$donoraliases
    donorcol <- resolve_column(colnames(mds_df), cfg$donorcol, aliases)
    has_donor <- !is.na(donorcol)
  }

  # Write MDS coordinates CSV
  utils::write.csv(mds_df, file.path(outdir, "snp_mds_coordinates.csv"),
                   row.names = FALSE)
  if (!is.null(logger)) {
    logger$log("snpcheck",
               sprintf("wrote snp_mds_coordinates.csv (%d samples)",
                       nrow(mds_df)))
  }

  # No donor column: simple plot
  if (!has_donor) {
    grDevices::pdf(file.path(outdir, "snp_identity_mds.pdf"),
                   width = width, height = height)
    plot(mds_df$MDS1, mds_df$MDS2, type = "n",
         xlab = "MDS1", ylab = "MDS2",
         main = "SNP probe MDS (no donor ID available)")
    graphics::text(mds_df$MDS1, mds_df$MDS2,
                   labels = mds_df$sample_id, cex = 0.6)
    grDevices::dev.off()
    if (!is.null(logger)) {
      logger$log("snpcheck",
                 "no donor column -- plotted sample IDs only")
    }
    return(invisible(list(mds_coords = mds_df, flagged = data.frame(),
                          donorcol = NA_character_)))
  }

  # --- Step 3: Per-donor centroids ---
  donor <- mds_df[[donorcol]]
  donors <- unique(donor[!is.na(donor)])

  centroids <- do.call(rbind, lapply(donors, function(pid) {
    rows <- which(donor == pid)
    data.frame(donor  = pid,
               c_MDS1 = mean(mds_df$MDS1[rows]),
               c_MDS2 = mean(mds_df$MDS2[rows]),
               stringsAsFactors = FALSE)
  }))

  # --- Step 4: For each sample, find nearest centroid ---
  mds_df$nearest_donor    <- NA_character_
  mds_df$dist_to_own      <- NA_real_
  mds_df$dist_to_nearest  <- NA_real_

  for (i in seq_len(nrow(mds_df))) {
    dists <- sqrt((centroids$c_MDS1 - mds_df$MDS1[i])^2 +
                  (centroids$c_MDS2 - mds_df$MDS2[i])^2)
    nearest_idx <- which.min(dists)
    mds_df$nearest_donor[i]   <- centroids$donor[nearest_idx]
    mds_df$dist_to_nearest[i] <- dists[nearest_idx]

    own_idx <- which(centroids$donor == donor[i])
    if (length(own_idx)) mds_df$dist_to_own[i] <- dists[own_idx]
  }

  mds_df$flagged <- !is.na(donor) &
                    !is.na(mds_df$nearest_donor) &
                    donor != mds_df$nearest_donor

  flagged_df <- mds_df[mds_df$flagged,
                       c("sample_id", donorcol,
                         "nearest_donor",
                         "dist_to_own", "dist_to_nearest"),
                       drop = FALSE]
  colnames(flagged_df)[2] <- "reported_donor"

  if (nrow(flagged_df) > 0) {
    utils::write.csv(flagged_df,
                     file.path(outdir, "snp_identity_flags.csv"),
                     row.names = FALSE)
    if (!is.null(logger)) {
      logger$log("snpcheck",
                 sprintf("WARNING: %d sample(s) cluster near wrong donor",
                         nrow(flagged_df)))
    }
  } else if (!is.null(logger)) {
    logger$log("snpcheck",
               "all samples cluster with their reported donor")
  }

  # --- Step 5: Plots ---
  grDevices::pdf(file.path(outdir, "snp_identity_mds.pdf"),
                 width = width, height = height)

  plot(mds_df$MDS1, mds_df$MDS2, type = "n",
       xlab = "MDS1", ylab = "MDS2",
       main = sprintf("SNP probe MDS - labelled by %s", donorcol),
       sub  = sprintf("%d samples, %d donors, %d flagged",
                      nrow(mds_df), length(donors), nrow(flagged_df)),
       cex.lab = 1.3, cex.main = 1.4)

  ok <- !mds_df$flagged & !is.na(donor)
  if (any(ok))
    graphics::text(mds_df$MDS1[ok], mds_df$MDS2[ok],
                   labels = donor[ok], col = "black", cex = 0.55)

  if (any(mds_df$flagged))
    graphics::text(mds_df$MDS1[mds_df$flagged],
                   mds_df$MDS2[mds_df$flagged],
                   labels = donor[mds_df$flagged],
                   col = "red", cex = 0.65, font = 2)

  no_pid <- is.na(donor)
  if (any(no_pid))
    graphics::text(mds_df$MDS1[no_pid], mds_df$MDS2[no_pid],
                   labels = mds_df$sample_id[no_pid],
                   col = "grey50", cex = 0.45)

  # Page 2: zoomed-in plot of involved donors
  if (nrow(flagged_df) > 0) {
    involved_pids <- unique(c(flagged_df$reported_donor,
                              flagged_df$nearest_donor))
    involved_idx  <- which(donor %in% involved_pids)

    if (length(involved_idx) >= 2) {
      sub_df    <- mds_df[involved_idx, ]
      sub_donor <- donor[involved_idx]

      x_range <- range(sub_df$MDS1)
      y_range <- range(sub_df$MDS2)
      x_pad <- diff(x_range) * 0.15
      y_pad <- diff(y_range) * 0.15

      plot(sub_df$MDS1, sub_df$MDS2, type = "n",
           xlim = x_range + c(-x_pad, x_pad),
           ylim = y_range + c(-y_pad, y_pad),
           xlab = "MDS1", ylab = "MDS2",
           main = sprintf("Flagged donors only (%d donors, %d samples)",
                          length(involved_pids), nrow(sub_df)),
           sub  = "Showing all samples from donors with flagged samples",
           cex.lab = 1.3, cex.main = 1.4)

      sub_ok <- !sub_df$flagged
      if (any(sub_ok))
        graphics::text(sub_df$MDS1[sub_ok], sub_df$MDS2[sub_ok],
                       labels = sub_donor[sub_ok],
                       col = "black", cex = 0.7)

      sub_flagged <- sub_df$flagged
      if (any(sub_flagged))
        graphics::text(sub_df$MDS1[sub_flagged], sub_df$MDS2[sub_flagged],
                       labels = sub_donor[sub_flagged],
                       col = "red", cex = 0.8, font = 2)

      involved_centroids <- centroids[centroids$donor %in% involved_pids, ]
      for (r in seq_len(nrow(involved_centroids))) {
        graphics::points(involved_centroids$c_MDS1[r],
                         involved_centroids$c_MDS2[r],
                         pch = 3, cex = 2, col = "grey40", lwd = 2)
        graphics::text(involved_centroids$c_MDS1[r],
                       involved_centroids$c_MDS2[r],
                       labels = involved_centroids$donor[r],
                       pos = 1, cex = 0.6, col = "grey40", font = 3)
      }
    }
  }

  grDevices::dev.off()

  if (!is.null(logger)) {
    logger$log("snpcheck",
               sprintf("wrote snp_identity_mds.pdf (%d pages)",
                       if (nrow(flagged_df) > 0) 2 else 1))
  }

  invisible(list(mds_coords = mds_df, flagged = flagged_df,
                 donorcol = donorcol))
}

#' Smart-match a sample ID column in a data.frame to a vector of IDs
#' @keywords internal
#' @noRd
find_matching_id_column <- function(df, ids) {
  best_col  <- NULL
  best_frac <- 0

  for (col in colnames(df)) {
    vals <- as.character(df[[col]])
    frac <- mean(ids %in% vals)
    if (frac > best_frac) {
      best_frac <- frac
      best_col  <- col
    }
  }
  if (best_frac >= 0.5) best_col else NULL
}
