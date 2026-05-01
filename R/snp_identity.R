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
#   4. If participant/donor IDs are available:
#      a. Compute per-participant centroids in MDS space
#      b. For each sample, find the nearest centroid
#      c. Flag samples whose nearest centroid belongs to a different
#         participant than their reported ID (potential sample swap)
#   5. Output: 2-page PDF (all samples + zoomed flagged participants),
#      MDS coordinates CSV, and flags CSV
#
# Key functions:
#   check_snp_identity()     — main entry point (exported)
#   find_matching_id_column() — auto-detect which column matches sample IDs
###############################################################################
#' Sample identity check via MDS on SNP (rs) probes
#'
#' Runs multidimensional scaling on the rs probe beta matrix, merges
#' with sample sheet metadata to identify participant IDs, plots each
#' sample as its participant label, and flags samples that cluster
#' near a different participant's group.
#'
#' @param rs_betas Numeric matrix with samples as rows, rs probes as
#'   columns (as returned by \code{\link{extract_snp_betas}}).
#'   Alternatively, an RDS file path containing such a matrix.
#' @param sample_sheet Optional data.frame or file path to a CSV sample
#'   sheet. If provided, merged with rs_betas row names to look up
#'   participant IDs. If NULL, only the MDS coordinates CSV is produced.
#' @param participant_col Optional column name in the sample sheet
#'   containing participant/donor IDs. If NULL, the function tries to
#'   auto-detect from aliases.
#' @param out_dir Output directory (created if needed).
#' @param pdf_width Width of output PDFs in inches (default 12).
#' @param pdf_height Height of output PDFs in inches (default 10).
#' @param logger Optional logger.
#' @return Invisibly, a list with mds_coords (data.frame), flagged
#'   (data.frame of mismatched samples), and participant_col (resolved).
#' @export
check_snp_identity <- function(rs_betas,
                                sample_sheet = NULL,
                                participant_col = NULL,
                                out_dir = ".",
                                pdf_width = 12,
                                pdf_height = 10,
                                logger = NULL) {
  # Load rs_betas from file if path given
  if (is.character(rs_betas) && length(rs_betas) == 1) {
    if (!file.exists(rs_betas)) stop("rs_betas file not found: ", rs_betas)
    rs_betas <- readRDS(rs_betas)
  }
  stopifnot(is.matrix(rs_betas) && is.numeric(rs_betas))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Impute NAs with probe median (rs probes cluster at 0/0.5/1)
  # --- Step 1: Impute NAs with probe medians ---
  # rs probes cluster at 0/0.5/1 (genotypes), so median is sensible
  na_count <- sum(is.na(rs_betas))
  if (na_count > 0) {
    col_medians <- apply(rs_betas, 2, stats::median, na.rm = TRUE)
    for (j in seq_len(ncol(rs_betas))) {
      na_rows <- is.na(rs_betas[, j])
      rs_betas[na_rows, j] <- col_medians[j]
    }
    if (!is.null(logger))
      logger$log("snp_identity",
                 sprintf("imputed %d NAs with probe medians", na_count))
  }

  # MDS
  # --- Step 2: MDS on Euclidean distances ---
  d <- stats::dist(rs_betas, method = "euclidean")
  mds <- stats::cmdscale(d, k = 2)
  colnames(mds) <- c("MDS1", "MDS2")
  rm(d); gc(verbose = FALSE)

  mds_df <- data.frame(
    sample_id = rownames(mds),
    MDS1 = mds[, 1], MDS2 = mds[, 2],
    stringsAsFactors = FALSE)
  rm(mds); gc(verbose = FALSE)

  # Load and merge sample sheet
  if (!is.null(sample_sheet)) {
    if (is.character(sample_sheet) && length(sample_sheet) == 1)
      sample_sheet <- utils::read.csv(sample_sheet, stringsAsFactors = FALSE)
    stopifnot(is.data.frame(sample_sheet))

    id_col <- find_matching_id_column(sample_sheet, mds_df$sample_id)
    if (!is.null(id_col)) {
      mds_df <- merge(mds_df, sample_sheet,
                      by.x = "sample_id", by.y = id_col,
                      all.x = TRUE, sort = FALSE)
      if (!is.null(logger))
        logger$log("snp_identity",
                   sprintf("merged sample sheet via column '%s'", id_col))
    } else if (!is.null(logger)) {
      logger$log("snp_identity",
                 "WARNING: could not find matching ID column in sample sheet")
    }
  }

  # Resolve participant column
  if (!is.null(participant_col)) {
    if (!participant_col %in% colnames(mds_df))
      stop("participant_col '", participant_col,
           "' not found. Available: ",
           paste(colnames(mds_df), collapse = ", "))
    has_participant <- TRUE
  } else {
    participant_aliases <- c("Donor", "donor", "Subject", "subject",
                             "Participant", "participant", "SubjectID",
                             "DonorID", "ID", "donor_id", "subject_id")
    participant_col <- resolve_column(colnames(mds_df), "Donor",
                                      participant_aliases)
    has_participant <- !is.na(participant_col)
  }

  # Write MDS coordinates CSV
  utils::write.csv(mds_df, file.path(out_dir, "snp_mds_coordinates.csv"),
                   row.names = FALSE)
  if (!is.null(logger))
    logger$log("snp_identity",
               sprintf("wrote snp_mds_coordinates.csv (%d samples)", nrow(mds_df)))

  # No participant column: simple plot
  if (!has_participant) {
    grDevices::pdf(file.path(out_dir, "snp_identity_mds.pdf"),
                   width = pdf_width, height = pdf_height)
    plot(mds_df$MDS1, mds_df$MDS2, type = "n",
         xlab = "MDS1", ylab = "MDS2",
         main = "SNP probe MDS (no participant ID available)")
    graphics::text(mds_df$MDS1, mds_df$MDS2,
                   labels = mds_df$sample_id, cex = 0.6)
    grDevices::dev.off()
    if (!is.null(logger))
      logger$log("snp_identity",
                 "no participant column -- plotted sample IDs only")
    return(invisible(list(mds_coords = mds_df, flagged = data.frame(),
                          participant_col = NA_character_)))
  }

  # --- Identity checking ---
  participant <- mds_df[[participant_col]]
  participants <- unique(participant[!is.na(participant)])

  # Per-participant centroids
  # --- Step 3: Compute per-participant centroids ---
  centroids <- do.call(rbind, lapply(participants, function(pid) {
    rows <- which(participant == pid)
    data.frame(participant = pid,
               c_MDS1 = mean(mds_df$MDS1[rows]),
               c_MDS2 = mean(mds_df$MDS2[rows]),
               stringsAsFactors = FALSE)
  }))

  # For each sample, find nearest centroid
  mds_df$nearest_participant <- NA_character_
  mds_df$dist_to_own <- NA_real_
  mds_df$dist_to_nearest <- NA_real_

  # --- Step 4: For each sample, find nearest centroid ---
  # Flag samples whose nearest centroid belongs to a different participant
  for (i in seq_len(nrow(mds_df))) {
    dists <- sqrt((centroids$c_MDS1 - mds_df$MDS1[i])^2 +
                  (centroids$c_MDS2 - mds_df$MDS2[i])^2)
    nearest_idx <- which.min(dists)
    mds_df$nearest_participant[i] <- centroids$participant[nearest_idx]
    mds_df$dist_to_nearest[i] <- dists[nearest_idx]

    own_idx <- which(centroids$participant == participant[i])
    if (length(own_idx))
      mds_df$dist_to_own[i] <- dists[own_idx]
  }

  mds_df$flagged <- !is.na(participant) &
                    !is.na(mds_df$nearest_participant) &
                    participant != mds_df$nearest_participant

  flagged_df <- mds_df[mds_df$flagged,
                       c("sample_id", participant_col,
                         "nearest_participant",
                         "dist_to_own", "dist_to_nearest"),
                       drop = FALSE]
  colnames(flagged_df)[2] <- "reported_participant"

  if (nrow(flagged_df) > 0) {
    utils::write.csv(flagged_df,
                     file.path(out_dir, "snp_identity_flags.csv"),
                     row.names = FALSE)
    if (!is.null(logger))
      logger$log("snp_identity",
                 sprintf("WARNING: %d sample(s) cluster near wrong participant",
                         nrow(flagged_df)))
  } else if (!is.null(logger)) {
    logger$log("snp_identity",
               "all samples cluster with their reported participant")
  }

  # --- Step 5: Generate plots ---
  # --- Plot 1: All samples ---
  grDevices::pdf(file.path(out_dir, "snp_identity_mds.pdf"),
                 width = pdf_width, height = pdf_height)

  plot(mds_df$MDS1, mds_df$MDS2, type = "n",
       xlab = "MDS1", ylab = "MDS2",
       main = sprintf("SNP probe MDS - labeled by %s", participant_col),
       sub  = sprintf("%d samples, %d participants, %d flagged",
                      nrow(mds_df), length(participants), nrow(flagged_df)),
       cex.lab = 1.3, cex.main = 1.4)

  # Non-flagged: black text
  ok <- !mds_df$flagged & !is.na(participant)
  if (any(ok))
    graphics::text(mds_df$MDS1[ok], mds_df$MDS2[ok],
                   labels = participant[ok], col = "black", cex = 0.55)

  # Flagged: red bold text
  if (any(mds_df$flagged))
    graphics::text(mds_df$MDS1[mds_df$flagged],
                   mds_df$MDS2[mds_df$flagged],
                   labels = participant[mds_df$flagged],
                   col = "red", cex = 0.65, font = 2)

  # No participant info: grey
  no_pid <- is.na(participant)
  if (any(no_pid))
    graphics::text(mds_df$MDS1[no_pid], mds_df$MDS2[no_pid],
                   labels = mds_df$sample_id[no_pid],
                   col = "grey50", cex = 0.45)

  # --- Plot 2: Flagged participants only (zoomed) ---
  if (nrow(flagged_df) > 0) {
    # Include ALL samples from any participant that has a flagged sample,
    # plus all samples from the participant they were misassigned to
    involved_pids <- unique(c(
      flagged_df$reported_participant,
      flagged_df$nearest_participant
    ))
    involved_idx <- which(participant %in% involved_pids)

    if (length(involved_idx) >= 2) {
      sub_df <- mds_df[involved_idx, ]
      sub_participant <- participant[involved_idx]

      # Compute axis limits with some padding
      x_range <- range(sub_df$MDS1)
      y_range <- range(sub_df$MDS2)
      x_pad <- diff(x_range) * 0.15
      y_pad <- diff(y_range) * 0.15

      plot(sub_df$MDS1, sub_df$MDS2, type = "n",
           xlim = x_range + c(-x_pad, x_pad),
           ylim = y_range + c(-y_pad, y_pad),
           xlab = "MDS1", ylab = "MDS2",
           main = sprintf("Flagged participants only (%d participants, %d samples)",
                          length(involved_pids), nrow(sub_df)),
           sub = sprintf("Showing all samples from participants with flagged samples"),
           cex.lab = 1.3, cex.main = 1.4)

      # Non-flagged samples from involved participants: black
      sub_ok <- !sub_df$flagged
      if (any(sub_ok))
        graphics::text(sub_df$MDS1[sub_ok], sub_df$MDS2[sub_ok],
                       labels = sub_participant[sub_ok],
                       col = "black", cex = 0.7)

      # Flagged samples: red bold
      sub_flagged <- sub_df$flagged
      if (any(sub_flagged))
        graphics::text(sub_df$MDS1[sub_flagged], sub_df$MDS2[sub_flagged],
                       labels = sub_participant[sub_flagged],
                       col = "red", cex = 0.8, font = 2)

      # Draw circles around participant centroids for context
      involved_centroids <- centroids[centroids$participant %in% involved_pids, ]
      for (r in seq_len(nrow(involved_centroids))) {
        graphics::points(involved_centroids$c_MDS1[r],
                         involved_centroids$c_MDS2[r],
                         pch = 3, cex = 2, col = "grey40", lwd = 2)
        graphics::text(involved_centroids$c_MDS1[r],
                       involved_centroids$c_MDS2[r],
                       labels = involved_centroids$participant[r],
                       pos = 1, cex = 0.6, col = "grey40", font = 3)
      }
    }
  }

  grDevices::dev.off()

  if (!is.null(logger))
    logger$log("snp_identity",
               sprintf("wrote snp_identity_mds.pdf (%d pages)",
                       if (nrow(flagged_df) > 0) 2 else 1))

  invisible(list(mds_coords = mds_df, flagged = flagged_df,
                 participant_col = participant_col))
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