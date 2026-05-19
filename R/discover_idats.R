###############################################################################
# discover_idats.R — IDAT file discovery and sample sheet parsing
#
# The user supplies a single top-level directory. All IDAT discovery
# and all sample-sheet discovery are confined to that directory tree.
#
# Workflow:
#   1. Recursively find every IDAT pair (_Grn.idat / _Red.idat) under
#      the root.
#   2. Recursively find every sample sheet under the root and read them
#      all (auto-detecting tab vs. comma delimiters).
#   3. Concatenate the sheets, taking the union of their columns and
#      filling absent columns with NA.
#   4. Reconcile sheet rows against the IDATs actually on disk:
#        - IDATs present on disk but absent from every sheet: a minimal
#          row is synthesized and appended (not written back to disk).
#        - Rows in a sheet whose IDATs are missing on disk: a WARNING is
#          emitted and the row is dropped.
#      A mismatch is always a warning, never an error.
#   5. Resolve duplicate sample IDs: when one sample_id appears in more
#      than one sheet row, the row with the most non-missing values is
#      kept (ties -> first occurrence) and a WARNING naming the
#      sample_id and the sheet file(s) is emitted.
#   6. Detect the array platform from IDAT headers and decompose the
#      Sentrix ID into Chip/Row/Col.
###############################################################################

#' Discover IDAT files and sample sheets under a single root directory
#'
#' @param idat_root Top-level directory. IDATs and sample sheets are
#'   searched for recursively, confined to this tree.
#' @param sample_sheet_pattern Regex for sheet filenames; NULL to skip
#'   sheet discovery entirely (all rows synthesized).
#' @param basename_col Preferred basename column name.
#' @param expected_platform Optional platform string for a cross-folder
#'   sanity check (mismatch is a warning).
#' @param logger Optional logger.
#' @return A data.frame with at least Basename, sample_id, batch_folder,
#'   sheet_path, detected_platform, Chip, Row, Col.
#' @export
discover_idats <- function(idat_root,
                           sample_sheet_pattern = NULL,
                           basename_col = NULL,
                           expected_platform = NULL,
                           logger = NULL) {
  cfg <- methylQC_options()
  if (is.null(sample_sheet_pattern)) sample_sheet_pattern <- cfg$sample_sheet_pattern
  if (is.null(basename_col))         basename_col         <- cfg$basename_col
  stopifnot(dir.exists(idat_root))
  loginfo <- function(...) if (!is.null(logger)) logger$log("discover", sprintf(...))

  # --- 1. Find every IDAT pair under the root ---
  all_grn <- list.files(idat_root, pattern = "_Grn\\.idat$",
                        recursive = TRUE, full.names = TRUE,
                        ignore.case = TRUE)
  if (!length(all_grn)) {
    stop("No IDAT files (*_Grn.idat) found under ", idat_root)
  }
  all_bn <- sub("_Grn\\.idat$", "", all_grn, ignore.case = TRUE)
  red_exists <- file.exists(paste0(all_bn, "_Red.idat"))
  if (any(!red_exists)) {
    loginfo("%d IDAT(s) have a Grn file but no matching Red file - skipped",
            sum(!red_exists))
    all_bn <- all_bn[red_exists]
  }
  disk <- data.frame(
    Basename     = all_bn,
    sample_id    = basename(all_bn),
    batch_folder = basename(dirname(all_bn)),
    stringsAsFactors = FALSE)
  idat_dirs <- unique(dirname(all_bn))
  loginfo("found %d IDAT pair(s) across %d folder(s)",
          nrow(disk), length(idat_dirs))

  # --- 2. Find and read every sample sheet under the root ---
  sheets <- list()
  if (!is.null(sample_sheet_pattern)) {
    sheet_paths <- list.files(idat_root, pattern = sample_sheet_pattern,
                              recursive = TRUE, full.names = TRUE,
                              ignore.case = TRUE)
    for (sp in sheet_paths) {
      df <- tryCatch(read_one_sheet(sp, basename_col, logger),
                     error = function(e) {
                       warning(sprintf("Could not read sample sheet '%s': %s",
                                        sp, conditionMessage(e)), call. = FALSE)
                       NULL
                     })
      if (!is.null(df) && nrow(df) > 0) sheets[[sp]] <- df
    }
    loginfo("read %d sample sheet(s)", length(sheets))
  }

  # --- 3. Concatenate sheets with a column union ---
  sheet_all <- if (length(sheets)) {
    all_cols <- unique(unlist(lapply(sheets, colnames)))
    sheets <- lapply(sheets, function(d) {
      for (m in setdiff(all_cols, colnames(d))) d[[m]] <- NA
      d[, all_cols, drop = FALSE]
    })
    do.call(rbind, sheets)
  } else {
    data.frame(Basename = character(0), sample_id = character(0),
               sheet_path = character(0), stringsAsFactors = FALSE)
  }

  # --- 4. Reconcile sheet rows against IDATs on disk ---
  # 4a. Sheet rows whose IDATs are missing on disk -> warn + drop
  if (nrow(sheet_all) > 0) {
    on_disk <- sheet_all$sample_id %in% disk$sample_id
    if (any(!on_disk)) {
      missing_ids <- sheet_all$sample_id[!on_disk]
      warning(sprintf(paste0(
        "%d sample(s) listed in sample sheet(s) have no IDAT files on ",
        "disk and were dropped: %s"),
        length(missing_ids),
        paste(utils::head(missing_ids, 20), collapse = ", ")),
        call. = FALSE)
      loginfo("%d sheet row(s) dropped: IDATs not found on disk",
              sum(!on_disk))
      sheet_all <- sheet_all[on_disk, , drop = FALSE]
    }
  }

  # 4b. Resolve duplicate sample IDs within the sheets
  if (nrow(sheet_all) > 0 && anyDuplicated(sheet_all$sample_id)) {
    sheet_all <- resolve_duplicate_rows(sheet_all, logger)
  }

  # 4c. IDATs on disk absent from every sheet -> synthesize rows
  matched <- disk$sample_id %in% sheet_all$sample_id
  ss <- sheet_all
  if (any(!matched)) {
    extra <- disk[!matched, , drop = FALSE]
    extra$sheet_path <- NA_character_
    loginfo("%d IDAT(s) not listed in any sheet - synthesizing rows",
            nrow(extra))
    if (nrow(ss) > 0) {
      for (m in setdiff(colnames(ss), colnames(extra))) extra[[m]] <- NA
      for (m in setdiff(colnames(extra), colnames(ss)))  ss[[m]]    <- NA
      extra <- extra[, colnames(ss), drop = FALSE]
      ss <- rbind(ss, extra)
    } else {
      ss <- extra
    }
  }

  # Attach authoritative Basename / batch_folder from disk
  di <- match(ss$sample_id, disk$sample_id)
  ss$Basename     <- disk$Basename[di]
  ss$batch_folder <- disk$batch_folder[di]
  if (!"sheet_path" %in% colnames(ss)) ss$sheet_path <- NA_character_

  # --- 5. Detect platform per folder ---
  plat_by_dir <- vapply(idat_dirs, function(d)
    detect_platform_from_folder(d, logger), character(1))
  names(plat_by_dir) <- basename(idat_dirs)
  ss$detected_platform <- plat_by_dir[ss$batch_folder]

  detected <- unique(stats::na.omit(ss$detected_platform))
  if (length(detected) > 1) {
    warning("Inconsistent platforms detected across folders: ",
            paste(detected, collapse = ", "),
            ". Proceeding with the most common one.", call. = FALSE)
    detected <- names(sort(table(ss$detected_platform), decreasing = TRUE))[1]
  }
  if (length(detected) == 1) {
    loginfo("detected platform: %s", detected)
    if (!is.null(expected_platform) && detected != expected_platform) {
      warning(sprintf("Detected platform '%s' differs from expected '%s'.",
                      detected, expected_platform), call. = FALSE)
    }
  }

  # --- 6. Decompose Sentrix ID into Chip/Row/Col ---
  ss <- decompose_basename(ss, logger = logger)

  by_batch <- table(ss$batch_folder)
  loginfo("final: %d samples across %d batch(es) (%s)",
          nrow(ss), length(by_batch),
          paste(names(by_batch), by_batch, sep = "=", collapse = "; "))
  ss
}

#' Resolve duplicate sample IDs by keeping the most complete row
#'
#' When a sample_id appears more than once, keep the row with the most
#' non-missing values (ties -> first occurrence). Warn, naming the
#' sample_id and the source sheet file(s).
#' @keywords internal
#' @noRd
resolve_duplicate_rows <- function(sheet_all, logger = NULL) {
  dup_ids <- unique(sheet_all$sample_id[duplicated(sheet_all$sample_id)])
  keep <- rep(TRUE, nrow(sheet_all))
  for (sid in dup_ids) {
    idx <- which(sheet_all$sample_id == sid)
    completeness <- vapply(idx, function(i)
      sum(!is.na(unlist(sheet_all[i, , drop = TRUE]))), integer(1))
    winner <- idx[which.max(completeness)]   # which.max -> first on ties
    losers <- setdiff(idx, winner)
    keep[losers] <- FALSE
    sheet_files <- unique(stats::na.omit(sheet_all$sheet_path[idx]))
    warning(sprintf(paste0(
      "Duplicate sample_id '%s' found in %d sheet rows (%s); kept the ",
      "most complete row, dropped %d."),
      sid, length(idx),
      paste(basename(sheet_files), collapse = ", "), length(losers)),
      call. = FALSE)
  }
  if (!is.null(logger)) {
    logger$log("discover",
               sprintf("resolved %d duplicate sample_id(s)", length(dup_ids)))
  }
  sheet_all[keep, , drop = FALSE]
}

#' Read a single sample sheet and resolve its Basename column
#' @keywords internal
#' @noRd
read_one_sheet <- function(sheet_path, basename_col, logger = NULL) {
  df <- read_sheet_autodelim(sheet_path)
  folder <- dirname(sheet_path)

  if (basename_col %in% colnames(df)) {
    bn <- df[[basename_col]]
    df$Basename <- ifelse(startsWith(bn, "/") | grepl("^[A-Za-z]:", bn),
                          bn, file.path(folder, bn))
  } else if ("Fname" %in% colnames(df)) {
    df$Basename <- file.path(folder, df$Fname)
  } else if (all(c("Sentrix_ID", "Sentrix_Position") %in% colnames(df))) {
    df$Basename <- file.path(folder,
                             paste0(df$Sentrix_ID, "_", df$Sentrix_Position))
  } else if ("Grn" %in% colnames(df)) {
    bn <- sub("_Grn\\.idat$", "", df$Grn, ignore.case = TRUE)
    df$Basename <- file.path(folder, bn)
  } else {
    stop("no recognized basename column (Basename / Fname / ",
         "Sentrix_ID+Position / Grn)")
  }

  df$sample_id  <- basename(df$Basename)
  df$sheet_path <- sheet_path
  if (!is.null(logger)) {
    logger$log("discover",
               sprintf("read sheet %s (%d rows)",
                       basename(sheet_path), nrow(df)))
  }
  df
}

#' Decompose sample_id (Sentrix) into Chip/Row/Col
#' @keywords internal
#' @noRd
decompose_basename <- function(ss, logger = NULL) {
  pattern <- "^(.+)_R(\\d+)C(\\d+)$"
  m <- regexec(pattern, ss$sample_id)
  parts <- regmatches(ss$sample_id, m)

  ss$Chip <- vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_,
                    character(1))
  ss$Row  <- vapply(parts, function(p) if (length(p) >= 3) p[3] else NA_character_,
                    character(1))
  ss$Col  <- vapply(parts, function(p) if (length(p) >= 4) p[4] else NA_character_,
                    character(1))

  if (!is.null(logger)) {
    logger$log("discover",
               sprintf("decomposed Sentrix: %d / %d samples -> Chip/Row/Col",
                       sum(!is.na(ss$Chip)), nrow(ss)))
  }
  ss
}

#' Detect array platform by reading the first IDAT header in a folder
#' @keywords internal
#' @noRd
detect_platform_from_folder <- function(folder, logger = NULL) {
  grn_files <- list.files(folder, pattern = "_Grn\\.idat$",
                          ignore.case = TRUE, full.names = TRUE)
  if (length(grn_files) == 0) return(NA_character_)
  basename_path <- sub("_Grn\\.idat$", "", grn_files[1], ignore.case = TRUE)

  platform <- tryCatch({
    sdf <- sesame::readIDATpair(basename_path)
    p <- attr(sdf, "platform")
    if (is.null(p) || is.na(p)) {
      if (exists("sdfPlatform", where = asNamespace("sesame"))) {
        p <- sesame::sdfPlatform(sdf)
      }
    }
    rm(sdf); gc(verbose = FALSE)
    if (is.null(p)) NA_character_ else as.character(p)
  }, error = function(e) {
    if (!is.null(logger)) {
      logger$log("discover",
                 sprintf("platform detection failed for %s: %s",
                         basename(grn_files[1]), conditionMessage(e)))
    }
    NA_character_
  })

  if (!is.null(logger) && !is.na(platform)) {
    logger$log("discover",
               sprintf("folder %s -> %s", basename(folder), platform))
  }
  platform
}

#' Auto-detect tab vs. comma delimiter and read a sample sheet
#' @keywords internal
#' @noRd
read_sheet_autodelim <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE)
  n_tabs   <- nchar(gsub("[^\t]", "", first_line))
  n_commas <- nchar(gsub("[^,]", "", first_line))
  sep <- if (n_tabs > n_commas) "\t" else ","
  utils::read.table(path, header = TRUE, sep = sep,
                    stringsAsFactors = FALSE,
                    check.names = FALSE, na.strings = c("", "NA"))
}
