###############################################################################
# discover_idats.R — IDAT file discovery and sample sheet parsing
#
# Recursively scans a directory tree for paired IDAT files (_Grn.idat /
# _Red.idat). For each subdirectory containing IDATs:
#   1. Searches for a sample sheet matching a configurable pattern
#   2. If found, reads it (auto-detecting tab vs. comma delimiters)
#   3. If not found, synthesizes a minimal sheet from IDAT filenames
#   4. Detects the array platform by reading the first IDAT header
#
# Sample sheet columns are resolved flexibly:
#   - Basename from: Basename, Fname, Sentrix_ID+Position, or Grn column
#   - Sex/Age/Cell via resolve_column() with configurable aliases
#
# After assembly, Sentrix IDs are decomposed into Chip/Row/Col for
# downstream batch association analysis in PCA.
###############################################################################
#' Discover IDAT files and optional sample sheets in nested folders
#'
#' Walks \code{dir} recursively, reads or synthesises a sample sheet per
#' folder, detects the array platform from IDAT headers, and decomposes
#' each basename into \code{Chip}, \code{Row}, and \code{Col} columns.
#'
#' @param dir Top-level directory.
#' @param sheetpattern Regex for sheet filenames; \code{NULL} to skip.
#'   Default: \code{sheetpattern} option.
#' @param basecol Preferred basename column name. Default: \code{basecol}
#'   option.
#' @param platform Optional platform string for a cross-folder consistency
#'   check.
#' @param logger Optional logger.
#' @return A data.frame with at least \code{Basename}, \code{sample_id},
#'   \code{batch_folder}, \code{sheet_path}, \code{detected_platform},
#'   \code{Chip}, \code{Row}, \code{Col}.
#' @export
discover <- function(dir,
                     sheetpattern = NULL,
                     basecol = NULL,
                     platform = NULL,
                     logger = NULL) {
  cfg <- mqcopts()
  if (is.null(sheetpattern)) sheetpattern <- cfg$sheetpattern
  if (is.null(basecol))      basecol      <- cfg$basecol
  stopifnot(dir.exists(dir))

  # Find all Green channel IDATs recursively
  all_grn <- list.files(dir, pattern = "_Grn\\.idat$",
                        recursive = TRUE, full.names = TRUE,
                        ignore.case = TRUE)
  if (!length(all_grn)) {
    stop("No IDAT files (*_Grn.idat) found under ", dir)
  }
  idat_dirs <- unique(dirname(all_grn))
  if (!is.null(logger)) {
    logger$log("discover",
               sprintf("found %d IDAT pair(s) across %d folder(s)",
                       length(all_grn), length(idat_dirs)))
  }

  # Process each folder independently: read/synthesize sheet + detect platform
  per_folder <- lapply(idat_dirs, function(folder) {
    df <- read_or_synthesize_sheet(folder, sheetpattern, basecol, logger)
    df$batch_folder <- basename(folder)
    detected <- detect_platform_from_folder(folder, logger)
    df$detected_platform <- detected
    if (!is.null(platform) && !is.na(detected) && detected != platform) {
      stop(sprintf("Platform mismatch in '%s': detected %s, expected %s",
                   folder, detected, platform))
    }
    df
  })

  detected_platforms <- unique(unlist(lapply(per_folder,
                                             function(d) d$detected_platform)))
  detected_platforms <- detected_platforms[!is.na(detected_platforms)]
  if (length(detected_platforms) > 1) {
    stop("Inconsistent platforms across folders: ",
         paste(detected_platforms, collapse = ", "))
  }
  if (!is.null(logger) && length(detected_platforms) == 1) {
    logger$log("discover",
               sprintf("detected platform: %s (consistent across all folders)",
                       detected_platforms))
  }

  # Harmonise columns across folders (sheets may differ in columns)
  all_cols <- unique(unlist(lapply(per_folder, colnames)))
  per_folder <- lapply(per_folder, function(d) {
    for (m in setdiff(all_cols, colnames(d))) d[[m]] <- NA
    d[, all_cols]
  })
  ss <- do.call(rbind, per_folder)

  # Verify that both Green and Red IDATs exist for each sample
  grn_ok <- file.exists(paste0(ss$Basename, "_Grn.idat"))
  red_ok <- file.exists(paste0(ss$Basename, "_Red.idat"))
  ok <- grn_ok & red_ok
  if (any(!ok) && !is.null(logger)) {
    logger$log("discover",
               sprintf("dropping %d row(s) with missing IDATs", sum(!ok)))
  }
  ss <- ss[ok, , drop = FALSE]

  # Disambiguate duplicate sample IDs across batches
  dups <- ss$sample_id[duplicated(ss$sample_id)]
  if (length(dups)) {
    ss$sample_id <- ifelse(ss$sample_id %in% dups,
                           paste(ss$batch_folder, ss$sample_id, sep = "__"),
                           ss$sample_id)
  }

  # Decompose basename into Chip/Row/Col
  ss <- decompose_basename(ss, logger = logger)

  if (!is.null(logger)) {
    by_batch <- table(ss$batch_folder)
    logger$log("discover",
               sprintf("final: %d samples across %d batches (%s)",
                       nrow(ss), length(by_batch),
                       paste(names(by_batch), by_batch,
                             sep = "=", collapse = "; ")))
  }
  ss
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

  n_decomposed <- sum(!is.na(ss$Chip))
  if (!is.null(logger)) {
    logger$log("discover",
               sprintf("decomposed Sentrix: %d / %d samples -> Chip/Row/Col",
                       n_decomposed, nrow(ss)))
  }
  ss
}

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

#' @keywords internal
#' @noRd
read_or_synthesize_sheet <- function(folder, pattern, basecol, logger) {
  sheet_path <- find_sample_sheet(folder, pattern)
  if (is.null(sheet_path)) {
    grn <- list.files(folder, pattern = "_Grn\\.idat$",
                      ignore.case = TRUE, full.names = TRUE)
    bn  <- sub("_Grn\\.idat$", "", grn, ignore.case = TRUE)
    if (!is.null(logger)) {
      logger$log("discover",
                 sprintf("no sheet in %s -- synthesising from %d IDAT pairs",
                         basename(folder), length(bn)))
    }
    return(data.frame(Basename = bn, sample_id = basename(bn),
                      sheet_path = NA_character_,
                      stringsAsFactors = FALSE))
  }
  df <- read_sheet_autodelim(sheet_path)

  # Resolve the IDAT basename from whichever column is available
  if (basecol %in% colnames(df)) {
    bn <- df[[basecol]]
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
    stop("Sample sheet ", sheet_path, " lacks any recognised basename column.")
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

#' @keywords internal
#' @noRd
find_sample_sheet <- function(folder, pattern) {
  if (is.null(pattern)) return(NULL)
  hits <- list.files(folder, pattern = pattern,
                     ignore.case = TRUE, full.names = TRUE)
  hits <- hits[dirname(hits) == folder]
  if (length(hits) == 0) return(NULL)
  if (length(hits) > 1) {
    stop("Multiple files matching sample sheet pattern in ", folder, ": ",
         paste(basename(hits), collapse = ", "))
  }
  hits[1]
}
