## ---------------------------------------------------------------------------
## discover.R -- find IDAT pairs and reconcile them with a sample sheet.
## ---------------------------------------------------------------------------

#' Discover IDAT pairs beneath a directory
#'
#' Recursively locates complete Grn/Red IDAT pairs, reads any sample sheets it
#' finds alongside them, and returns one row per sample. Sample identifiers are
#' de-duplicated by prefixing the batch folder when the same Sentrix barcode
#' appears in more than one batch.
#'
#' @param indir directory to search.
#' @param recursive search subdirectories.
#' @param sheet optional explicit path to a sample sheet; when \code{NULL},
#'   sheets are auto-detected using the \code{sheetpattern} option.
#' @param logger optional logger.
#' @return a data.frame with at least \code{sample_id}, \code{Basename},
#'   \code{batch_folder} and \code{detected_platform}.
#' @export
#' @examples
#' \dontrun{
#' ss <- discover("~/idats")
#' }
discover <- function(indir, recursive = TRUE, sheet = NULL, logger = NULL) {
  logger <- logger %||% nulllog()
  cfg <- mqcopts()
  if (!dir.exists(indir)) stop("input directory does not exist: ", indir, call. = FALSE)

  grn <- list.files(indir, pattern = "_Grn\\.idat(\\.gz)?$",
                    recursive = recursive, full.names = TRUE, ignore.case = TRUE)
  if (!length(grn)) stop("no *_Grn.idat files found under ", indir, call. = FALSE)

  base <- sub("_Grn\\.idat(\\.gz)?$", "", grn, ignore.case = TRUE)
  red_ok <- vapply(base, function(b)
    any(file.exists(paste0(b, c("_Red.idat", "_Red.idat.gz")))), logical(1))

  if (any(!red_ok)) {
    logger$log("discover",
               sprintf("%d Grn file(s) have no matching Red file and were skipped",
                       sum(!red_ok)), warn = TRUE)
  }
  base <- base[red_ok]
  if (!length(base)) stop("no complete Grn/Red IDAT pairs found under ", indir,
                          call. = FALSE)

  batch <- basename(dirname(base))
  sentrix <- basename(base)

  ## De-duplicate: the same Sentrix barcode can legitimately appear in two
  ## batch folders. Prefix with the batch so downstream joins stay one-to-one.
  sid <- sentrix
  dup <- sentrix %in% sentrix[duplicated(sentrix)]
  if (any(dup)) {
    sid[dup] <- paste(batch[dup], sentrix[dup], sep = "_")
    logger$log("discover",
               sprintf("%d duplicated Sentrix ID(s) disambiguated by batch folder",
                       sum(dup)))
  }
  if (anyDuplicated(sid))
    stop("sample identifiers are still duplicated after batch prefixing; ",
         "check for repeated IDAT files.", call. = FALSE)

  out <- data.frame(sample_id = sid,
                    Basename = base,
                    sentrix = sentrix,
                    batch_folder = batch,
                    stringsAsFactors = FALSE)

  out$detected_platform <- .detect_platform(base, batch, logger)

  ss <- .read_sheets(indir, sheet, cfg, logger)
  if (!is.null(ss) && nrow(ss)) {
    key <- .resolve_col(ss, c(cfg$idcol, cfg$idaliases))
    if (!is.na(key)) {
      ss$.join <- as.character(ss[[key]])
      m <- match(out$sentrix, ss$.join)
      m2 <- match(out$sample_id, ss$.join)
      m[is.na(m)] <- m2[is.na(m)]
      extra <- setdiff(names(ss), c(".join", names(out)))
      for (cl in extra) out[[cl]] <- ss[[cl]][m]
      hit <- sum(!is.na(m))
      logger$log("discover",
                 sprintf("sample sheet matched %d of %d samples on column '%s'",
                         hit, nrow(out), key))
      if (hit == 0L)
        logger$log("discover",
                   "sample sheet matched NO samples; check the identifier column",
                   warn = TRUE)
    } else {
      logger$log("discover",
                 "no usable identifier column found in the sample sheet",
                 warn = TRUE)
    }
  }

  logger$log("discover", sprintf("%d samples across %d batch folder(s)",
                                 nrow(out), length(unique(out$batch_folder))))
  out
}

## Platform inference, one probe per batch folder.
##
## 3.0.0 called sesame::inferPlatformFromTango(), which is internal to sesame
## and not exported, so the call always errored into the fallback; and it
## inferred once and did rep(p, length(base)), which made every sample report
## the same platform by construction. qc()'s "this batch mixes array
## platforms" guard could therefore never fire. Probing one file per batch
## folder makes that guard real at the cost of one IDAT read per folder, which
## is where a genuine platform mix actually shows up.
.detect_platform <- function(base, batch, logger = NULL) {
  logger <- logger %||% nulllog()
  probe <- function(pfx) tryCatch({
    sdf <- sesame::readIDATpair(pfx)
    p <- attr(sdf, "platform")
    if (is.null(p) || !nzchar(p)) NA_character_ else as.character(p)
  }, error = function(e) NA_character_)

  folders <- unique(batch)
  first <- vapply(folders, function(f) base[which(batch == f)[1]], character(1))
  got <- vapply(unname(first), probe, character(1))
  names(got) <- folders

  out <- unname(got[batch])
  seen <- unique(stats::na.omit(out))
  if (!length(seen))
    logger$log("discover",
               "could not infer the array platform; pass platform= explicitly",
               warn = TRUE)
  else if (length(seen) > 1L)
    logger$log("discover", sprintf(
      "batch folders report more than one platform (%s); process each separately",
      paste(seen, collapse = ", ")), warn = TRUE)
  else
    logger$log("discover", sprintf("platform %s across %d batch folder(s)",
                                   seen, length(folders)))
  out
}

.read_sheets <- function(indir, sheet, cfg, logger) {
  files <- if (!is.null(sheet)) sheet else
    list.files(indir, pattern = cfg$sheetpattern, recursive = TRUE,
               full.names = TRUE)
  files <- files[file.exists(files)]
  if (!length(files)) return(NULL)

  frames <- lapply(files, function(f) {
    ext <- tolower(tools::file_ext(f))
    df <- tryCatch({
      if (ext %in% c("xlsx", "xls")) {
        if (!have_pkg("readxl", "reading Excel sample sheets", logger)) return(NULL)
        as.data.frame(readxl::read_excel(f), stringsAsFactors = FALSE)
      } else if (ext %in% c("tsv", "txt")) {
        utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
      } else {
        .read_illumina_csv(f)
      }
    }, error = function(e) {
      logger$log("discover", sprintf("could not read %s: %s", basename(f),
                                     conditionMessage(e)), warn = TRUE)
      NULL
    })
    df
  })
  frames <- Filter(function(d) !is.null(d) && nrow(d) > 0, frames)
  if (!length(frames)) return(NULL)

  common <- Reduce(intersect, lapply(frames, names))
  if (!length(common)) return(frames[[1]])
  do.call(rbind, lapply(frames, function(d) d[, common, drop = FALSE]))
}

## Illumina sample sheets carry a [Data] header block before the real table.
.read_illumina_csv <- function(f) {
  ln <- readLines(f, warn = FALSE)
  hit <- grep("^\\s*\\[Data\\]", ln)
  skip <- if (length(hit)) hit[1] else 0L
  utils::read.csv(f, skip = skip, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Find the first matching column name in a data frame
#'
#' @param df data frame
#' @param candidates candidate column names, in priority order
#' @return the matched name, or \code{NA_character_}
#' @keywords internal
#' @noRd
.resolve_col <- function(df, candidates) {
  nm <- names(df)
  for (c1 in candidates) {
    hit <- which(tolower(nm) == tolower(c1))
    if (length(hit)) return(nm[hit[1]])
  }
  NA_character_
}

#' Validate that a sample sheet has the columns downstream stages need
#'
#' @param ss sample sheet
#' @param logger optional logger
#' @return a list naming the resolved columns, invisibly
#' @export
#' @examples
#' \dontrun{
#' checkmeta(ss)
#' }
checkmeta <- function(ss, logger = NULL) {
  logger <- logger %||% nulllog()
  cfg <- mqcopts()
  got <- list(
    id    = .resolve_col(ss, c(cfg$idcol, cfg$idaliases)),
    sex   = .resolve_col(ss, c(cfg$sexcol, cfg$sexaliases)),
    age   = .resolve_col(ss, c(cfg$agecol, cfg$agealiases)),
    batch = .resolve_col(ss, c(cfg$batchcol, cfg$batchaliases)),
    cell  = .resolve_col(ss, c(cfg$cellcol, cfg$cellaliases)),
    donor = .resolve_col(ss, c(cfg$donorcol, cfg$donoraliases)))
  for (k in names(got)) {
    if (is.na(got[[k]]))
      logger$log("meta", sprintf("no %s column found; related checks are skipped", k))
  }
  invisible(got)
}
