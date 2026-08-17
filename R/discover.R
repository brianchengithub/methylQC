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
  if (!is.null(ss) && nrow(ss)) out <- .join_sheet(out, ss, cfg, logger)

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
  if (!length(files)) {
    logger$log("discover", paste(
      "no sample sheet found. Looked for files matching the 'sheetpattern'",
      "option; pass sheet= explicitly if yours is named differently."),
      warn = TRUE)
    return(NULL)
  }

  frames <- lapply(files, function(f) {
    df <- tryCatch(.read_table(f, logger), error = function(e) {
      logger$log("discover", sprintf("could not read %s: %s", basename(f),
                                     conditionMessage(e)), warn = TRUE)
      NULL
    })
    if (is.null(df) || !nrow(df)) return(NULL)
    ## A broad filename pattern will pick up unrelated text files, so keep only
    ## tables that actually look like a sample annotation.
    if (is.na(.resolve_col(df, c(cfg$idcol, cfg$idaliases)))) {
      logger$log("discover", sprintf(
        "ignoring %s: no recognisable sample identifier column (%s)",
        basename(f), paste(utils::head(names(df), 6), collapse = ", ")))
      return(NULL)
    }
    logger$log("discover", sprintf("read %s: %d row(s), %d column(s)",
                                   basename(f), nrow(df), ncol(df)))
    df
  })
  frames <- Filter(Negate(is.null), frames)
  if (!length(frames)) return(NULL)
  if (length(frames) == 1L) return(frames[[1]])

  common <- Reduce(intersect, lapply(frames, names))
  if (!length(common)) return(frames[[1]])
  do.call(rbind, lapply(frames, function(d) d[, common, drop = FALSE]))
}

## Read one annotation table, whatever shape it arrives in: Excel, or delimited
## text with the separator inferred rather than assumed. 3.0.x sent every .txt
## through read.delim(), so a space- or semicolon-separated file parsed as a
## single column and was silently useless.
.read_table <- function(f, logger = NULL) {
  logger <- logger %||% nulllog()
  ext <- tolower(tools::file_ext(f))
  if (ext %in% c("xlsx", "xls")) {
    if (!have_pkg("readxl", "reading Excel sample sheets", logger)) return(NULL)
    return(as.data.frame(readxl::read_excel(f), stringsAsFactors = FALSE))
  }

  ln <- readLines(f, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  if (!length(ln)) return(NULL)

  ## Illumina sheets carry a [Data] block before the real table.
  hit <- grep("^\\s*\\[Data\\]", ln)
  skip <- if (length(hit)) hit[1] else 0L
  hdr <- ln[skip + 1L]

  sep <- .guess_sep(hdr)
  utils::read.table(f, sep = sep, header = TRUE, skip = skip,
                    stringsAsFactors = FALSE, check.names = FALSE,
                    quote = "\"'", comment.char = "", fill = TRUE,
                    na.strings = c("NA", "", "NaN", "n/a", "N/A", "#N/A"))
}

## Pick the delimiter that splits the header into the most fields.
.guess_sep <- function(header) {
  cand <- c("\t", ",", ";", "|")
  n <- vapply(cand, function(s) length(strsplit(header, s, fixed = TRUE)[[1]]),
              integer(1))
  if (max(n) > 1L) return(cand[which.max(n)])
  if (length(strsplit(trimws(header), "[[:space:]]+")[[1]]) > 1L) return("")
  "\t"
}

## Attach the sheet's columns to the discovered samples.
##
## The identifier in a lab-made sheet is rarely the Sentrix barcode: it may be
## the IDAT file name, a path, or Sentrix_ID and Sentrix_Position in separate
## columns. Every plausible key is tried and the one that matches most samples
## wins, so the user does not have to know which convention their sheet uses.
.join_sheet <- function(out, ss, cfg, logger) {
  key <- .resolve_col(ss, c(cfg$idcol, cfg$idaliases))
  keys <- list()
  if (!is.na(key)) keys[[key]] <- as.character(ss[[key]])

  ## Sentrix_ID + Sentrix_Position, the Illumina convention.
  sid <- .resolve_col(ss, c("Sentrix_ID", "SentrixID", "Slide", "Sentrix_Barcode"))
  spo <- .resolve_col(ss, c("Sentrix_Position", "SentrixPosition", "Array", "Well"))
  if (!is.na(sid) && !is.na(spo))
    keys[["Sentrix_ID+Position"]] <- paste0(ss[[sid]], "_", ss[[spo]])

  if (!length(keys)) {
    logger$log("discover",
               "no usable identifier column found in the sample sheet",
               warn = TRUE)
    return(out)
  }

  ## Candidate targets on our side, and a normaliser that strips directories,
  ## IDAT suffixes and case so "…/200607130026_R06C01_Grn.idat" matches
  ## "200607130026_R06C01".
  norm <- function(x) {
    x <- basename(as.character(x))
    x <- sub("_(Grn|Red)\\.idat(\\.gz)?$", "", x, ignore.case = TRUE)
    tolower(trimws(x))
  }
  targets <- list(sentrix = norm(out$sentrix),
                  sample_id = norm(out$sample_id),
                  basename = norm(out$Basename))

  best <- list(n = -1L, m = NULL, key = NA_character_, target = NA_character_)
  for (kn in names(keys)) {
    kv <- norm(keys[[kn]])
    for (tn in names(targets)) {
      m <- match(targets[[tn]], kv)
      n <- sum(!is.na(m))
      if (n > best$n) best <- list(n = n, m = m, key = kn, target = tn)
    }
  }

  if (best$n == 0L) {
    logger$log("discover", sprintf(paste(
      "sample sheet matched NO samples. Tried column(s) %s against the Sentrix",
      "barcode, the sample id and the IDAT file name. Sheet ids look like '%s';",
      "IDAT names look like '%s'."),
      paste(names(keys), collapse = ", "),
      as.character(keys[[1]])[1], out$sentrix[1]), warn = TRUE)
    return(out)
  }

  drop <- c(names(out), sid, spo)
  extra <- setdiff(names(ss), stats::na.omit(drop))
  for (cl in extra) out[[cl]] <- ss[[cl]][best$m]
  logger$log("discover", sprintf(
    "sample sheet matched %d of %d samples (sheet column '%s' vs %s); added: %s",
    best$n, nrow(out), best$key, best$target,
    if (length(extra)) paste(extra, collapse = ", ") else "no new columns"))
  if (best$n < nrow(out))
    logger$log("discover", sprintf(
      "%d sample(s) have no sample sheet row and will carry NA metadata",
      nrow(out) - best$n), warn = TRUE)
  out
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
