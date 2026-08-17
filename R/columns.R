## ---------------------------------------------------------------------------
## columns.R -- work out which sample sheet column is which.
##
## Goal:    resolve the id / sex / age / donor / batch / cell columns of an
##          arbitrary sample sheet, without requiring the user to rename
##          anything.
## Approach: match on name first, because a name match is predictable and fast,
##          then CONFIRM it against the column's actual values. A confirmed
##          name match wins. A name match whose contents contradict it is
##          rejected. When no name matches, the values alone are used to find a
##          candidate, so a sheet with idiosyncratic headers still resolves.
## Inputs:  a data.frame, plus the alias lists from mqcopts().
## Outputs: a named list of column names, or NA where nothing was found.
## Usage:   cols <- checkmeta(ss)
##
## Name matching alone is fragile in both directions. It misses a column called
## something the alias list never anticipated, and -- worse, because it is
## silent -- it accepts a column whose name looks right but whose contents are
## not what the name implies, such as a "Sex" column holding a chromosome count
## or a free-text note.
## ---------------------------------------------------------------------------

## How many rows to inspect. Enough to characterise a column, cheap on a sheet
## of any size.
.MQC_PROFILE_N <- 200L

## Fraction of non-missing values that must fit the expected pattern.
.MQC_PROFILE_FRAC <- 0.80

## A Sentrix barcode with its position, e.g. 200607130026_R06C01.
.MQC_SENTRIX_RE <- "^[0-9]{8,14}_R[0-9]{2}C[0-9]{2}$"

#' Sample the leading non-missing values of a column
#' @keywords internal
#' @noRd
.peek <- function(x, n = .MQC_PROFILE_N) {
  x <- x[!is.na(x)]
  if (!length(x)) return(character(0))
  as.character(utils::head(x, n))
}

## ---------------------------------------------------------------------------
## Content tests. Each answers: could this column be a <kind> column?
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.looks_sex <- function(x) {
  v <- .peek(x)
  if (length(v) < 2L) return(FALSE)
  u <- unique(toupper(trimws(v)))
  ## A sex column has few distinct values and nearly all of them normalise.
  if (length(u) > 6L) return(FALSE)
  mean(!is.na(normsex(v))) >= .MQC_PROFILE_FRAC
}

#' @keywords internal
#' @noRd
.looks_age <- function(x) {
  v <- .peek(x)
  if (length(v) < 3L) return(FALSE)
  a <- suppressWarnings(parseage(v))
  ok <- is.finite(a)
  if (mean(ok) < .MQC_PROFILE_FRAC) return(FALSE)
  a <- a[ok]
  ## Human ages: plausible range, and not a constant or an identifier. A
  ## sentrix barcode parses as a number too, hence the upper bound.
  length(unique(a)) > 1L && stats::median(a) >= 0 && max(a) <= 130
}

#' @keywords internal
#' @noRd
.looks_id <- function(x, targets = NULL) {
  v <- .peek(x, .Machine$integer.max)
  if (length(v) < 2L) return(FALSE)
  if (anyDuplicated(v)) return(FALSE)             # an identifier is unique
  if (!is.null(targets) && length(targets))
    return(mean(.normid(v) %in% .normid(targets)) >= 0.5)
  ## With nothing to match against, only a Sentrix barcode is self-evidently an
  ## identifier. Accepting any unique non-numeric column here would let a
  ## two-line note file masquerade as a sample sheet.
  mean(grepl(.MQC_SENTRIX_RE, v)) >= .MQC_PROFILE_FRAC
}

#' @keywords internal
#' @noRd
.looks_grouping <- function(x, maxlev = 50L) {
  v <- .peek(x)
  if (length(v) < 2L) return(FALSE)
  nu <- length(unique(v))
  nu >= 1L && nu <= min(maxlev, length(v))
}

## Normalise an identifier for comparison: drop directories, IDAT suffixes and
## case, so "…/200607130026_R06C01_Grn.idat" and "200607130026_R06C01" agree.
#' @keywords internal
#' @noRd
.normid <- function(x) {
  x <- basename(as.character(x))
  x <- sub("_(Grn|Red)\\.idat(\\.gz)?$", "", x, ignore.case = TRUE)
  tolower(trimws(x))
}

.MQC_COL_TESTS <- list(
  id    = function(x, targets) .looks_id(x, targets),
  sex   = function(x, targets) .looks_sex(x),
  age   = function(x, targets) .looks_age(x),
  donor = function(x, targets) .looks_grouping(x, maxlev = .Machine$integer.max),
  batch = function(x, targets) .looks_grouping(x, maxlev = 200L),
  cell  = function(x, targets) .looks_grouping(x, maxlev = 20L))

#' Resolve one column by name, confirmed against its contents
#'
#' @param df the sheet.
#' @param aliases candidate names, in priority order.
#' @param kind one of the names of \code{.MQC_COL_TESTS}.
#' @param targets identifiers to match against, for \code{kind = "id"}.
#' @param logger optional logger.
#' @return a list with \code{col} and \code{how} (\code{"name"},
#'   \code{"content"}, \code{"name-unconfirmed"} or \code{NA}).
#' @keywords internal
#' @noRd
.pick_col <- function(df, aliases, kind, targets = NULL, logger = NULL) {
  logger <- logger %||% nulllog()
  test <- .MQC_COL_TESTS[[kind]]
  nm <- names(df)
  none <- list(col = NA_character_, how = NA_character_)
  if (!length(nm)) return(none)

  named <- nm[tolower(nm) %in% tolower(aliases)]
  ## Preserve the caller's priority order among the names that matched.
  named <- named[order(match(tolower(named), tolower(aliases)))]

  for (k in named) {
    if (isTRUE(test(df[[k]], targets)))
      return(list(col = k, how = "name"))
  }

  ## No name matched, or every one that did was contradicted by its contents.
  ## Fall back to whichever remaining column fits the pattern.
  for (k in setdiff(nm, named)) {
    if (isTRUE(test(df[[k]], targets))) {
      logger$log("meta", sprintf(
        "no recognised %s column name; using '%s' because its values look like %s",
        kind, k, kind))
      return(list(col = k, how = "content"))
    }
  }

  if (length(named)) {
    ## A name matched but nothing confirmed it. Say so and use it anyway --
    ## silently discarding the user's own label would be worse -- but the log
    ## records that the contents did not agree.
    logger$log("meta", sprintf(
      "column '%s' is named like a %s column but its values do not look like %s; using it anyway",
      named[1], kind, kind), warn = TRUE)
    return(list(col = named[1], how = "name-unconfirmed"))
  }
  none
}

#' Validate that a sample sheet has the columns downstream stages need
#'
#' Each column is matched on name and then confirmed against its own values,
#' so a sheet with unexpected headers still resolves and a misleading header is
#' not taken at face value.
#'
#' @param ss sample sheet.
#' @param logger optional logger.
#' @return a list naming the resolved columns, invisibly.
#' @export
#' @examples
#' ss <- data.frame(Fname = c("a", "b", "c"), Age = c(30, 40, 50),
#'                  Sex = c("M", "F", "F"), Donor = c("d1", "d2", "d3"))
#' checkmeta(ss)
checkmeta <- function(ss, logger = NULL) {
  logger <- logger %||% nulllog()
  cfg <- mqcopts()

  ## Discovery columns are the ground truth for what an identifier looks like.
  targets <- unique(c(if ("sentrix" %in% names(ss)) as.character(ss$sentrix),
                      if ("sample_id" %in% names(ss)) as.character(ss$sample_id)))

  spec <- list(
    id    = list(cfg$idcol,    cfg$idaliases,    targets),
    sex   = list(cfg$sexcol,   cfg$sexaliases,   NULL),
    age   = list(cfg$agecol,   cfg$agealiases,   NULL),
    batch = list(cfg$batchcol, cfg$batchaliases, NULL),
    cell  = list(cfg$cellcol,  cfg$cellaliases,  NULL),
    donor = list(cfg$donorcol, cfg$donoraliases, NULL))

  got <- list(); how <- list()
  ## Once a column is claimed it is not offered to a later kind, so a single
  ## column cannot be reported as both the donor and the batch.
  claimed <- character(0)
  for (k in names(spec)) {
    s <- spec[[k]]
    sub <- ss[, setdiff(names(ss), claimed), drop = FALSE]
    r <- .pick_col(sub, c(s[[1]], s[[2]]), k, targets = s[[3]], logger = logger)
    got[[k]] <- r$col
    how[[k]] <- r$how
    if (!is.na(r$col)) claimed <- c(claimed, r$col)
  }

  for (k in names(got)) {
    if (is.na(got[[k]]))
      logger$log("meta", sprintf("no %s column found; related checks are skipped", k))
    else if (identical(how[[k]], "name"))
      logger$log("meta", sprintf("%s column: '%s'", k, got[[k]]))
  }
  attr(got, "how") <- unlist(how)
  invisible(got)
}
