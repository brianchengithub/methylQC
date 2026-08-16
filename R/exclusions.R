## ---------------------------------------------------------------------------
## exclusions.R -- sample-level flags: call rate, intensity, sex, age.
## ---------------------------------------------------------------------------

#' Normalise a reported sex column to M / F / NA
#'
#' @param x character or factor vector
#' @return a character vector of \code{"M"}, \code{"F"} or \code{NA}
#' @export
#' @examples
#' normsex(c("male", "F", "Female", "1", "unknown"))
normsex <- function(x) {
  s <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(s))
  out[s %in% c("M", "MALE", "BOY", "1")] <- "M"
  out[s %in% c("F", "FEMALE", "GIRL", "2")] <- "F"
  out
}

#' Parse an age column to numeric years
#'
#' @param x vector of ages
#' @return numeric vector; unparseable entries become \code{NA}
#' @export
#' @examples
#' parseage(c("42", "7.5 years", "unknown"))
parseage <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  suppressWarnings(as.numeric(gsub("[^0-9.\\-]", "", as.character(x))))
}

#' Flag samples failing quality thresholds
#'
#' Applies the call-rate, intensity, sex-concordance and age-outlier checks
#' and returns the sample sheet with flag columns added.
#'
#' Sex is called by \code{\link{sexcall}} from the sex-chromosome intensities
#' measured during Stage 1, not by \pkg{sesame}. When the cohort does not
#' support a call -- because it is single-sex, too small, or too noisy --
#' \code{inferred_sex} is \code{NA} for every sample and \code{sex_mismatch} is
#' \code{FALSE}, rather than a split being invented.
#'
#' @param qc per-sample metrics from \code{\link{qcmetrics}}.
#' @param sexcol resolved reported-sex column name, or \code{NA}.
#' @param agecol resolved age column name, or \code{NA}.
#' @param sexinfo a \code{\link{sexcall}} result, or \code{NULL} to derive one
#'   from the \code{sex_chrX_intensity} / \code{sex_chrY_intensity} columns.
#' @param agesd age outlier cutoff, in cohort standard deviations.
#' @param logger optional logger.
#' @return \code{qc} with added flag columns and a \code{flagged} /
#'   \code{flag_reason} summary.
#' @export
#' @examples
#' \dontrun{
#' flagged <- flagsamples(qc, sexcol = "Reported_Sex")
#' }
flagsamples <- function(qc, sexcol = NA_character_, agecol = NA_character_,
                        sexinfo = NULL, agesd = 3, logger = NULL) {

  logger <- logger %||% nulllog()
  out <- qc

  ## ---- sex ---------------------------------------------------------------
  out$reported_sex <- if (!is.na(sexcol) && sexcol %in% names(out))
    normsex(out[[sexcol]]) else NA_character_

  if (is.null(sexinfo) &&
      all(c("sex_chrX_intensity", "sex_chrY_intensity") %in% names(out)))
    sexinfo <- sexcall(out$sex_chrX_intensity, out$sex_chrY_intensity,
                       out$sample_id)

  if (is.null(sexinfo)) {
    out$inferred_sex <- NA_character_
    out$sex_unclear <- TRUE
    out$sex_confidence <- NA_real_
  } else {
    m <- match(out$sample_id, sexinfo$sample_id)
    out$inferred_sex <- sexinfo$inferred_sex[m]
    out$sex_unclear <- sexinfo$sex_unclear[m]
    out$sex_confidence <- sexinfo$sex_confidence[m]
    if (!isTRUE(attr(sexinfo, "called")))
      logger$log("sex", sprintf("no sex called: %s", attr(sexinfo, "reason")),
                 warn = TRUE)
    else
      logger$log("sex", sprintf(
        "sex called at chrY intensity %.0f (separation %.2f); %d F, %d M, %d unclear",
        attr(sexinfo, "threshold"), attr(sexinfo, "separation"),
        sum(out$inferred_sex == "F", na.rm = TRUE),
        sum(out$inferred_sex == "M", na.rm = TRUE),
        sum(out$sex_unclear, na.rm = TRUE)))
  }

  both <- !is.na(out$reported_sex) & !is.na(out$inferred_sex)
  out$sex_mismatch <- both & (out$reported_sex != out$inferred_sex)
  out$sex_unknown <- is.na(out$reported_sex) | is.na(out$inferred_sex)

  ## ---- age ---------------------------------------------------------------
  out$reported_age <- if (!is.na(agecol) && agecol %in% names(out))
    parseage(out[[agecol]]) else NA_real_
  out$age_outlier <- FALSE
  a <- out$reported_age
  if (sum(is.finite(a)) >= 5L) {
    mu <- mean(a, na.rm = TRUE); s <- stats::sd(a, na.rm = TRUE)
    if (is.finite(s) && s > 0)
      out$age_outlier <- is.finite(a) & abs(a - mu) > agesd * s
  }

  ## ---- combine -----------------------------------------------------------
  reasons <- list(
    low_callrate  = out$low_callrate,
    low_intensity = out$low_intensity,
    sex_mismatch  = out$sex_mismatch)
  for (k in names(reasons)) reasons[[k]][is.na(reasons[[k]])] <- FALSE

  out$flagged <- Reduce(`|`, reasons)
  out$flag_reason <- vapply(seq_len(nrow(out)), function(i) {
    hit <- names(reasons)[vapply(reasons, function(v) isTRUE(v[i]), logical(1))]
    if (!length(hit)) "" else paste(hit, collapse = ";")
  }, character(1))

  logger$log("flags", sprintf(
    "%d of %d sample(s) flagged (call rate %d, intensity %d, sex mismatch %d)",
    sum(out$flagged), nrow(out), sum(reasons$low_callrate),
    sum(reasons$low_intensity), sum(reasons$sex_mismatch)))
  if (any(out$age_outlier, na.rm = TRUE))
    logger$log("flags", sprintf(
      "%d age outlier(s) beyond %g SD (recorded, not used to exclude)",
      sum(out$age_outlier, na.rm = TRUE), agesd))
  out
}

#' Samples that passed every check
#'
#' @param flagged output of \code{\link{flagsamples}}
#' @return character vector of sample IDs
#' @export
#' @examples
#' \dontrun{
#' ids <- retained(flagged)
#' }
retained <- function(flagged) {
  keep <- !is.na(flagged$flagged) & !flagged$flagged
  flagged$sample_id[keep]
}
