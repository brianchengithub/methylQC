## ---------------------------------------------------------------------------
## clocks.R -- epigenetic age.
##
## Goal:    predict Horvath (2013) epigenetic age for every sample.
## Approach: apply the package's own bundled 353-probe model to unmasked betas,
##          imputing any absent probe from the cohort itself rather than from a
##          constant, and report per-sample coverage.
## Inputs:  the unmasked beta matrix.
## Outputs: a data.frame of sample_id, age and coverage columns.
## Usage:   ages <- predictages(betas)
##
## TWO DECISIONS WORTH STATING.
##
## 1. NO MASKING IS APPLIED. The published clocks were trained on beta values
##    that had neither a design mask nor a detection mask applied, so applying
##    one here would feed the model something different from what it was fitted
##    on. betas_all.rds is stored unmasked and is passed straight through.
##
## 2. THE MODEL IS BUNDLED, NOT FETCHED. methylQC 3.0.0 delegated to
##    sesame::predictAge() with a model from sesameData under the keys
##    "Clock_Horvath353" and "Anno/HM450/Clock_Horvath353.rds". Neither key
##    exists -- the only clock in sesameData is MM285.clock347, for the mouse
##    array -- so clockmodel() always returned NULL, predictages() always
##    returned NULL, and epigenetic age was never computed on any human array.
##    The coefficients ship with the package instead (horvath_clock.R, from
##    Horvath 2013 Additional File 3), which also makes results reproducible
##    across annotation releases.
## ---------------------------------------------------------------------------

## Coverage below which the prediction is reported but marked unreliable.
## Horvath's own transform is calibrated on the full 353; dropping a tenth of
## them shifts predictions by a year or more.
.MQC_CLOCK_MINCOVER <- 0.90

#' The bundled epigenetic clock model
#'
#' @param clock model name; only \code{"Horvath353"} is bundled.
#' @param logger optional logger.
#' @return a list with \code{probes}, \code{coef}, \code{intercept} and
#'   \code{defaults}, or \code{NULL} for an unknown name.
#' @export
#' @examples
#' m <- clockmodel("Horvath353")
#' length(m$probes)
clockmodel <- function(clock = "Horvath353", logger = NULL) {
  logger <- logger %||% nulllog()
  if (!identical(clock, "Horvath353")) {
    logger$log("clock", sprintf(
      "unknown clock '%s'; only 'Horvath353' is bundled", clock), warn = TRUE)
    return(NULL)
  }
  list(name = "Horvath353",
       probes = .horvath_probes,
       coef = .horvath_coefficients,
       intercept = .horvath_intercept,
       defaults = .horvath_zeroshot_betas)
}

## Horvath's inverse age transform, with adult age 20.
.horvath_untransform <- function(x) {
  if (!is.finite(x)) return(NA_real_)
  if (x < 0) 21 * exp(x) - 1 else 21 * x + 20
}

#' Predict epigenetic age for a cohort
#'
#' @param betas probes-by-samples beta matrix, unmasked.
#' @param clock clock name.
#' @param model an explicit model; overrides \code{clock}.
#' @param logger optional logger.
#' @return a data.frame with \code{sample_id}, the predicted age, the number of
#'   model probes present and the model size, or \code{NULL} if no model is
#'   usable.
#' @export
#' @examples
#' \dontrun{
#' ages <- predictages(betas)
#' }
predictages <- function(betas, clock = "Horvath353", model = NULL,
                        logger = NULL) {
  logger <- logger %||% nulllog()
  model <- model %||% clockmodel(clock, logger)
  if (is.null(model) || is.null(betas) || !ncol(betas)) return(NULL)

  ## EPICv2 identifiers carry replicate suffixes; the model keys on bare cg
  ## identifiers, so match on stripped names when the matrix has not been
  ## collapsed.
  rn <- rownames(betas)
  if (has_v2_suffix(rn)) rn <- stripv2(rn)

  want <- model$probes
  hit <- match(want, rn)
  present <- which(!is.na(hit))
  ncover <- length(present)
  logger$log("clock", sprintf("%d of %d %s probes are present on this array",
                              ncover, length(want), model$name))
  if (ncover < 10L) {
    logger$log("clock",
               "fewer than 10 model probes present; age prediction skipped",
               warn = TRUE)
    return(NULL)
  }
  if (ncover < length(want) * .MQC_CLOCK_MINCOVER)
    logger$log("clock", sprintf(
      paste("only %.1f%% of model probes are present; ages are computed on the",
            "reduced set and are likely biased"),
      100 * ncover / length(want)), warn = TRUE)

  sub <- betas[hit[present], , drop = FALSE]

  ## Missing values are filled from the cohort's own mean for that probe, and
  ## only from Horvath's published defaults when the probe is absent from the
  ## array entirely. 3.0.0 let sesame substitute zero while keeping the full
  ## intercept, which biases every prediction in the same direction.
  rmean <- rowMeans(sub, na.rm = TRUE)
  dflt <- model$defaults[want[present]]
  rmean[!is.finite(rmean)] <- unname(dflt[!is.finite(rmean)])
  rmean[!is.finite(rmean)] <- 0.5
  na_idx <- which(is.na(sub), arr.ind = TRUE)
  if (nrow(na_idx)) sub[na_idx] <- rmean[na_idx[, 1]]

  cf <- model$coef[present]
  x <- model$intercept + as.vector(crossprod(sub, cf))
  ages <- vapply(x, .horvath_untransform, numeric(1))

  key <- paste0("age_", tolower(model$name))
  out <- data.frame(sample_id = colnames(betas), stringsAsFactors = FALSE)
  out[[key]] <- ages
  out[[paste0(key, "_nprobes")]] <- ncover
  out[[paste0(key, "_nmodel")]] <- length(want)
  out
}
