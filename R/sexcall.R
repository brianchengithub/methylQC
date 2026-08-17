## ---------------------------------------------------------------------------
## sexcall.R -- cohort-relative sex calling from sex-chromosome intensity.
##
## Goal:    call sex per sample without depending on sesame::inferSex(), and
##          without inventing a split when the cohort has only one sex in it.
## Method:  per-sample median total intensity over curated chrX-x-linked and
##          chrY-clean probe sets (measured in the worker, see preprocess.R),
##          then a cohort-level two-cluster split on log2 chrY intensity.
## Inputs:  the sex_chrX_intensity / sex_chrY_intensity columns of the Stage 1
##          diagnostics.
## Outputs: inferred_sex, sex_unclear, sex_confidence, plus the chosen
##          threshold and separation as attributes.
## Usage:   sexcall(qc$sex_chrX_intensity, qc$sex_chrY_intensity, qc$sample_id)
##
## WHY NOT sesame::inferSex(): its signature is inferSex(betas, platform) and it
## wants a named beta vector. methylQC 3.0.0 handed it a SigDF, so every call
## errored inside a tryCatch and every sample came back NA.
##
## WHY NOT THE v2 RULE AS WRITTEN: v2 chose the chrY cut by minimising the sum
## of within-half |residuals| from lm(chrY ~ chrX). That quantity is minimised
## by cutting anywhere through a homogeneous cloud, so the search could never
## discover that there was nothing to split. It required only >= 3 samples on
## each side, which any continuous distribution satisfies. Simulated over 200
## replicates per cohort shape, the consequence was:
##
##     cohort              spurious "sex mismatch" calls
##     20F +  0M           50.0% of the cohort, in 100% of runs
##      0F + 20M           49.8% of the cohort, in 100% of runs
##     38F +  2M           13.2% of the cohort, in 100% of runs
##     50F + 50M            0.2%
##
## A single-sex cohort is entirely ordinary -- a female breast series, a male
## prostate series -- so the rule failed hardest exactly where a sample swap
## matters most. v3.0.1 replaces the cut criterion with a SEPARATION criterion
## and refuses to call when the separation is not there.
## ---------------------------------------------------------------------------

## Minimum samples per cluster before a split is considered real.
.MQC_SEX_MINCLUST <- 3L

## Minimum ratio between the two clusters' median chrY intensity. Males carry
## roughly three to four times the female chrY signal, so a genuine split has a
## ratio well above this. A chance split of one cloud does not: the widest gap
## between order statistics of a normal sample can reach one within-group SD
## purely at random when n is small, which is enough to satisfy a gap test on
## its own, but it moves the two medians apart by only a few per cent.
.MQC_SEX_MINRATIO <- 1.5

## Smallest female cluster in which the high-confidence trim is applied.
.MQC_SEX_TRIMMIN <- 10L

#' Median sex-chromosome intensity for one sample
#'
#' Total (M+U, both channels) intensity summarised over the curated chrY and
#' X-linked probe sets. Called inside the worker while the \code{SigDF} is
#' still in hand, so the sex check does not depend on \code{savesdf}.
#'
#' Probe identifiers are stripped of EPICv2 replicate suffixes before matching.
#' The curated lists are bare \code{cg} identifiers, so on EPICv2 -- where the
#' \code{SigDF} carries \code{cg00004963_TC21} style identifiers and this runs
#' before any collapse -- an unstripped match returns nothing and every sample
#' reports \code{NA}.
#'
#' @param sdf a \code{SigDF}.
#' @return a list with \code{chrX} and \code{chrY} median intensities, each
#'   \code{NA_real_} when fewer than 10 probes of that set are present.
#' @keywords internal
#' @noRd
sexintensity <- function(sdf) {
  na <- list(chrX = NA_real_, chrY = NA_real_)
  pid <- sdf$Probe_ID
  if (is.null(pid)) return(na)
  pid <- stripv2(pid)

  z <- function(v) if (is.null(v)) 0 else ifelse(is.na(v), 0, v)
  ti <- z(sdf$MG) + z(sdf$MR) + z(sdf$UG) + z(sdf$UR)

  xi <- which(pid %in% .sesame_chrX_xlinked)
  yi <- which(pid %in% .sesame_chrY_clean)
  list(
    chrX = if (length(xi) >= 10L) stats::median(ti[xi], na.rm = TRUE) else NA_real_,
    chrY = if (length(yi) >= 10L) stats::median(ti[yi], na.rm = TRUE) else NA_real_)
}

#' Call sex for a cohort from sex-chromosome intensity
#'
#' Splits the cohort on log2 chrY intensity at the point of greatest
#' separation, then refuses the split unless the two groups are separated by a
#' clear band of empty space. Samples that sit far from both cluster axes are
#' returned as unclear rather than being forced into the nearer group.
#'
#' The separation statistic is the empty gap between the two groups divided by
#' the larger within-group standard deviation, all on the log2 scale. Log2
#' first because intensity is multiplicative in scanner gain and DNA input, so
#' a cohort-relative rule on the raw scale is not comparable across batches. A
#' genuine male/female split scores about 10; splitting a single-sex cohort at
#' its own widest internal gap scores about 0.3, so the default floor of 1.0
#' separates the two cases by an order of magnitude in either direction.
#'
#' @param chrX,chrY per-sample median intensities.
#' @param ids sample identifiers.
#' @param sep separation floor; below this no sex is called.
#' @param band orthogonal distance, in within-cluster standard deviations,
#'   outside which a sample is reported unclear.
#' @param minn cohort size below which no sex is called.
#' @param cutsd trim width for the high-confidence female set, in MADs about
#'   the low chrY cluster's median. It selects which females define the
#'   boundary; it does not itself set the cut-off.
#' @return a data.frame with \code{sample_id}, \code{inferred_sex}
#'   (\code{"M"}/\code{"F"}/\code{NA}), \code{sex_unclear} and
#'   \code{sex_confidence}. Attributes \code{called}, \code{reason},
#'   \code{threshold} and \code{separation} record what the cohort supported.
#' @export
#' @examples
#' set.seed(1)
#' x <- c(rnorm(10, 6000, 300), rnorm(10, 3200, 300))
#' y <- c(rnorm(10, 800, 80), rnorm(10, 3000, 250))
#' table(sexcall(x, y, paste0("s", 1:20))$inferred_sex)
sexcall <- function(chrX, chrY, ids, sep = NULL, band = NULL, minn = NULL,
                    cutsd = NULL) {
  cfg <- mqcopts()
  sep <- .opt(sep, "sexsep", cfg)
  band <- .opt(band, "sexband", cfg)
  minn <- .opt(minn, "sexmin", cfg)
  cutsd <- .opt(cutsd, "sexcutsd", cfg)

  n <- length(ids)
  out <- data.frame(sample_id = as.character(ids),
                    inferred_sex = NA_character_,
                    sex_unclear = TRUE,
                    sex_confidence = NA_real_,
                    stringsAsFactors = FALSE)
  refuse <- function(reason) {
    attr(out, "called") <- FALSE
    attr(out, "reason") <- reason
    attr(out, "threshold") <- NA_real_
    attr(out, "break") <- NA_real_
    attr(out, "separation") <- NA_real_
    out
  }

  ok <- is.finite(chrX) & is.finite(chrY) & chrX > 0 & chrY > 0
  if (n < minn) return(refuse(sprintf("cohort of %d is below sexmin = %d", n, minn)))
  if (sum(ok) < max(minn, 2L * .MQC_SEX_MINCLUST))
    return(refuse(sprintf("only %d sample(s) have usable sex-chromosome intensity",
                          sum(ok))))

  lx <- log2(chrX[ok]); ly <- log2(chrY[ok])

  ## ---- find the cut of greatest separation -------------------------------
  ## Candidates are the midpoints of consecutive sorted chrY values that leave
  ## at least .MQC_SEX_MINCLUST samples on each side. The score is the empty
  ## gap between the groups over the larger within-group SD, which is small for
  ## any cut through a single cloud and large only when a real gap exists.
  o <- order(ly)
  sorted <- ly[o]
  m <- length(sorted)
  lo_i <- .MQC_SEX_MINCLUST
  hi_i <- m - .MQC_SEX_MINCLUST
  if (hi_i < lo_i) return(refuse("too few samples to form two clusters"))

  best <- list(score = -Inf, cut = NA_real_, ratio = NA_real_)
  for (i in seq.int(lo_i, hi_i)) {
    a <- sorted[seq_len(i)]
    b <- sorted[seq.int(i + 1L, m)]
    sd_a <- stats::sd(a); sd_b <- stats::sd(b)
    spread <- max(sd_a, sd_b, na.rm = TRUE)
    if (!is.finite(spread) || spread <= 0) next
    score <- (min(b) - max(a)) / spread
    if (score > best$score)
      best <- list(score = score, cut = (max(a) + min(b)) / 2,
                   ratio = 2^(stats::median(b) - stats::median(a)))
  }
  if (!is.finite(best$score))
    return(refuse("sex-chromosome intensity has no usable spread"))

  ## Two independent conditions, because either alone is foolable. The gap test
  ## alone accepts a chance tail split in a small single-sex cohort; the ratio
  ## test alone accepts a smooth intensity gradient with no real grouping.
  if (best$score < sep)
    return(refuse(sprintf(
      paste("chrY intensity is unimodal (separation %.2f < %.2f); the cohort",
            "looks single-sex or the signal is too noisy, so no sex was called"),
      best$score, sep)))
  if (!is.finite(best$ratio) || best$ratio < .MQC_SEX_MINRATIO)
    return(refuse(sprintf(
      paste("the two chrY groups differ by only %.2fx (a real male/female split",
            "is 3-4x), so the cohort looks single-sex and no sex was called"),
      best$ratio)))

  ## ---- hard cut-off at the natural break ---------------------------------
  ## The break locates the two clusters. The low one is female: XX carries no Y
  ## chromosome, so its chrY readings are background. Measured on reference
  ## arrays the low cluster sits near 800 intensity units and the high near
  ## 3000.
  ##
  ## The cut-off is then placed HALFWAY between the highest high-confidence
  ## female and the next reading above her, which is taken to be the first
  ## male. Both terms come from this dataset, so the cut-off is dataset
  ## specific and needs no distributional assumption. A cut-off derived from
  ## the female cluster's MAD instead does assume normality, and MAD is 0.674
  ## sigma for a normal, so a k-MAD rule is really a 0.674k-sigma rule and a
  ## genuine female crosses it from time to time.
  ##
  ## "High confidence" trims stragglers before taking the maximum: a female
  ## with elevated chrY -- contamination, a mixed sample, XXY -- would
  ## otherwise drag the cut-off up behind her and let true males through under
  ## it. The trim is deliberately generous, so it removes only samples that sit
  ## clearly apart from their own cluster.
  low <- ly[ly <= best$cut]
  centre <- stats::median(low)
  spread_low <- stats::mad(low)
  if (!is.finite(spread_low) || spread_low <= 0) spread_low <- stats::sd(low)
  if (!is.finite(spread_low) || spread_low <= 0)
    return(refuse("the low chrY cluster has no measurable spread"))

  ## Trimming needs a cluster big enough for its own MAD to be stable. Below
  ## that, trimming removes a legitimate high female often enough to cost more
  ## than the contamination it catches.
  hiconf <- if (length(low) >= .MQC_SEX_TRIMMIN)
    low[low <= centre + cutsd * spread_low] else low
  if (length(hiconf) < .MQC_SEX_MINCLUST) hiconf <- low
  ntrim <- length(low) - length(hiconf)

  ## The first reading strictly above the high-confidence females, whether that
  ## is a trimmed straggler or the first male.
  above <- sorted[sorted > max(hiconf)]
  if (!length(above)) return(refuse("no reading above the female cluster"))
  cut <- (max(hiconf) + min(above)) / 2

  ## Hard cut-off: every sample with usable intensity gets a call, with no
  ## ambiguous band. Distance from the line, in female-cluster MADs, is
  ## reported so a borderline call is still visible in the sample sheet.
  call <- ifelse(ly > cut, "M", "F")
  conf <- abs(ly - cut) / spread_low

  out$inferred_sex[ok] <- call
  out$sex_unclear[ok] <- FALSE
  out$sex_confidence[ok] <- conf
  attr(out, "called") <- TRUE
  attr(out, "reason") <- "cohort is bimodal in chrY intensity"
  attr(out, "threshold") <- 2^cut
  attr(out, "break") <- 2^best$cut
  attr(out, "separation") <- best$score
  attr(out, "n_trimmed") <- ntrim
  attr(out, "female_median") <- 2^centre
  out
}
