## ---------------------------------------------------------------------------
## celltypes.R -- cell composition deconvolution and SNP identity checks.
##
## Goal:    estimate cell type proportions, and detect sample swaps from the
##          rs genotyping probes.
## Approach: EpiDISH for deconvolution, on the samples that have usable data;
##          a vectorised pairwise genotype concordance for identity.
## Inputs:  the beta matrix.
## Outputs: a per-sample proportion frame, and a filtered pair list.
## Usage:   props <- rundish(betas); pairs <- snpcheck(snpbetas(betas))
## ---------------------------------------------------------------------------

## The EpiDISH reference panels that actually exist in the installed package.
## 3.0.0 offered a "breast" panel mapping to centBreastSub.m, which is not an
## EpiDISH dataset; the lookup threw, the error was swallowed, and cell
## composition was silently skipped for anyone who selected it.
.MQC_DISH_REFS <- c(
  blood      = "centDHSbloodDMC.m",
  bloodsub   = "centBloodSub.m",
  epithelial = "centEpiFibIC.m",
  epifibfat  = "centEpiFibFatIC.m",
  `12ct`     = "cent12CT.m")

#' Estimate cell type proportions
#'
#' Wraps \pkg{EpiDISH}. Reference probe sets are keyed on bare \code{cg}
#' identifiers, so EPICv2 data must have its replicate suffixes resolved
#' first (see \code{\link{collapsev2}}); otherwise no reference probe matches.
#'
#' @param betas probes-by-samples beta matrix.
#' @param ref reference panel; see \code{mqcdefaults()} for the valid names.
#' @param method deconvolution method passed to EpiDISH.
#' @param logger optional logger.
#' @return a data.frame with \code{sample_id} and one column per cell type, or
#'   \code{NULL} if EpiDISH is unavailable or nothing matched.
#' @export
#' @examples
#' \dontrun{
#' props <- rundish(betas, ref = "blood")
#' }
rundish <- function(betas, ref = NULL, method = NULL, logger = NULL) {
  logger <- logger %||% nulllog()
  cfg <- mqcopts()
  ref <- .opt(ref, "dishref", cfg)
  method <- .opt(method, "dishmethod", cfg)

  if (!have_pkg("EpiDISH", "cell composition estimation", logger)) return(NULL)

  refset <- tryCatch(.dish_ref(ref), error = function(e) NULL)
  if (is.null(refset)) {
    logger$log("cells", sprintf("unknown reference panel '%s'; valid panels are %s",
                                ref, paste(names(.MQC_DISH_REFS), collapse = ", ")),
               warn = TRUE)
    return(NULL)
  }

  ## EPICv2 identifiers carry replicate suffixes; match on stripped names when
  ## the matrix has not been collapsed.
  rn <- rownames(betas)
  if (has_v2_suffix(rn)) rn <- stripv2(rn)
  hit <- match(rownames(refset), rn)
  keep <- which(!is.na(hit))
  logger$log("cells", sprintf("%d reference probe(s) of %d matched the array",
                              length(keep), nrow(refset)))
  if (length(keep) < 20L) {
    logger$log("cells", sprintf(
      paste("only %d reference probe(s) matched; skipping deconvolution.",
            "On EPICv2 this usually means replicate suffixes were not",
            "resolved -- see collapsev2()."), length(keep)), warn = TRUE)
    return(NULL)
  }

  bm <- betas[hit[keep], , drop = FALSE]
  rownames(bm) <- rownames(refset)[keep]
  rm_ <- refset[keep, , drop = FALSE]

  ## A single failed IDAT is an all-NA column, and EpiDISH fails on the whole
  ## matrix if any column is unusable -- so in 3.0.0 one bad array destroyed
  ## cell composition for the entire cohort. Deconvolve the usable samples and
  ## return NA rows for the rest.
  usable <- colSums(!is.na(bm)) >= 20L
  if (!any(usable)) {
    logger$log("cells", "no sample has enough non-missing reference probes",
               warn = TRUE)
    return(NULL)
  }
  if (any(!usable))
    logger$log("cells", sprintf(
      "%d sample(s) lack usable reference probes and get NA proportions",
      sum(!usable)), warn = TRUE)

  sub <- bm[, usable, drop = FALSE]
  ok <- stats::complete.cases(sub)
  if (sum(ok) < 20L) {
    logger$log("cells", "too few complete reference probes across samples",
               warn = TRUE)
    return(NULL)
  }

  res <- tryCatch(
    EpiDISH::epidish(beta.m = sub[ok, , drop = FALSE],
                     ref.m = rm_[ok, , drop = FALSE], method = method),
    error = function(e) {
      logger$log("cells", paste("EpiDISH failed:", conditionMessage(e)),
                 warn = TRUE)
      NULL
    })
  if (is.null(res) || is.null(res$estF)) return(NULL)

  props <- as.data.frame(res$estF, stringsAsFactors = FALSE)
  names(props) <- paste0("cell_", names(props))
  out <- data.frame(sample_id = colnames(betas), stringsAsFactors = FALSE)
  m <- match(out$sample_id, rownames(res$estF))
  for (cl in names(props)) out[[cl]] <- props[[cl]][m]

  logger$log("cells", sprintf(
    "estimated %d cell type(s) for %d of %d sample(s) using %s over %d probes",
    ncol(props), sum(!is.na(m)), nrow(out), method, sum(ok)))
  attr(out, "n_probes") <- sum(ok)
  attr(out, "n_types") <- ncol(props)
  out
}

.dish_ref <- function(ref) {
  nm <- .MQC_DISH_REFS[[as.character(ref)]]
  if (is.null(nm)) stop("unknown reference '", ref, "'", call. = FALSE)
  e <- new.env(parent = emptyenv())
  utils::data(list = nm, package = "EpiDISH", envir = e)
  get(nm, envir = e)
}


## ---------------------------------------------------------------------------
## SNP identity
## ---------------------------------------------------------------------------

#' Extract SNP (rs) probe beta values
#'
#' Infinium arrays carry a small panel of rs probes that genotype the sample.
#' They are useful for detecting sample swaps and for confirming that repeated
#' measurements come from the same donor.
#'
#' @param betas probes-by-samples beta matrix.
#' @return the rs-probe submatrix, or \code{NULL} if none are present.
#' @export
#' @examples
#' b <- matrix(runif(4), 2, 2, dimnames = list(c("rs001", "cg1"), c("a", "b")))
#' dim(snpbetas(b))
snpbetas <- function(betas) {
  idx <- grep("^rs", rownames(betas))
  if (!length(idx)) return(NULL)
  betas[idx, , drop = FALSE]
}

#' Check donor identity using SNP probes
#'
#' Genotypes are called by thresholding the rs beta values into three classes
#' and compared pairwise. Pairs at or above \code{minmatch} are reported as the
#' same individual.
#'
#' The comparison is quadratic in samples, so at cohort scale the full pair
#' list is not returned: only pairs that are interesting are kept, meaning
#' those at or above \code{keepmin} concordance plus every pair the donor
#' labels say should match. A swap shows up in the second group -- a pair
#' expected to match that does not -- so filtering on concordance alone would
#' hide exactly the case worth finding.
#'
#' @param snps rs-probe beta matrix from \code{\link{snpbetas}}.
#' @param donors optional named vector mapping sample IDs to donor IDs.
#' @param minmatch concordance above which two samples are called identical.
#' @param keepmin concordance below which an unrelated pair is not reported.
#' @param logger optional logger.
#' @return a data.frame of sample pairs with their concordance, or \code{NULL}.
#' @export
#' @examples
#' \dontrun{
#' pairs <- snpcheck(snpbetas(betas))
#' }
snpcheck <- function(snps, donors = NULL, minmatch = 0.90, keepmin = NULL,
                     logger = NULL) {
  logger <- logger %||% nulllog()
  cfg <- mqcopts()
  keepmin <- .opt(keepmin, "snpmin", cfg)

  if (is.null(snps) || nrow(snps) < 10L || ncol(snps) < 2L) {
    logger$log("snp", "too few rs probes or samples for an identity check")
    return(NULL)
  }

  ## 0 / 1 / 2 genotype classes from beta values
  g <- matrix(NA_integer_, nrow(snps), ncol(snps), dimnames = dimnames(snps))
  g[snps < 0.25] <- 0L
  g[snps >= 0.25 & snps <= 0.75] <- 1L
  g[snps > 0.75] <- 2L
  n <- ncol(g)
  ids <- colnames(g)

  ## Concordance for every pair at once. 3.0.0 looped over combn(n, 2) with a
  ## closure per pair, which at 1500 samples is 1,124,250 closure calls and a
  ## 1.1-million-row CSV. Three crossproducts give the same numbers in one
  ## pass and let the filtering happen before any frame is built.
  ind <- function(k) { m <- (g == k); m[is.na(m)] <- FALSE; m * 1 }
  obs <- (!is.na(g)) * 1
  match_n <- crossprod(ind(0L)) + crossprod(ind(1L)) + crossprod(ind(2L))
  valid_n <- crossprod(obs)
  conc <- ifelse(valid_n > 0, match_n / valid_n, NA_real_)

  ut <- which(upper.tri(conc), arr.ind = TRUE)
  a <- ids[ut[, 1]]; b <- ids[ut[, 2]]
  cv <- conc[ut]

  exp_same <- rep(NA, length(cv))
  if (!is.null(donors)) {
    da <- unname(donors[match(a, names(donors))])
    db <- unname(donors[match(b, names(donors))])
    exp_same <- !is.na(da) & !is.na(db) & da == db
  }

  keep <- (!is.na(cv) & cv >= keepmin) | (!is.na(exp_same) & exp_same)
  n_all <- length(cv)
  if (!any(keep)) {
    logger$log("snp", sprintf(
      "checked %d pair(s) over %d rs probes; none reached %.2f concordance",
      n_all, nrow(g), keepmin))
    return(NULL)
  }

  out <- data.frame(
    sample_a = a[keep], sample_b = b[keep], concordance = cv[keep],
    same_individual = !is.na(cv[keep]) & cv[keep] >= minmatch,
    stringsAsFactors = FALSE)
  if (!is.null(donors)) {
    out$expected_same <- exp_same[keep]
    out$discordant <- !is.na(out$expected_same) &
      (out$expected_same != out$same_individual)
    nd <- sum(out$discordant, na.rm = TRUE)
    if (nd > 0)
      logger$log("snp", sprintf(
        "%d sample pair(s) disagree with the donor labels; possible sample swap",
        nd), warn = TRUE)
  }

  out <- out[order(-out$concordance), , drop = FALSE]
  rownames(out) <- NULL
  logger$log("snp", sprintf(
    "checked %d pair(s) over %d rs probes; %d matched at >= %.2f, %d row(s) reported",
    n_all, nrow(g), sum(out$same_individual, na.rm = TRUE), minmatch, nrow(out)))
  out
}
