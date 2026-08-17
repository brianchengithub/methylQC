## ---------------------------------------------------------------------------
## redetect.R -- rebuild the detection matrix without re-reading IDATs.
##
## Goal:    let a run made before the ELBAR ordering fix be repaired in minutes
##          rather than reprocessed from raw data.
## Approach: the retained SigDFs are post-noob, which is exactly what ELBAR
##          needs, so the detection p-values can simply be recomputed from them.
## Inputs:  an output directory holding sdfs_all.rds.
## Outputs: a rewritten detP_all.rds, and an invalidated QC cache.
## Usage:   redetect("~/qcout"); prep("~/qcout")
##
## Why this is sound. The stored SigDFs have had C, D and B applied before
## being retained -- that has been true of every version that wrote them -- and
## ELBAR(return.pval = TRUE) reads a SigDF and returns p-values without
## modifying it. Re-running it on the stored object therefore produces exactly
## the p-values a fresh run would, and the beta matrix is not involved at all.
##
## What this does NOT repair: anything derived from the raw IDATs that is not
## the detection matrix. Per-sample instrument statistics, the design mask and
## the betas are all unchanged by the ordering fix, so they do not need
## repairing; but a run predating the sex-chromosome intensity columns or the
## chip-position columns still needs a full pipeline() to acquire those.
## ---------------------------------------------------------------------------

#' Recompute detection p-values from retained SigDF objects
#'
#' Rebuilds \code{detP_all.rds} from \code{sdfs_all.rds}, then invalidates the
#' QC cache so that \code{\link{prep}} recomputes everything downstream. The
#' beta matrix is untouched, because it does not depend on when detection ran.
#'
#' The intended use is repairing a run made before methylQC 3.0.2, where ELBAR
#' was called before rather than after \code{noob} and the resulting p-values
#' were unusable (see \code{METHODS.md} section 2.2). It needs the run to have
#' been made with \code{savesdf = TRUE}, which is the default; without the
#' SigDFs there is nothing to recompute from and the run has to be repeated.
#'
#' @param dir an output directory, or a parent containing several.
#' @param logger optional logger.
#' @return the path(s) written, invisibly.
#' @export
#' @examples
#' \dontrun{
#' redetect("~/qcout")   # rewrite detP_all.rds from the stored SigDFs
#' prep("~/qcout")       # then redo Stage 2 against it
#' }
redetect <- function(dir, logger = NULL) {
  dirs <- .find_outdirs(dir, key = "betas")
  many <- length(dirs) > 1L
  out <- character(0)
  for (k in seq_along(dirs)) {
    if (many) message(sprintf("  [%d/%d] %s", k, length(dirs), dirs[k]))
    out <- c(out, .redetect_one(dirs[k], logger))
  }
  invisible(if (length(out) == 1L) out[[1]] else out)
}

.redetect_one <- function(dir, logger = NULL) {
  own <- is.null(logger)
  lg <- logger %||% makelog(mqcpath(dir, "log", create = TRUE))
  if (own) on.exit(lg$close(), add = TRUE)

  sp <- mqcpath(dir, "sdfs")
  if (!file.exists(sp))
    stop("no ", basename(sp), " in ", dir, ". Detection p-values can only be ",
         "recomputed from retained SigDF objects; without them the run has to ",
         "be repeated with pipeline().", call. = FALSE)

  betas <- readRDS(mqcpath(dir, "betas"))
  sdfs <- readRDS(sp)
  if (is.null(names(sdfs))) names(sdfs) <- colnames(betas)
  miss <- setdiff(colnames(betas), names(sdfs))
  if (length(miss))
    stop(length(miss), " sample(s) in the beta matrix have no retained SigDF ",
         "(e.g. '", miss[1], "'); rerun pipeline() instead.", call. = FALSE)

  lg$log("redetect", sprintf("recomputing detection for %d sample(s) from %s",
                             ncol(betas), basename(sp)))

  ids <- rownames(betas)
  detp <- matrix(NA_real_, length(ids), ncol(betas),
                 dimnames = list(ids, colnames(betas)))
  warned <- 0L
  for (j in seq_len(ncol(betas))) {
    s <- sdfs[[colnames(betas)[j]]]
    if (is.null(s)) next                       # failed sample: stays all-NA
    w <- NULL
    pv <- withCallingHandlers(
      tryCatch(sesame::ELBAR(s, return.pval = TRUE), error = function(e) NULL),
      warning = function(x) { w <<- c(w, conditionMessage(x))
                              invokeRestart("muffleWarning") })
    if (any(grepl("dichotomous", w %||% ""))) warned <- warned + 1L
    if (is.null(pv)) next
    pv[is.na(pv)] <- 1                          # fail closed, as Stage 1 does
    detp[, j] <- align_to(pv, ids, paste0("detection p-values for ",
                                          colnames(betas)[j]), fill = 1)
  }
  ## A failed sample has no p-values; it must look maximally suspect.
  detp[is.na(detp)] <- 1

  ## The dichotomy warning here would mean the stored SigDFs are pre-noob,
  ## which no released version produced -- worth saying loudly rather than
  ## writing a matrix that is wrong in the same way as the one being replaced.
  if (warned > 0L)
    lg$log("redetect", sprintf(
      paste("ELBAR reported a dichotomous background for %d sample(s). The",
            "stored SigDFs may not be noob-corrected; check the result before",
            "trusting it."), warned), warn = TRUE)

  save_rds_atomic(detp, mqcpath(dir, "detp"))
  lg$log("redetect", sprintf("wrote %s", mqcpath(dir, "detp")))

  ## The cache holds call rates and probe-failure counts derived from the old
  ## matrix; leaving it in place would let qcplots() serve stale numbers.
  cp <- mqcpath(dir, "cache")
  if (file.exists(cp)) {
    unlink(cp)
    lg$log("redetect", "removed the stale QC cache; run prep() to rebuild it")
  }
  mqcpath(dir, "detp")
}
