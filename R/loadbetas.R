## ---------------------------------------------------------------------------
## loadbetas.R -- the downstream front door.
##
## Goal:    hand an analyst a masked beta matrix from a path, without their
##          having to know the output layout.
## Approach: find the Stage 1 matrices under whatever directory is given, read
##          the three files that belong together, and apply the masking policy
##          the caller names.
## Inputs:  an output directory, or a parent holding several.
## Outputs: a beta matrix, or a named list of them.
## Usage:   b <- loadbetas("~/qcout", maskuse = "both")
##
## cleanmat() is the in-memory operation and stays available for anyone who has
## already read the matrices themselves. This wraps it so the common case --
## "give me the betas from this run, masked" -- is one call rather than four
## readRDS() invocations against paths the caller has to remember.
## ---------------------------------------------------------------------------

#' Load a masked beta matrix from a methylQC output directory
#'
#' Finds the Stage 1 matrices beneath \code{dir}, reads the beta values, the
#' detection p-values and the design mask together, and applies the masking
#' policy given. Nothing on disk is altered.
#'
#' \code{dir} may be a single output directory or a parent holding several. In
#' the second case every run beneath it is loaded and a named list is returned,
#' so a set of projects can be pulled into one analysis without naming each
#' path.
#'
#' @param dir an output directory, or a parent containing several.
#' @param maskuse which masks to apply: \code{"both"}, \code{"detection"},
#'   \code{"design"} or \code{"none"}. Defaults to the \code{maskuse} option.
#' @param pthresh detection p-value threshold; defaults to the \code{detp}
#'   option.
#' @param mvals return M-values instead of beta values.
#' @param dropflagged drop the samples the QC report flagged. Requires the
#'   sample sheet, so it is skipped with a warning if Stage 2 has not run.
#' @param logger optional logger.
#' @return a probes-by-samples matrix, or a named list of them when \code{dir}
#'   holds more than one run.
#' @export
#' @examples
#' \dontrun{
#' b <- loadbetas("~/qcout")                              # design + detection
#' b <- loadbetas("~/qcout", maskuse = "detection")
#' m <- loadbetas("~/qcout", mvals = TRUE, dropflagged = TRUE)
#' all <- loadbetas("~/projects")                         # every run beneath
#' }
loadbetas <- function(dir, maskuse = NULL, pthresh = NULL, mvals = FALSE,
                      dropflagged = FALSE, logger = NULL) {
  logger <- logger %||% nulllog()
  dirs <- .find_outdirs(dir, key = "betas")
  many <- length(dirs) > 1L
  if (many)
    message(sprintf("methylQC: %d run(s) found under %s", length(dirs), dir))

  out <- lapply(seq_along(dirs), function(k) {
    d <- dirs[k]
    if (many) message(sprintf("  [%d/%d] %s", k, length(dirs), d))
    .loadbetas_one(d, maskuse, pthresh, mvals, dropflagged, logger)
  })
  ## Name by the path relative to what the caller asked for: three projects
  ## whose output directories are all called "out" would otherwise come back
  ## with three identical names.
  root <- sub(paste0(.Platform$file.sep, "?$"), "", normalizePath(dir, mustWork = FALSE))
  nm <- sub(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", root),
                   .Platform$file.sep, "?"), "",
            normalizePath(dirs, mustWork = FALSE))
  names(out) <- ifelse(nzchar(nm), nm, basename(dirs))
  if (!many) return(out[[1]])
  out
}

.loadbetas_one <- function(dir, maskuse, pthresh, mvals, dropflagged, logger) {
  betas <- readRDS(mqcpath(dir, "betas"))

  ## The design mask and the detection matrix are optional only in the sense
  ## that a policy may not need them; asking for a mask whose input is absent
  ## is an error rather than a silent no-op.
  detp <- if (file.exists(mqcpath(dir, "detp"))) readRDS(mqcpath(dir, "detp")) else NULL
  dm <- if (file.exists(mqcpath(dir, "designmask")))
    readRDS(mqcpath(dir, "designmask")) else NULL

  use <- .opt(maskuse, "maskuse")
  if (use %in% c("both", "detection") && is.null(detp))
    stop("maskuse = '", use, "' needs data/matrices/detP_all.rds, which is ",
         "absent from ", dir, call. = FALSE)
  if (use %in% c("both", "design") && is.null(dm))
    stop("maskuse = '", use, "' needs data/matrices/design_mask.rds, which is ",
         "absent from ", dir, call. = FALSE)

  b <- cleanmat(betas, detp = detp, design = dm, maskuse = use,
                pthresh = pthresh, logger = logger)

  if (isTRUE(dropflagged)) {
    p <- mqcpath(dir, "sheet")
    ss <- if (file.exists(p))
      tryCatch(utils::read.csv(p, stringsAsFactors = FALSE),
               error = function(e) NULL) else NULL
    if (is.null(ss) || !all(c("sample_id", "flagged") %in% names(ss))) {
      warning("dropflagged = TRUE but ", dir, " has no sample sheet with a ",
              "'flagged' column; run prep() first. Keeping every sample.",
              call. = FALSE)
    } else {
      drop <- as.character(ss$sample_id[isTRUE_vec(ss$flagged)])
      keep <- setdiff(colnames(b), drop)
      logger$log("load", sprintf("dropping %d flagged sample(s), keeping %d",
                                 length(colnames(b)) - length(keep), length(keep)))
      b <- b[, keep, drop = FALSE]
    }
  }

  if (isTRUE(mvals)) b <- mvals(b)
  b
}
