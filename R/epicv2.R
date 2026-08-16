## ---------------------------------------------------------------------------
## epicv2.R -- EPICv2 replicate probe resolution.
##
## EPICv2 measures some CpG sites with more than one physical probe design,
## distinguished by a suffix: cg00004963_TC21 and cg00004963_TC22 are two
## different probes measuring the same CpG. Older arrays have no suffix.
##
## Different probe designs at the same site give systematically different beta
## values because probe design affects hybridisation, so averaging replicates
## produces a value matching neither EPICv1 nor whole-genome bisulphite
## sequencing. Choosing the replicate with the best detection p-value picks
## the brightest probe, which answers a different question from which probe is
## most accurate.
##
## Peters et al. (2024, BMC Genomics 25:251) evaluated every replicate group
## against matched EPICv1 and WGBS data and published recommendations as to
## the superior probe. That recommendation is static, empirical and validated
## against an independent sequencing ground truth, so methylQC uses it in
## preference to either alternative.
##
## The table is frozen into inst/extdata at build time (see
## data-raw/build_epicv2_replicates.R) rather than read from the annotation
## package at run time, so that an annotation update cannot silently change
## results for an existing analysis.
##
## IMPORTANT: methylQC does not ship the table itself, because it is derived
## from the EPICv2manifest annotation package. Run build_epicv2_table() once
## to generate it. Until then collapsemethod = "peters" falls back to
## "minpval" with a warning.
## ---------------------------------------------------------------------------

#' Path to the frozen EPICv2 replicate table
#' @keywords internal
#' @noRd
.epicv2_table_path <- function() {
  system.file("extdata", "epicv2_replicates.rds", package = "methylQC")
}

#' Load the frozen EPICv2 replicate selection table
#'
#' @return a data.frame with columns \code{group} (stripped probe ID),
#'   \code{chosen} (full EPICv2 IlmnID to keep) and \code{criterion}; or
#'   \code{NULL} if the table has not been built.
#' @export
#' @examples
#' str(epicv2table())
epicv2table <- function() {
  p <- .epicv2_table_path()
  if (!nzchar(p) || !file.exists(p)) return(NULL)
  tab <- tryCatch(readRDS(p), error = function(e) NULL)
  if (is.null(tab) || !all(c("group", "chosen") %in% names(tab))) return(NULL)
  tab
}

#' Build the frozen EPICv2 replicate selection table
#'
#' Derives, for each EPICv2 replicate group, the probe recommended by Peters
#' et al. (2024) and writes the result to \code{inst/extdata} so that it ships
#' with the installed package. Requires the \pkg{EPICv2manifest} annotation
#' package.
#'
#' Run this once, from the package source directory, then rebuild the package.
#'
#' @param out destination path for the frozen table.
#' @param criterion which Peters column to prefer. \code{"auto"} takes the
#'   first available of the sensitivity/precision recommendation columns.
#' @return the table, invisibly.
#' @export
#' @examples
#' \dontrun{
#' build_epicv2_table("inst/extdata/epicv2_replicates.rds")
#' }
build_epicv2_table <- function(out = "inst/extdata/epicv2_replicates.rds",
                               criterion = "auto") {
  if (!requireNamespace("EPICv2manifest", quietly = TRUE))
    stop("EPICv2manifest is required to build the table.\n",
         "  install.packages(\"pak\")\n",
         "  pak::pkg_install(\"bioc::EPICv2manifest\")", call. = FALSE)

  e <- new.env(parent = emptyenv())
  utils::data("EPICv2manifest", package = "EPICv2manifest", envir = e)
  mf <- get(ls(e)[1], envir = e)
  mf <- as.data.frame(mf, stringsAsFactors = FALSE)
  if (is.null(mf$IlmnID)) mf$IlmnID <- rownames(mf)

  if (is.null(mf$Name))
    stop("EPICv2manifest has no 'Name' column; cannot group replicates.",
         call. = FALSE)

  ## Peters' per-probe recommendation columns. The exact names have moved
  ## between annotation releases, so probe for them rather than hard-coding.
  cand <- grep("(?i)(sensitivity|precision|rep_?result|posrep|superior)",
               names(mf), value = TRUE, perl = TRUE)
  if (identical(criterion, "auto")) {
    crit_col <- if (length(cand)) cand[1] else NA_character_
  } else {
    if (!criterion %in% names(mf))
      stop("column '", criterion, "' not found. Available candidates: ",
           paste(cand, collapse = ", "), call. = FALSE)
    crit_col <- criterion
  }

  ## Restrict to probes that actually have replicates. `namerep` marks these.
  if (!is.null(mf$namerep)) {
    mf <- mf[!is.na(mf$namerep) & mf$namerep == "Y", , drop = FALSE]
  } else {
    dup <- mf$Name[duplicated(mf$Name)]
    mf <- mf[mf$Name %in% dup, , drop = FALSE]
  }
  if (!nrow(mf)) stop("no replicate groups found in the manifest", call. = FALSE)

  sp <- split(seq_len(nrow(mf)), mf$Name)
  pick <- vapply(sp, function(idx) {
    if (length(idx) == 1L) return(mf$IlmnID[idx])
    if (!is.na(crit_col)) {
      v <- mf[[crit_col]][idx]
      if (is.numeric(v) && any(!is.na(v))) return(mf$IlmnID[idx[which.max(v)]])
      hit <- which(!is.na(v) & toupper(as.character(v)) %in% c("Y", "TRUE", "1"))
      if (length(hit)) return(mf$IlmnID[idx[hit[1]]])
    }
    mf$IlmnID[idx[1]]          # manifest order, sesame's own fallback
  }, character(1))

  tab <- data.frame(
    group     = names(sp),
    chosen    = unname(pick),
    criterion = if (is.na(crit_col)) "manifest_order" else crit_col,
    built     = as.character(Sys.Date()),
    stringsAsFactors = FALSE)

  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(tab, out)
  message(sprintf("wrote %d replicate groups to %s (criterion: %s)",
                  nrow(tab), out, tab$criterion[1]))
  invisible(tab)
}

#' Strip the EPICv2 replicate suffix from probe IDs
#'
#' \code{cg00004963_TC21} becomes \code{cg00004963}. IDs without a suffix are
#' returned unchanged, so this is safe to apply to any platform.
#'
#' @param ids character vector of probe IDs
#' @return character vector with suffixes removed
#' @export
#' @examples
#' stripv2(c("cg00004963_TC21", "cg00000029"))
stripv2 <- function(ids) sub("_[A-Z]{2}[0-9]{2}$", "", ids)

#' Does this probe set carry EPICv2 replicate suffixes?
#'
#' @param ids probe IDs
#' @return \code{TRUE} if any ID carries a replicate suffix
#' @keywords internal
#' @noRd
has_v2_suffix <- function(ids) any(grepl("_[A-Z]{2}[0-9]{2}$", ids))

#' Collapse EPICv2 replicate probes
#'
#' Reduces suffixed EPICv2 probe IDs to one row per CpG. Betas, detection
#' p-values and the design mask are collapsed together in a single operation
#' keyed on the same group membership, so they cannot disagree afterwards.
#'
#' This is the fix for a silent defect in methylQC v2: there, \code{collapse =
#' TRUE} collapsed only the beta matrix, leaving the mask and detection lookups
#' keyed on suffixed IDs that no longer matched, so both came back entirely
#' \code{NA} and the run applied no quality control at all while producing
#' files of the correct shape.
#'
#' @param betas probes-by-samples numeric matrix with suffixed row names.
#' @param detp matching detection p-value matrix, or \code{NULL}.
#' @param design named logical design mask, or \code{NULL}.
#' @param method \code{"peters"} (keep the recommended superior probe),
#'   \code{"minpval"} (keep the replicate with the lowest median detection
#'   p-value) or \code{"mean"} (average the replicates).
#' @param logger optional logger.
#' @return a list with collapsed \code{betas}, \code{detp}, \code{design} and a
#'   \code{report} data.frame recording how many groups took each path.
#' @export
#' @examples
#' b <- matrix(runif(6), 3, 2,
#'             dimnames = list(c("cg1_TC21", "cg1_TC22", "cg2_BC11"), c("s1", "s2")))
#' collapsev2(b, method = "mean")$betas
collapsev2 <- function(betas, detp = NULL, design = NULL,
                       method = c("peters", "minpval", "mean"),
                       logger = NULL) {
  method <- match.arg(method)
  logger <- logger %||% nulllog()
  stopifnot(is.matrix(betas), !is.null(rownames(betas)))

  ids <- rownames(betas)
  grp <- stripv2(ids)

  if (!is.null(detp)) detp <- conform_to(detp, betas, "detP matrix")
  if (!is.null(design)) design <- align_to(design, ids, "design mask")

  if (identical(method, "peters")) {
    tab <- epicv2table()
    if (is.null(tab)) {
      logger$log("epicv2",
        paste("the frozen Peters replicate table is not installed;",
              "falling back to method = 'minpval'. Run build_epicv2_table()",
              "to enable 'peters'."), warn = TRUE)
      method <- "minpval"
    }
  }

  ## ---- decide, per group, which rows to keep ------------------------------
  keep_idx <- NULL
  path <- character(0)

  if (identical(method, "mean")) {
    ord <- order(grp)
    b <- rowsum(betas[ord, , drop = FALSE], group = grp[ord], reorder = TRUE)
    n <- as.vector(table(grp[ord]))
    b <- b / n
    out_betas <- b
    out_detp <- if (is.null(detp)) NULL else {
      d <- rowsum(detp[ord, , drop = FALSE], group = grp[ord], reorder = TRUE)
      d / n
    }
    out_design <- if (is.null(design)) NULL else {
      ## a collapsed probe is design-masked only if EVERY replicate is
      s <- rowsum(as.numeric(design)[ord], group = grp[ord], reorder = TRUE)
      stats::setNames(as.vector(s) == n, rownames(b))
    }
    report <- data.frame(path = "mean", n_groups = nrow(b),
                         stringsAsFactors = FALSE)
    logger$log("epicv2", sprintf("collapsed %d probes to %d CpGs by mean",
                                 length(ids), nrow(b)))
    return(list(betas = out_betas, detp = out_detp, design = out_design,
                report = report, method = "mean"))
  }

  ## Selection methods: pick exactly one row per group.
  sp <- split(seq_along(ids), grp)
  singles <- lengths(sp) == 1L

  chosen <- integer(length(sp))
  path <- rep(NA_character_, length(sp))
  chosen[singles] <- vapply(sp[singles], function(i) i[1], integer(1))
  path[singles] <- "single"

  multi <- which(!singles)
  if (length(multi)) {
    if (identical(method, "peters")) {
      tab <- epicv2table()
      lut <- stats::setNames(tab$chosen, tab$group)
      gn <- names(sp)[multi]
      want <- unname(lut[gn])
      for (k in seq_along(multi)) {
        i <- sp[[multi[k]]]
        hit <- if (!is.na(want[k])) match(want[k], ids[i]) else NA_integer_
        if (!is.na(hit)) {
          chosen[multi[k]] <- i[hit]; path[multi[k]] <- "peters"
        } else {
          chosen[multi[k]] <- .pick_minpval(i, detp)
          path[multi[k]] <- "minpval_fallback"
        }
      }
    } else {
      for (k in multi) {
        chosen[k] <- .pick_minpval(sp[[k]], detp)
        path[k] <- "minpval"
      }
    }
  }

  keep_idx <- chosen
  out_betas <- betas[keep_idx, , drop = FALSE]
  rownames(out_betas) <- names(sp)
  out_detp <- if (is.null(detp)) NULL else {
    d <- detp[keep_idx, , drop = FALSE]; rownames(d) <- names(sp); d
  }
  out_design <- if (is.null(design)) NULL else
    stats::setNames(unname(design[keep_idx]), names(sp))

  report <- as.data.frame(table(path), stringsAsFactors = FALSE)
  names(report) <- c("path", "n_groups")
  logger$log("epicv2",
             sprintf("collapsed %d probes to %d CpGs (%s)",
                     length(ids), length(sp),
                     paste(sprintf("%s=%d", report$path, report$n_groups),
                           collapse = ", ")))

  list(betas = out_betas, detp = out_detp, design = out_design,
       report = report, method = method)
}

## Lowest median detection p-value across samples wins. With no detP matrix
## available, fall back to manifest order (sesame's own default).
.pick_minpval <- function(i, detp) {
  if (is.null(detp)) return(i[1])
  med <- .rowMedians(detp[i, , drop = FALSE])
  if (all(is.na(med))) return(i[1])
  i[which.min(med)]
}
