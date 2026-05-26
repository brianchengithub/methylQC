###############################################################################
# epidish.R — EpiDISH cell-type deconvolution
#
# Single user-facing function: rundish()
#
# Dispatches on the first argument:
#   rundish(<matrix>, ...)        in-memory  -- returns the proportions
#   rundish("<dir>",   ...)       on-disk    -- runs in a Stage 2 output dir
#                                                writes cell_proportions.csv
#                                                merges into sample_sheet.csv
#                                                returns the proportions
#
# Reference resolution:
#   dishref can be a single string ("centDHSbloodDMC.m") or a named list
#   keyed by platform (list(EPIC = "cent12CT.m", "450k" = "cent12CT450k.m")).
#   When a list is given, the function picks the entry whose name matches
#   the `platform` argument, falling back to case-insensitive match and
#   then to contains-match ("EPICv2" -> "EPIC").
#
# Tissue gate:
#   There isn't one. rundish() runs whatever you point it at; the only
#   refusal condition is too few probes overlapping the reference (< 50),
#   which is a data-level check, not a label-level guess.
###############################################################################

#' Resolve dishref against a platform
#'
#' Accepts a single string (returned as-is, no platform needed) or a
#' named list keyed by platform. Resolution order: exact name match,
#' case-insensitive name match, then contains-match
#' (e.g. "EPICv2" matches a key named "EPIC").
#'
#' @param dishref The dishref config value (string or list).
#' @param platform Platform string. Required iff dishref is a list.
#' @return Single character reference name.
#' @keywords internal
#' @noRd
resolve_dishref <- function(dishref, platform = NULL) {
  if (is.character(dishref) && length(dishref) == 1) return(dishref)
  if (is.list(dishref)) {
    if (is.null(platform))
      stop("dishref is a platform-keyed list but no 'platform' was supplied.",
           call. = FALSE)
    if (platform %in% names(dishref)) return(as.character(dishref[[platform]]))
    idx <- which(tolower(names(dishref)) == tolower(platform))
    if (length(idx) == 1) return(as.character(dishref[[idx]]))
    idx <- which(vapply(names(dishref),
                        function(k) grepl(k, platform, ignore.case = TRUE),
                        logical(1)))
    if (length(idx) >= 1) return(as.character(dishref[[idx[1]]]))
    stop(sprintf("No dishref entry for platform '%s'. Available: %s.",
                 platform, paste(names(dishref), collapse = ", ")),
         call. = FALSE)
  }
  stop("dishref must be a single string or a named list keyed by platform.",
       call. = FALSE)
}

#' Load EpiDISH reference into a local environment
#' @keywords internal
#' @noRd
load_dishref <- function(ref) {
  ref_env <- new.env()
  utils::data(list = ref, package = "EpiDISH", envir = ref_env)
  if (!exists(ref, envir = ref_env))
    stop(sprintf("EpiDISH reference '%s' not found. Update EpiDISH: ",
                 ref),
         "BiocManager::install(\"EpiDISH\")", call. = FALSE)
  get(ref, envir = ref_env)
}

#' Run EpiDISH on a beta matrix (internal in-memory primitive)
#'
#' @keywords internal
#' @noRd
rundish_matrix <- function(betas, platform = NULL, ref = NULL,
                           method = NULL, logger = NULL) {
  cfg <- mqcopts()
  if (is.null(ref))    ref    <- resolve_dishref(cfg$dishref, platform)
  if (is.null(method)) method <- cfg$dishmethod

  ref_data <- load_dishref(ref)
  overlap <- intersect(rownames(betas), rownames(ref_data))
  if (length(overlap) < 50)
    stop(sprintf(
      "Only %d probes overlap with reference '%s'. Check that 'betas' contains the reference's CpGs (typical: ~600 for cent12CT.m, ~333 for centDHSbloodDMC.m).",
      length(overlap), ref), call. = FALSE)

  if (!is.null(logger))
    logger$log("rundish",
               sprintf("reference '%s' (method=%s): %d / %d probes present (%.1f%%)",
                       ref, method, length(overlap), nrow(ref_data),
                       100 * length(overlap) / nrow(ref_data)))

  res <- EpiDISH::epidish(beta.m = betas[overlap, , drop = FALSE],
                          ref.m  = ref_data[overlap, , drop = FALSE],
                          method = method)$estF

  if (!is.null(logger)) {
    rs <- rowSums(res)
    logger$log("rundish",
               sprintf("row-sum check: median=%.3f range=[%.3f, %.3f] (~1 expected)",
                       stats::median(rs), min(rs), max(rs)))
    for (ct in colnames(res)) {
      qs <- stats::quantile(res[, ct], c(0.25, 0.5, 0.75), na.rm = TRUE)
      logger$log("rundish",
                 sprintf("  %-12s median=%.3f IQR=[%.3f, %.3f]",
                         ct, qs[2], qs[1], qs[3]))
    }
  }
  rm(ref_data); gc(verbose = FALSE)
  res
}

#' Merge cell proportions into an existing sample sheet on disk
#'
#' Reads \code{file.path(dir, "sample_sheet.csv")}, removes any existing
#' cell-proportion columns (i.e. columns named exactly as in
#' \code{colnames(props)}), \code{left_join}s the new proportions on
#' \code{sample_id}, and writes the file back.
#'
#' @keywords internal
#' @noRd
merge_props_into_sheet <- function(dir, props, logger = NULL) {
  ss_path <- file.path(dir, "sample_sheet.csv")
  if (!file.exists(ss_path)) {
    if (!is.null(logger))
      logger$log("rundish",
                 sprintf("no sample_sheet.csv in %s; skipping merge", dir))
    return(invisible(FALSE))
  }
  ss <- utils::read.csv(ss_path, stringsAsFactors = FALSE, check.names = FALSE)
  drop <- intersect(colnames(ss), colnames(props))
  if (length(drop) > 0) {
    ss <- ss[, setdiff(colnames(ss), drop), drop = FALSE]
    if (!is.null(logger))
      logger$log("rundish",
                 sprintf("removed %d stale cell-prop column(s) from sample_sheet: %s",
                         length(drop), paste(drop, collapse = ", ")))
  }
  props_df <- data.frame(sample_id = rownames(props), props,
                         stringsAsFactors = FALSE, check.names = FALSE)
  ss <- dplyr::left_join(ss, props_df, by = "sample_id")
  utils::write.csv(ss, ss_path, row.names = FALSE)
  if (!is.null(logger))
    logger$log("rundish",
               sprintf("merged %d cell-prop columns into %s",
                       ncol(props), ss_path))
  invisible(TRUE)
}

#' Run EpiDISH cell-type deconvolution
#'
#' Single user-facing function with two calling forms. Dispatches on
#' the first argument:
#'
#' \describe{
#'   \item{\code{rundish(<matrix>, ...)}}{In-memory. Takes a beta
#'     matrix (probes x samples) and returns the
#'     sample-by-cell-type proportion matrix. No files written.}
#'   \item{\code{rundish("<dir>", ...)}}{On-disk standalone. Loads
#'     \code{betas_all.rds} (and \code{metadata.rds},
#'     \code{sample_sheet.csv}) from \code{dir}, optionally restricts
#'     to a subset of samples or to one cell-type label, runs the
#'     in-memory primitive, writes \code{cell_proportions.csv}, and
#'     left-joins the proportions into \code{sample_sheet.csv}. The
#'     proportion matrix is returned invisibly.}
#' }
#'
#' Reference: by default \code{mqcopts()$dishref} resolves to
#' \code{cent12CT.m} on EPIC and \code{cent12CT450k.m} on 450k (the
#' Salas et al. 2022 12-cell reference). Pass \code{ref =
#' "centDHSbloodDMC.m"} for the legacy 7-cell reference.
#'
#' There is no built-in tissue allowlist. The caller decides which
#' samples to run on. The function refuses (errors) only when the
#' reference's CpGs do not sufficiently overlap the input matrix
#' (< 50 probes), which is a data-level check.
#'
#' @param betas Either a numeric beta matrix (probes x samples) or
#'   a length-1 string giving a Stage 2 output directory path.
#' @param platform Platform string (e.g. \code{"EPIC"}, \code{"450k"}).
#'   Required when \code{ref} is \code{NULL} and the configured
#'   \code{dishref} is a list. For the on-disk form, defaults to the
#'   value in \code{metadata.rds}.
#' @param ref Reference dataset name (single string). Overrides
#'   \code{mqcopts()$dishref}.
#' @param method EpiDISH method. Default: \code{mqcopts()$dishmethod}
#'   (\code{"RPC"}).
#' @param samples (On-disk form only.) Character vector of sample IDs
#'   to restrict to. \code{NULL} = every sample in
#'   \code{betas_all.rds}.
#' @param celltype (On-disk form only.) Cell-type label to restrict to,
#'   matched case-insensitively against the cell-type column resolved
#'   via the \code{cellcol}/\code{cellaliases} options. \code{NULL} =
#'   no filter.
#' @param logger Optional logger from \code{\link{makelog}}. For the
#'   on-disk form, if NULL one is opened in the output directory.
#' @return The cell-proportion matrix (samples x cell types). For the
#'   on-disk form, returned invisibly along with the file side
#'   effects.
#' @examples
#' \dontrun{
#' # ---- In-memory ----
#' betas <- readRDS("results/PBMC/betas_all.rds")
#' props <- rundish(betas, platform = "EPIC")
#'
#' # ---- Standalone on a Stage 2 output directory ----
#' # Default: 12-cell reference, all samples
#' rundish("results/PBMC")
#'
#' # Restrict to one cell-type label within a top-level Stage 1 dir
#' rundish("results", celltype = "BM")        # "BM" here = B-memory
#'
#' # Restrict to a sample subset
#' rundish("results", samples = c("S1", "S5", "S12"))
#'
#' # Override reference (legacy 7-cell)
#' rundish("results/PBMC", ref = "centDHSbloodDMC.m")
#' }
#' @export
rundish <- function(betas,
                    platform = NULL,
                    ref      = NULL,
                    method   = NULL,
                    samples  = NULL,
                    celltype = NULL,
                    logger   = NULL) {

  # --- Form A: in-memory matrix -------------------------------------------
  if (is.matrix(betas) || is.data.frame(betas)) {
    if (!is.null(samples) || !is.null(celltype))
      warning("'samples' / 'celltype' are ignored for the in-memory form. ",
              "Subset the matrix before calling rundish().", call. = FALSE)
    return(rundish_matrix(as.matrix(betas),
                          platform = platform, ref = ref,
                          method = method, logger = logger))
  }

  # --- Form B: directory path ---------------------------------------------
  if (!(is.character(betas) && length(betas) == 1))
    stop("rundish() first argument must be a numeric matrix or a length-1 ",
         "directory path.", call. = FALSE)

  dir <- betas
  if (!dir.exists(dir))
    stop("rundish(): directory does not exist: ", dir, call. = FALSE)
  betas_path <- file.path(dir, "betas_all.rds")
  if (!file.exists(betas_path))
    stop("rundish(): no betas_all.rds in: ", dir, call. = FALSE)

  meta_path <- file.path(dir, "metadata.rds")
  ss_path   <- file.path(dir, "sample_sheet.csv")
  metadata  <- if (file.exists(meta_path)) readRDS(meta_path) else list()
  if (is.null(platform)) platform <- metadata$platform
  ss <- if (file.exists(ss_path))
          utils::read.csv(ss_path, stringsAsFactors = FALSE,
                          check.names = FALSE) else NULL

  if (is.null(logger)) logger <- makelog(dir)
  logger$log("rundish",
             sprintf("dir=%s, platform=%s",
                     dir, if (is.null(platform)) "<none>" else platform))

  mat <- readRDS(betas_path)

  keep_ids <- colnames(mat)
  if (!is.null(samples)) {
    keep_ids <- intersect(keep_ids, samples)
    logger$log("rundish",
               sprintf("samples filter: %d of %d requested matched",
                       length(keep_ids), length(samples)))
  }
  if (!is.null(celltype)) {
    if (is.null(ss))
      stop("celltype = '", celltype,
           "' requested but no sample_sheet.csv in ", dir, call. = FALSE)
    cfg <- mqcopts()
    cell_col <- resolve_column(colnames(ss), cfg$cellcol, cfg$cellaliases)
    if (is.na(cell_col))
      stop("No cell-type column resolved in sample sheet (looked for: ",
           cfg$cellcol, " + aliases).", call. = FALSE)
    ct_ids <- ss$sample_id[!is.na(ss[[cell_col]]) &
                           tolower(trimws(ss[[cell_col]])) ==
                           tolower(trimws(celltype))]
    keep_ids <- intersect(keep_ids, ct_ids)
    logger$log("rundish",
               sprintf("celltype='%s' filter via column '%s': %d samples",
                       celltype, cell_col, length(keep_ids)))
  }
  if (length(keep_ids) < 2)
    stop("rundish(): fewer than 2 samples after filtering.", call. = FALSE)
  mat <- mat[, keep_ids, drop = FALSE]

  if (is.null(ref)) ref <- resolve_dishref(mqcopts()$dishref, platform)
  logger$log("rundish",
             sprintf("running EpiDISH on %d samples, ref=%s",
                     ncol(mat), ref))

  props <- rundish_matrix(mat, platform = platform, ref = ref,
                          method = method, logger = logger)

  # Write proportions CSV
  out_csv <- file.path(dir, "cell_proportions.csv")
  utils::write.csv(data.frame(sample_id = rownames(props), props,
                              check.names = FALSE,
                              stringsAsFactors = FALSE),
                   out_csv, row.names = FALSE)
  logger$log("rundish", sprintf("wrote %s", out_csv))

  # Merge into sample sheet on disk (same as the in-pipeline path)
  merge_props_into_sheet(dir, props, logger = logger)

  rm(mat); gc(verbose = FALSE)
  invisible(props)
}
