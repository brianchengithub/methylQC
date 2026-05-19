###############################################################################
# config.R — methylQC configuration system
#
# All tunable parameters stored as R options with prefix "methylQC.".
# Column names are configurable with case-insensitive alias matching.
###############################################################################

#' Default configuration for methylQC
#' @return A named list of defaults.
#' @export
methylQC_defaults <- function() {
  list(
    # ---- Sample sheet discovery ----
    sample_sheet_pattern = "sample.*\\.(csv|tsv|txt)$",
    basename_col         = "Basename",

    # ---- Key column names (user-configurable with aliases) ----
    sample_id_col        = "sample_id",
    sample_id_aliases    = c("sample_id", "Sample_ID", "SampleID",
                             "Sample_Name", "Sample", "sample_name"),
    donor_col            = "Donor",
    donor_aliases        = c("Donor", "donor", "Subject", "subject",
                             "Participant", "participant", "SubjectID",
                             "DonorID", "ID", "donor_id", "subject_id"),
    reported_sex_col     = "Reported_Sex",
    reported_sex_aliases = c("Reported_Sex", "Sex", "gender", "Gender",
                             "Female", "Male", "reported_sex", "sex"),
    reported_age_col     = "Age",
    reported_age_aliases = c("Age", "age", "age_years", "AgeYears",
                             "age_at_sample", "reported_age"),
    batch_col            = "batch_folder",
    batch_aliases        = c("batch_folder", "Batch", "batch", "Plate",
                             "plate", "Slide", "slide", "Chip", "chip",
                             "Array_Plate", "SentrixBarcode_A"),
    cell_col             = "Cell",
    cell_aliases         = c("Cell", "cell", "Cell_Type", "cell_type",
                             "CellType", "Tissue", "tissue",
                             "Tissue_Type", "tissue_type",
                             "Sample_Type", "sample_type",
                             "SampleType", "Source", "source"),

    # ---- Processing ----
    n_cores   = 1L,
    save_sdfs = TRUE,

    # ---- SeSAMe preprocessing ----
    # prep code passed to openSesame. "QCDPB" = Quality mask, Channel
    # inference, Dye-bias correction, detection P-values (pOOBAH), and
    # noob Background correction. pOOBAH (P) runs BEFORE noob (B), so
    # detection p-values are computed on non-noob signal, while the
    # output betas/M-values are noob-corrected.
    prep_code = "QCDPB",

    # ---- QC thresholds ----
    detP_thresh     = 0.05,
    sample_call_min = 0.95,
    probe_call_min  = 0.95,
    intensity_min   = 1300,

    # ---- QC plots ----
    # Number of most-variable probes used for MDS and PCA
    n_top_variable = 100000L,

    # ---- k-NN imputation (apply_mask) ----
    # Imputation uses the top knn_var_probes most-variable probes to
    # compute sample-to-sample Euclidean distance, then imputes each
    # missing value as the inverse-distance-weighted mean of the
    # knn_k nearest neighbours.
    knn_k          = 50L,
    knn_var_probes = 30000L,

    # ---- EpiDISH ----
    # Default reference is the 12 immune cell-type panel (cent12CT.m),
    # used as a purity check for both whole blood and isolated/sorted
    # blood cell populations.
    epidish_reference  = "cent12CT.m",
    epidish_method     = "RPC",
    epidish_cell_types = c("PBMC", "WB", "WBC", "buffy", "buffy coat",
                           "whole blood", "blood",
                           "CD4CM", "CD4EM", "CD4N", "CD8CM", "CD8EM",
                           "CD8N", "BN", "BM", "NK", "Mono", "Granulo",
                           "B", "B cell", "T cell", "monocyte",
                           "granulocyte", "neutrophil", "eosinophil",
                           "basophil", "lymphocyte", "leukocyte"),

    # ---- SeSAMe replicate probe collapsing (legacy openSesame path) ----
    collapse_to_pfx = FALSE,
    collapse_method = "mean",

    # ---- EPICv2 de-duplication + probe ID harmonization ----
    # When TRUE: for EPICv2 data, replicate probes (sharing a cg/ch
    # prefix) are de-duplicated by selecting the probe with the fewest
    # cross-sample detection failures (ties -> first probe), then probe
    # IDs are harmonized to the target platform via sesameData mappings.
    epicv2_harmonize = FALSE,
    epicv2_target    = "EPIC"
  )
}

#' Get current methylQC configuration
#' @param ... Named overrides to apply temporarily.
#' @return A named list of all current options.
#' @export
methylQC_options <- function(...) {
  defaults <- methylQC_defaults()
  current <- lapply(names(defaults), function(k) {
    getOption(paste0("methylQC.", k), defaults[[k]])
  })
  names(current) <- names(defaults)
  overrides <- list(...)
  if (length(overrides)) current[names(overrides)] <- overrides
  current
}

#' Set methylQC options globally
#' @param ... Named options (e.g., \code{donor_col = "SubjectID"}).
#' @return Invisibly, the previous values.
#' @export
methylQC_set <- function(...) {
  overrides <- list(...)
  if (length(overrides) == 0) return(invisible(list()))
  if (is.null(names(overrides)) || any(names(overrides) == "")) {
    stop("All arguments to methylQC_set() must be named.", call. = FALSE)
  }
  valid <- names(methylQC_defaults())
  unknown <- setdiff(names(overrides), valid)
  if (length(unknown)) {
    warning("Unknown methylQC options: ", paste(unknown, collapse = ", "),
            call. = FALSE)
    overrides <- overrides[intersect(names(overrides), valid)]
  }
  if (length(overrides) == 0) return(invisible(list()))
  prefixed <- stats::setNames(overrides, paste0("methylQC.", names(overrides)))
  old <- do.call(options, prefixed)
  invisible(old)
}

#' Reset all methylQC options to defaults
#' @export
methylQC_reset <- function() {
  defaults <- methylQC_defaults()
  prefixed <- stats::setNames(vector("list", length(defaults)),
                              paste0("methylQC.", names(defaults)))
  do.call(options, prefixed)
  invisible(NULL)
}

#' Resolve a column name via preferred + aliases (case-insensitive)
#' @keywords internal
#' @noRd
resolve_column <- function(columns, preferred, aliases = character(0)) {
  candidates <- unique(c(preferred, aliases))
  for (cand in candidates) {
    hit <- columns[tolower(columns) == tolower(cand)]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

#' Normalize sex encodings to M/F/NA
#' @keywords internal
#' @noRd
normalize_sex <- function(x, col_name = "") {
  if (is.null(x)) return(character(0))
  result <- rep(NA_character_, length(x))
  is_numeric_binary <- all(stats::na.omit(as.character(x)) %in% c("0", "1"))
  lower_name <- tolower(col_name)
  if (is_numeric_binary && grepl("female", lower_name)) {
    result <- ifelse(x == 1, "F", ifelse(x == 0, "M", NA_character_))
  } else if (is_numeric_binary && grepl("^male$|_male$", lower_name)) {
    result <- ifelse(x == 1, "M", ifelse(x == 0, "F", NA_character_))
  } else {
    upper <- toupper(as.character(x))
    result <- ifelse(upper %in% c("M", "MALE"), "M",
              ifelse(upper %in% c("F", "FEMALE"), "F", NA_character_))
  }
  result
}

.onLoad <- function(libname, pkgname) {
  defaults <- methylQC_defaults()
  for (k in names(defaults)) {
    opt_name <- paste0("methylQC.", k)
    if (is.null(getOption(opt_name))) {
      do.call(options, stats::setNames(list(defaults[[k]]), opt_name))
    }
  }
}

.onAttach <- function(libname, pkgname) {
  v <- tryCatch(as.character(utils::packageVersion(pkgname)),
                error = function(e) "")
  packageStartupMessage(
    sprintf("methylQC %s -- run check_dependencies() to verify your environment.", v))
}
