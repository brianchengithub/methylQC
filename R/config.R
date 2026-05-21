###############################################################################
# config.R — methylQC configuration system
#
# All tunable parameters stored as R options with prefix "methylQC.".
# Column names are configurable with case-insensitive alias matching.
#
# Option keys are short, separator-free identifiers (see mqcdefaults()).
###############################################################################

#' Default configuration for methylQC
#'
#' Returns the full set of tunable options. Option keys are short,
#' separator-free identifiers. Override globally with [mqcset()] or
#' temporarily by passing named arguments to most pipeline functions.
#'
#' @return A named list of defaults.
#' @export
mqcdefaults <- function() {
  list(
    # ---- Sample sheet discovery ----
    sheetpattern = "sample.*\\.(csv|tsv|txt)$",
    basecol      = "Basename",

    # ---- Key column names (user-configurable with aliases) ----
    idcol        = "sample_id",
    idaliases    = c("sample_id", "Sample_ID", "SampleID",
                     "Sample_Name", "Sample", "sample_name"),
    donorcol     = "Donor",
    donoraliases = c("Donor", "donor", "Subject", "subject",
                     "Participant", "participant", "SubjectID",
                     "DonorID", "ID", "donor_id", "subject_id"),
    sexcol       = "Reported_Sex",
    sexaliases   = c("Reported_Sex", "Sex", "gender", "Gender",
                     "Female", "Male", "reported_sex", "sex"),
    agecol       = "Age",
    agealiases   = c("Age", "age", "age_years", "AgeYears",
                     "age_at_sample", "reported_age"),
    batchcol     = "batch_folder",
    batchaliases = c("batch_folder", "Batch", "batch", "Plate",
                     "plate", "Slide", "slide", "Chip", "chip",
                     "Array_Plate", "SentrixBarcode_A"),
    cellcol      = "Cell",
    cellaliases  = c("Cell", "cell", "Cell_Type", "cell_type",
                     "CellType", "Tissue", "tissue",
                     "Tissue_Type", "tissue_type",
                     "Sample_Type", "sample_type",
                     "SampleType", "Source", "source"),

    # ---- Processing ----
    cores   = 1L,
    savesdf = FALSE,

    # ---- QC thresholds ----
    # detp      : detection p-value cutoff. A probe PASSES at detP <= detp.
    # samplemin : minimum per-sample call rate (frac_dt) before low-detection flag.
    # intmin    : minimum per-sample mean intensity before low-intensity flag.
    # probemin  : per-probe pass-rate used ONLY for the console "Probe breakdown"
    #             summary. It is informational and drives no file or exclusion.
    # failmin   : per-probe sample-FAILURE fraction at/above which a probe is
    #             written to failed_probes.csv. The Page-2 plot draws a dashed
    #             vertical line at this value. Default 0.10 matches meffil's
    #             detectionp.cpgs.threshold; liberal end of the defensible
    #             EWAS-package range (ChAMP 0%, minfi 0%, DNAmArray 5%,
    #             meffil 10%). Appropriate when k-NN imputation is in the
    #             downstream workflow and cohort size is large enough that
    #             10% of samples still provides robust neighbours.
    detp      = 0.05,
    samplemin = 0.95,
    intmin    = 1300,
    probemin  = 0.95,
    failmin   = 0.10,

    # inclqual : if TRUE, quality-masked probes are kept in the probe-failure
    #            plot/CSV. Default FALSE excludes them (they fail by design).
    inclqual  = FALSE,

    # ---- QC plots ----
    # ntop : number of most-variable probes used for PCA (MDS uses all probes).
    ntop = 100000L,

    # ---- EpiDISH ----
    dishref    = "centDHSbloodDMC.m",
    dishmethod = "RPC",
    bloodtypes = c("PBMC", "WB", "WBC", "buffy", "buffy coat",
                   "whole blood", "blood"),

    # ---- SeSAMe ----
    collapse       = FALSE,
    collapsemethod = "mean"
  )
}

#' Get current methylQC configuration
#' @param ... Named overrides to apply temporarily to the returned list.
#' @return A named list of all current options.
#' @export
mqcopts <- function(...) {
  defaults <- mqcdefaults()
  current <- lapply(names(defaults), function(k) {
    getOption(paste0("methylQC.", k), defaults[[k]])
  })
  names(current) <- names(defaults)
  overrides <- list(...)
  if (length(overrides)) current[names(overrides)] <- overrides
  current
}

#' Set methylQC options globally
#' @param ... Named options (e.g., \code{donorcol = "SubjectID"}).
#' @return Invisibly, the previous values.
#' @export
mqcset <- function(...) {
  overrides <- list(...)
  if (length(overrides) == 0) return(invisible(list()))
  if (is.null(names(overrides)) || any(names(overrides) == "")) {
    stop("All arguments to mqcset() must be named.", call. = FALSE)
  }
  valid <- names(mqcdefaults())
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
#' @return Invisibly NULL.
#' @export
mqcreset <- function() {
  defaults <- mqcdefaults()
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

# NOTE: there is intentionally NO .onLoad() block.
#
# An earlier version of this file pre-populated every option via
# options(methylQC.* = ...) on package load. That meant any subsequent
# change to mqcdefaults() (e.g. tightening a threshold for a new
# methylQC release) was IGNORED for the rest of the R session, because
# the option was already non-NULL. mqcopts() returned the stale value,
# and the change went silently into the void.
#
# Without .onLoad(), session options stay NULL until the user explicitly
# calls mqcset(). mqcopts() always reads
#   getOption("methylQC.<key>", mqcdefaults()$<key>)
# so unset keys cleanly fall through to the current default each call.
