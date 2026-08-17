## ---------------------------------------------------------------------------
## config.R -- package-wide options.
## ---------------------------------------------------------------------------

.mqc_env <- new.env(parent = emptyenv())

.mqc_defaults <- list(

  ## ---- sample sheet column resolution ------------------------------------
  ## Any delimited file whose name looks like a sample annotation. The 3.0.x
  ## pattern required the literal string "sample?sheet", so a file named
  ## samples.BM.PH1.txt -- a perfectly ordinary name -- was never found and the
  ## run silently proceeded with no metadata at all.
  sheetpattern  = "(?i)(sample|target|pheno|manifest|annot).*\\.(csv|tsv|txt|xlsx?)$",
  idcol         = "Sample_Name",
  idaliases     = c("Sample_Name", "sample_name", "Sample_ID", "sample_id",
                    "SampleID", "Name", "ID",
                    ## file-name style identifiers, common in lab-made sheets
                    "Fname", "FileName", "File_Name", "File", "Array",
                    "Basename", "Sentrix", "Sentrix_Barcode", "Barcode"),
  basecol       = "Basename",
  sexcol        = "Reported_Sex",
  sexaliases    = c("Reported_Sex", "Sex", "sex", "gender", "Gender",
                    "Reported_Gender"),
  agecol        = "Age",
  agealiases    = c("Age", "age", "age_years", "AgeYears", "Age_Years"),
  batchcol      = "batch_folder",
  batchaliases  = c("batch_folder", "Batch", "batch", "Plate", "plate",
                    "Sentrix_ID", "Slide"),
  donorcol      = "Donor",
  donoraliases  = c("Donor", "donor", "Subject", "subject", "Patient"),
  cellcol       = "Cell",
  cellaliases   = c("Cell", "cell", "Cell_Type", "cell_type", "CellType",
                    "Tissue", "tissue"),

  ## ---- compute -----------------------------------------------------------
  ## NULL means "inherit from qcplan()". Set explicitly to override.
  workers       = NULL,
  batch         = NULL,
  memfrac       = 0.80,     # fraction of available RAM the run may use
  extreme       = FALSE,    # extreme memory conservation (per-batch files)
  savesdf       = TRUE,     # retain SigDFs; downgraded automatically if large

  ## ---- detection and masking ---------------------------------------------
  detp          = 0.05,     # ELBAR detection p-value threshold
  maskuse       = "both",   # "both" | "detection" | "design" | "none"

  ## ---- sample-level QC thresholds ----------------------------------------
  samplemin     = 0.95,     # minimum call rate
  intmad        = 3,        # MAD multiplier for intensity outliers
  intfloor      = 1300,     # absolute cohort-level sanity floor (warning only)

  ## ---- probe-level QC thresholds -----------------------------------------
  failmin       = 0.10,     # probe failure rate above which a probe is listed
  inclqual      = FALSE,    # include design-masked probes in the failure plot

  ## ---- sex calling ---------------------------------------------------------
  ## See sexcall.R. `sexsep` is the separation floor below which the cohort is
  ## declared unimodal and NO sex is called; `sexband` is the orthogonal
  ## distance, in within-cluster SDs, outside which a sample is "unclear".
  sexsep        = 1.0,
  sexband       = 5.0,
  sexmin        = 8,        # cohort size below which sex is not called

  ## ---- EPICv2 replicate handling -----------------------------------------
  collapse       = NA,        # NA = automatic (TRUE on EPICv2)
  collapsemethod = "peters",  # "peters" | "minpval" | "mean"

  ## ---- cell composition ---------------------------------------------------
  dishref       = "blood",
  dishmethod    = "RPC",

  ## ---- SNP identity --------------------------------------------------------
  ## Pairwise concordance is quadratic in samples. Above this cohort size only
  ## pairs at or above `snpmin` are retained, so the output stays bounded.
  snpmin        = 0.70,
  snpmaxpairs   = 2e6
)

## The EpiDISH reference panels that actually exist in the installed package.
## "breast" was listed here in 3.0.0 and maps to centBreastSub.m, which is not
## an EpiDISH dataset, so the option silently disabled deconvolution.
.mqc_valid <- list(
  maskuse        = c("both", "detection", "design", "none"),
  collapsemethod = c("peters", "minpval", "mean"),
  dishmethod     = c("RPC", "CBS", "CP"),
  dishref        = c("blood", "bloodsub", "epithelial", "epifibfat", "12ct")
)

#' Default methylQC options
#'
#' @return a named list of the factory defaults.
#' @export
#' @examples
#' mqcdefaults()$detp
mqcdefaults <- function() .mqc_defaults

#' Get current methylQC options
#'
#' @param ... optional option names; if given, only those are returned.
#' @return a named list of the currently active options.
#' @export
#' @examples
#' mqcopts("detp", "maskuse")
mqcopts <- function(...) {
  ## Merged with `[k] <- list(v)` rather than utils::modifyList(), because
  ## modifyList() DELETES any element whose replacement is NULL. Several
  ## options (workers, batch) default to NULL, so mqcset(workers = NULL) --
  ## the restore idiom shown in mqcset()'s own example -- removed the option
  ## from the table entirely and made mqcopts("workers") throw "unknown
  ## option".
  cur <- .mqc_defaults
  ov <- as.list(.mqc_env, all.names = TRUE)
  for (k in names(ov)) cur[k] <- list(ov[[k]])

  nm <- c(...)
  if (length(nm) == 0L) return(cur)
  unknown <- setdiff(nm, names(cur))
  if (length(unknown))
    stop("unknown methylQC option(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  cur[nm]
}

#' Set methylQC options
#'
#' @param ... named options to set, e.g. \code{mqcset(detp = 0.01)}. Passing
#'   \code{NULL} restores that option to its default.
#' @return the previous values of the changed options, invisibly.
#' @export
#' @examples
#' old <- mqcset(detp = 0.01)
#' mqcopts("detp")$detp
#' mqcset(detp = old$detp)
mqcset <- function(...) {
  new <- list(...)
  if (!length(new)) return(invisible(list()))
  if (is.null(names(new)) || any(!nzchar(names(new))))
    stop("all arguments to mqcset() must be named", call. = FALSE)

  unknown <- setdiff(names(new), names(.mqc_defaults))
  if (length(unknown))
    stop("unknown methylQC option(s): ", paste(unknown, collapse = ", "),
         ". See mqcdefaults().", call. = FALSE)

  for (k in intersect(names(new), names(.mqc_valid))) {
    v <- new[[k]]
    if (!is.null(v) && !all(v %in% .mqc_valid[[k]]))
      stop("option '", k, "' must be one of: ",
           paste(.mqc_valid[[k]], collapse = ", "), call. = FALSE)
  }
  .mqc_check_numeric(new)

  old <- mqcopts()[names(new)]
  for (k in names(new)) assign(k, new[[k]], envir = .mqc_env)
  invisible(old)
}

#' Reset methylQC options to their defaults
#'
#' @return \code{NULL}, invisibly.
#' @export
#' @examples
#' mqcreset()
mqcreset <- function() {
  rm(list = ls(.mqc_env, all.names = TRUE), envir = .mqc_env)
  invisible(NULL)
}

## Range checks for the numeric options. Catches, e.g., detp = 5 (which would
## mask nothing) or memfrac = 8 (which would plan a run into swap).
## `allow_na` marks the options whose documented "disable this check" value is
## NA; without it, mqcset(intfloor = NA) was rejected while intflag()'s own
## documentation told the user to pass it.
.mqc_check_numeric <- function(new) {
  chk <- function(k, lo, hi, int = FALSE, allow_na = FALSE) {
    if (!k %in% names(new) || is.null(new[[k]])) return(invisible(NULL))
    v <- new[[k]]
    if (length(v) != 1L)
      stop("option '", k, "' must be a single number", call. = FALSE)
    ## A bare NA is logical, so the is.numeric() test has to come after the
    ## missingness test or mqcset(intfloor = NA) -- which intflag() documents --
    ## is rejected as "not a number".
    if (is.na(v)) {
      if (isTRUE(allow_na)) return(invisible(NULL))
      stop("option '", k, "' must not be missing", call. = FALSE)
    }
    if (!is.numeric(v))
      stop("option '", k, "' must be a single number", call. = FALSE)
    if (int && v != round(v))
      stop("option '", k, "' must be a whole number", call. = FALSE)
    if (v < lo || v > hi)
      stop("option '", k, "' must be between ", lo, " and ", hi, call. = FALSE)
    invisible(NULL)
  }
  chk("detp",        0,  1)
  chk("samplemin",   0,  1)
  chk("failmin",     0,  1)
  chk("memfrac",     0.05, 0.95)
  chk("intmad",      0,  20)
  chk("intfloor",    0,  Inf, allow_na = TRUE)
  chk("sexsep",      0,  Inf)
  chk("sexband",     0,  Inf)
  chk("sexmin",      2,  1e6, int = TRUE)
  chk("snpmin",      0,  1)
  chk("snpmaxpairs", 1,  Inf)
  chk("workers",     1,  1024, int = TRUE)
  chk("batch",       1,  1e6,  int = TRUE)
  invisible(NULL)
}

## Resolve an explicit argument against the option table: an explicit non-NULL
## argument always wins, otherwise fall back to the option.
.opt <- function(x, name, cfg = NULL) {
  if (!is.null(x)) return(x)
  cfg <- cfg %||% mqcopts()
  cfg[[name]]
}
