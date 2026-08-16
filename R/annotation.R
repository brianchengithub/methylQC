## ---------------------------------------------------------------------------
## annotation.R -- platform annotation, resolved ONCE in the parent process.
##
## This file exists because of a specific defect in sesame's parallel path.
## sesame::readIDATpair() has:
##
##     if (is.null(manifest)) {
##         df_address <- sesameDataGet(paste0(attr(dm,'platform'), '.address'))
##         manifest <- df_address$ordering
##         controls <- df_address$controls
##     }
##
## openSesame() defaults manifest = NULL and passes it through to every
## bplapply worker, so each forked worker independently calls sesameDataGet().
## That routes into ExperimentHub / BiocFileCache, which is SQLite, and N
## forked processes contending for one SQLite file deadlock. The symptom is a
## freeze with no error and no CPU activity.
##
## The fix is to resolve the address object once here, in the parent, and hand
## both `manifest` and `controls` to every worker explicitly. Note that they
## must be passed TOGETHER: `controls` is only assigned inside that same
## `if (is.null(manifest))` block, so supplying the manifest alone leaves
## controls NULL, attr(sdf, "controls") is never set, and noob() fails at
## negControls() because combine.neg defaults to TRUE.
## ---------------------------------------------------------------------------

.anno_cache <- new.env(parent = emptyenv())

#' Resolve platform annotation once
#'
#' Fetches the address object (manifest ordering plus control probe table) for
#' a platform and memoises it for the session. Call this in the parent process
#' before dispatching any parallel work.
#'
#' @param platform platform string, e.g. \code{"EPIC"}, \code{"EPICv2"},
#'   \code{"HM450"}.
#' @param logger optional logger.
#' @return a list with \code{ordering} and \code{controls}.
#' @export
#' @examples
#' \dontrun{
#' addr <- mqcanno("EPIC")
#' }
mqcanno <- function(platform, logger = NULL) {
  stopifnot(is.character(platform), length(platform) == 1L, nzchar(platform))
  key <- paste0("addr:", platform)
  if (!is.null(.anno_cache[[key]])) return(.anno_cache[[key]])

  logger <- logger %||% nulllog()
  if (!requireNamespace("sesameData", quietly = TRUE))
    stop("sesameData is required. install.packages(\"pak\"); ",
         "pak::pkg_install(\"bioc::sesameData\")", call. = FALSE)

  logger$log("anno", sprintf("resolving %s annotation (once, in parent)", platform))
  addr <- tryCatch(
    sesameData::sesameDataGet(paste0(platform, ".address")),
    error = function(e)
      stop("could not load annotation for platform '", platform, "': ",
           conditionMessage(e),
           "\nIf this is the first run, cache the data with ",
           "sesameData::sesameDataCache() and try again.", call. = FALSE))

  if (is.null(addr$ordering))
    stop("annotation for '", platform, "' has no $ordering component",
         call. = FALSE)
  if (is.null(addr$controls) || !nrow(addr$controls))
    warning("annotation for '", platform, "' has no control probes; ",
            "noob background correction may fail.", call. = FALSE)

  out <- list(ordering = addr$ordering, controls = addr$controls)
  assign(key, out, envir = .anno_cache)
  out
}

#' Pre-warm the sesameData cache
#'
#' Calls \code{sesameData::sesameDataCache()} once. Harmless if already cached.
#'
#' @param logger optional logger
#' @return \code{TRUE} on success, \code{FALSE} otherwise
#' @keywords internal
#' @noRd
mqcprewarm <- function(logger = NULL) {
  logger <- logger %||% nulllog()
  ok <- tryCatch({
    suppressMessages(sesameData::sesameDataCache())
    TRUE
  }, error = function(e) {
    logger$log("anno", paste("sesameDataCache() failed:", conditionMessage(e)))
    FALSE
  })
  ok
}

#' Read one IDAT pair with annotation supplied explicitly
#'
#' Thin wrapper over \code{sesame::readIDATpair} that always passes both
#' \code{manifest} and \code{controls}, so no worker ever touches
#' ExperimentHub.
#'
#' @param pfx IDAT basename (path without the _Grn.idat / _Red.idat suffix)
#' @param platform platform string
#' @param addr the list returned by \code{\link{mqcanno}}
#' @return a \code{SigDF}
#' @keywords internal
#' @noRd
read_idat <- function(pfx, platform, addr) {
  sesame::readIDATpair(pfx,
                       platform = platform,
                       manifest = addr$ordering,
                       controls = addr$controls)
}


## ---------------------------------------------------------------------------
## Design mask
## ---------------------------------------------------------------------------

#' Probe design mask for a platform
#'
#' The design mask flags probes with known design problems: multi-mapping
#' probes, probes overlapping common SNPs, and cross-hybridising probes. It is
#' a platform lookup, so it is identical for every sample on that platform. v2
#' stored it fused with the per-sample detection mask in a probes-by-samples
#' logical matrix; v3 stores it once as a named logical vector, which is about
#' 3.5 MB for EPIC instead of 4.8 GB at 1500 samples, and keeps the two mask
#' types independently addressable.
#'
#' @param sdf a \code{SigDF}, used only to supply the probe set and platform.
#' @return a named logical vector, \code{TRUE} where the probe is masked by
#'   design.
#' @export
#' @examples
#' \dontrun{
#' dm <- designmask(sdf)
#' sum(dm)
#' }
designmask <- function(sdf) {
  m <- sesame::qualityMask(sesame::resetMask(sdf))
  stats::setNames(as.logical(m$mask), m$Probe_ID)
}


## ---------------------------------------------------------------------------
## Sex chromosome probes  (memoised; v2 re-fetched the manifest at 4 sites)
## ---------------------------------------------------------------------------

#' Sex chromosome probe IDs for a platform
#'
#' Memoised for the session. v2 called \code{sesameDataGet()} afresh from
#' masking.R, pipeline.R and plots.R, re-fetching and re-parsing the manifest
#' each time.
#'
#' @param platform platform string
#' @param logger optional logger
#' @return a character vector of probe IDs on chrX or chrY; \code{character(0)}
#'   if the manifest cannot be read.
#' @export
#' @examples
#' \dontrun{
#' length(sexprobes("EPIC"))
#' }
sexprobes <- function(platform, logger = NULL) {
  key <- paste0("sex:", platform)
  if (!is.null(.anno_cache[[key]])) return(.anno_cache[[key]])

  logger <- logger %||% nulllog()
  ## GenomeInfoDb is only needed for seqnames() on the manifest GRanges. It is
  ## a Suggests, so say plainly when it is absent: 3.0.0 swallowed the failure
  ## and returned character(0), which silently left every chrX and chrY probe
  ## in the PCA input while the report claimed they had been excluded.
  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    logger$log("anno", paste(
      "GenomeInfoDb is not installed, so sex-chromosome probes cannot be",
      "identified and will NOT be excluded from the PCA. Install with",
      "pak::pkg_install(\"bioc::GenomeInfoDb\")."), warn = TRUE)
    return(character(0))
  }
  ids <- tryCatch({
    mf <- sesameData::sesameData_getManifestGRanges(platform)
    chr <- as.character(GenomeInfoDb_seqnames(mf))
    names(mf)[chr %in% c("chrX", "chrY", "X", "Y")]
  }, error = function(e) {
    logger$log("anno",
               paste("could not read sex probes for", platform, ":",
                     conditionMessage(e)), warn = TRUE)
    character(0)
  })
  ids <- unique(stats::na.omit(ids))
  assign(key, ids, envir = .anno_cache)
  ids
}

## GenomeInfoDb is pulled in transitively by sesameData; avoid a hard Import
## just for seqnames().
#' @keywords internal
#' @noRd
GenomeInfoDb_seqnames <- function(x) {
  f <- get("seqnames", envir = asNamespace("GenomeInfoDb"))
  f(x)
}

#' Clear the annotation cache
#'
#' @return \code{NULL}, invisibly.
#' @export
#' @examples
#' mqcannoreset()
mqcannoreset <- function() {
  rm(list = ls(.anno_cache, all.names = TRUE), envir = .anno_cache)
  invisible(NULL)
}
