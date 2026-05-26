###############################################################################
# dependencies.R — Environment and dependency verification
#
# methylQC does not silently install or update packages: doing so is
# invasive, breaks on no-admin / HPC environments, and is disallowed by
# CRAN/Bioconductor policy. Instead, checkdeps() verifies that the
# required packages are present at adequate versions and prints the
# exact install command for anything missing or stale.
###############################################################################

#' Minimum required versions for key dependencies
#' @keywords internal
#' @noRd
.methylqc_min_versions <- list(
  sesame     = "1.20.0",
  sesameData = "1.20.0",
  EpiDISH    = "2.18.0",
  ggplot2    = "3.4.0",
  matrixStats = "0.62.0"
)

#' Verify the methylQC runtime environment
#'
#' Checks R version, required package availability and versions, the
#' SeSAMe data cache, and the presence of the default EpiDISH reference
#' panel. Prints a report and, for anything missing or out of date, the
#' exact command to fix it. Nothing is installed automatically.
#'
#' @param quiet If TRUE, suppress the printed report and only return the
#'   status list invisibly.
#' @return Invisibly, a list with elements \code{ok} (logical) and
#'   \code{problems} (character vector).
#' @export
checkdeps <- function(quiet = FALSE) {
  problems <- character(0)
  note <- function(...) if (!quiet) message(...)

  note("methylQC dependency check")
  note("-------------------------")

  # --- R version ---
  rv <- getRversion()
  note(sprintf("  R version: %s", rv))
  if (rv < "4.1.0") {
    problems <- c(problems, "R >= 4.1.0 is required.")
  }

  # --- Required packages + versions ---
  for (pkg in names(.methylqc_min_versions)) {
    minv <- .methylqc_min_versions[[pkg]]
    if (!requireNamespace(pkg, quietly = TRUE)) {
      note(sprintf("  %-12s MISSING (need >= %s)", pkg, minv))
      problems <- c(problems, sprintf(
        "%s is not installed. Install with: BiocManager::install(\"%s\")",
        pkg, pkg))
      next
    }
    have <- utils::packageVersion(pkg)
    if (have < minv) {
      note(sprintf("  %-12s %s  (need >= %s)  OUTDATED", pkg, have, minv))
      problems <- c(problems, sprintf(
        "%s %s is older than the required %s. Update with: BiocManager::install(\"%s\")",
        pkg, have, minv, pkg))
    } else {
      note(sprintf("  %-12s %s  OK", pkg, have))
    }
  }

  # --- SeSAMe data cache ---
  cache_ok <- tryCatch({
    if (requireNamespace("sesameData", quietly = TRUE)) {
      # A successful get of a small manifest indicates the cache is live
      invisible(sesameData::sesameDataGet("genomeInfo.hg38"))
      TRUE
    } else FALSE
  }, error = function(e) FALSE)
  if (cache_ok) {
    note("  sesameData cache: OK")
  } else {
    note("  sesameData cache: NOT INITIALIZED")
    problems <- c(problems, paste0(
      "The SeSAMe data cache is not initialized. Run once:\n",
      "    sesameData::sesameDataCacheAll()"))
  }

  # --- Default EpiDISH reference panel ---
  if (requireNamespace("EpiDISH", quietly = TRUE)) {
    dr <- mqcopts()$dishref
    # Flatten to a character vector of reference names to check.
    refs <- if (is.list(dr)) unique(unlist(dr, use.names = FALSE))
            else             unique(as.character(dr))
    for (ref1 in refs) {
      ref_ok <- tryCatch({
        e <- new.env()
        suppressWarnings(utils::data(list = ref1, package = "EpiDISH",
                                     envir = e))
        exists(ref1, envir = e)
      }, error = function(e) FALSE)
      if (ref_ok) {
        note(sprintf("  EpiDISH reference '%s': OK", ref1))
      } else {
        note(sprintf("  EpiDISH reference '%s': NOT FOUND", ref1))
        problems <- c(problems, sprintf(paste0(
          "EpiDISH reference '%s' is not available in your EpiDISH ",
          "install. Update EpiDISH: BiocManager::install(\"EpiDISH\")"),
          ref1))
      }
    }
  }

  ok <- length(problems) == 0
  if (!quiet) {
    note("-------------------------")
    if (ok) {
      note("All checks passed.")
    } else {
      note(sprintf("%d issue(s) found:", length(problems)))
      for (p in problems) note(paste0("  - ", p))
    }
  }
  invisible(list(ok = ok, problems = problems))
}
