## ---------------------------------------------------------------------------
## dependencies.R -- runtime dependency and platform checks.
## ---------------------------------------------------------------------------

.MQC_MIN_SESAME <- "1.20.0"   # first Bioconductor release carrying ELBAR
                              # (ELBAR and prep code "I" appear from 1.18.x;
                              #  1.20.0 is the first release we test against)

#' Check that required packages and the platform are usable
#'
#' Verifies that sesame is installed at a version carrying ELBAR, that
#' sesameData is present, and that the operating system supports forked
#' parallelism. Optional packages are reported but not required.
#'
#' @param parallel if \code{TRUE}, also check that forking is available.
#' @param stop_on_missing raise an error rather than a warning when a required
#'   package is absent.
#' @return a list describing what was found, invisibly.
#' @export
#' @examples
#' \dontrun{
#' checkdeps()
#' }
checkdeps <- function(parallel = TRUE, stop_on_missing = TRUE) {
  out <- list(ok = TRUE, missing = character(0), notes = character(0))

  need <- function(pkg, minver = NULL) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      out$missing <<- c(out$missing, pkg)
      out$ok <<- FALSE
      return(invisible(FALSE))
    }
    if (!is.null(minver)) {
      have <- as.character(utils::packageVersion(pkg))
      if (utils::compareVersion(have, minver) < 0) {
        out$ok <<- FALSE
        out$notes <<- c(out$notes,
          sprintf("%s %s is installed but >= %s is required", pkg, have, minver))
      }
    }
    invisible(TRUE)
  }

  need("sesame", .MQC_MIN_SESAME)
  need("sesameData")
  need("BiocParallel")
  need("ggplot2")
  # matrixStats is optional: utils.R falls back to base R

  if (isTRUE(parallel) && identical(.Platform$OS.type, "windows")) {
    out$notes <- c(out$notes,
      "Windows detected: forked parallelism is unavailable, running serially.")
  }

  ## A package present but below its minimum sets out$ok FALSE and adds a note.
  ## 3.0.0 only ever raised on an ABSENT package, so an old sesame without
  ## ELBAR passed the check and failed later inside a worker.
  if (!isTRUE(out$ok) && !length(out$missing) && length(out$notes)) {
    msg <- paste(c("methylQC dependency check failed:", out$notes), collapse = "\n  ")
    if (isTRUE(stop_on_missing)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }

  if (length(out$missing)) {
    msg <- paste0(
      "missing required package(s): ", paste(out$missing, collapse = ", "),
      "\nInstall with:\n",
      "  install.packages(\"pak\")\n",
      "  pak::pkg_install(c(",
      paste(sprintf("\"%s\"", .mqc_pakname(out$missing)), collapse = ", "), "))")
    if (isTRUE(stop_on_missing)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  if (length(out$notes)) for (n in out$notes) message("methylQC: ", n)

  invisible(out)
}

## Bioconductor packages need the bioc:: prefix for pak.
.MQC_BIOC <- c("sesame", "sesameData", "BiocParallel", "EpiDISH", "EPICv2manifest")
.mqc_pakname <- function(pkgs) {
  ifelse(pkgs %in% .MQC_BIOC, paste0("bioc::", pkgs), pkgs)
}

#' Report whether an optional package is available
#'
#' @param pkg package name
#' @param what what the package is needed for, used in the message
#' @param logger optional logger
#' @return \code{TRUE} if usable
#' @keywords internal
#' @noRd
have_pkg <- function(pkg, what = NULL, logger = NULL) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) {
    msg <- sprintf("optional package '%s' is not installed%s; install with pak::pkg_install(\"%s\")",
                   pkg, if (is.null(what)) "" else paste0(", skipping ", what),
                   .mqc_pakname(pkg))
    (logger %||% nulllog())$log("deps", msg)
  }
  ok
}

#' Are we on a platform that supports forked parallelism?
#'
#' @return \code{TRUE} on macOS and Linux.
#' @keywords internal
#' @noRd
can_fork <- function() !identical(.Platform$OS.type, "windows")
