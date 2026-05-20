###############################################################################
# logging.R — Simple file + console logger
#
# Provides timestamped logging to both the R console (via message()) and
# a persistent log file (pipeline_diagnostics.log). Each log entry includes
# a timestamp, a step label (e.g. "runsesame", "qcflags"), and the message
# text. The log file is append-only and persists across pipeline runs.
###############################################################################
#' Create a file + console logger
#'
#' @param outdir Directory where the log file will be written.
#' @return A list with \code{log(step, msg)} function and \code{path}.
#' @export
makelog <- function(outdir) {
  log_path <- file.path(outdir, "pipeline_diagnostics.log")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(log_path)) {
    cat("=== Pipeline log started:", as.character(Sys.time()), "===\n",
        file = log_path)
  }

  log_fn <- function(step, msg) {
    line <- sprintf("[%s] %s :: %s",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), step, msg)
    message(line)
    cat(line, "\n", file = log_path, append = TRUE)
  }

  list(log = log_fn, path = log_path)
}

#' Append a captured object to the logger's file
#' @param logger A logger from \code{\link{makelog}}.
#' @param obj Any object; printed via \code{print()}.
#' @export
logcapture <- function(logger, obj) {
  utils::capture.output(print(obj), file = logger$path, append = TRUE)
}
