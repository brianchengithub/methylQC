## ---------------------------------------------------------------------------
## logging.R -- a tiny append-only logger.
## ---------------------------------------------------------------------------

#' Create a pipeline logger
#'
#' Returns an object with a \code{$log(stage, message, warn = FALSE)} method
#' that appends a timestamped line to the log file and echoes it to the
#' console. When \code{warn = TRUE} an R \code{warning()} is also raised, so
#' failures are visible in an interactive session rather than only in the file.
#'
#' @param path log file path; the containing directory is created.
#' @param echo print each line to the console as well.
#' @return a list with \code{$log()}, \code{$path} and \code{$close()}.
#' @export
#' @examples
#' lg <- makelog(file.path(tempdir(), "demo.log"))
#' lg$log("demo", "hello")
#' lg$close()
makelog <- function(path, echo = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "at")

  log_fn <- function(stage, message, warn = FALSE) {
    line <- sprintf("[%s] %-12s %s",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stage, message)
    try(writeLines(line, con), silent = TRUE)
    try(flush(con), silent = TRUE)
    if (isTRUE(echo)) message(line)
    if (isTRUE(warn)) warning(sprintf("[%s] %s", stage, message), call. = FALSE)
    invisible(line)
  }

  list(
    log   = log_fn,
    path  = path,
    close = function() try(close(con), silent = TRUE)
  )
}

## A logger that discards everything, so that internal functions can take
## `logger = NULL` without every call site guarding with is.null().
#' @keywords internal
#' @noRd
nulllog <- function() {
  list(log = function(stage, message, warn = FALSE) {
        if (isTRUE(warn)) warning(sprintf("[%s] %s", stage, message), call. = FALSE)
        invisible(NULL)
       },
       path = NA_character_,
       close = function() invisible(NULL))
}

#' Capture output, messages and warnings from an expression into a log
#'
#' @param expr expression to evaluate
#' @param logger a logger from \code{\link{makelog}}
#' @param stage stage label recorded against each captured line
#' @return the value of \code{expr}
#' @export
#' @examples
#' lg <- makelog(file.path(tempdir(), "demo2.log"))
#' logcapture(message("working"), lg, "demo")
#' lg$close()
logcapture <- function(expr, logger, stage = "capture") {
  logger <- logger %||% nulllog()
  withCallingHandlers(
    expr,
    message = function(m) {
      logger$log(stage, sub("\n$", "", conditionMessage(m)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      logger$log(stage, paste("WARNING:", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
}
