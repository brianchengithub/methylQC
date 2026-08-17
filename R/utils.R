## ---------------------------------------------------------------------------
## utils.R -- small helpers shared across the package.
##
## Everything here is dependency-free base R so that it can be used by any
## other file without worrying about load order. There is exactly ONE
## definition of each helper in the package (v2 defined `%||%` twice, in
## preprocess.R and pipeline.R, which R CMD check flags).
## ---------------------------------------------------------------------------

#' Null-coalescing operator
#'
#' Returns \code{a} unless it is \code{NULL}, in which case \code{b}.
#'
#' @param a value to test
#' @param b fallback used when \code{a} is \code{NULL}
#' @return \code{a} if not \code{NULL}, otherwise \code{b}
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Null-or-missing coalescing operator
#'
#' Like \code{\%||\%} but also falls back when \code{a} is a length-one
#' \code{NA}. Several base R probes (\code{parallel::detectCores()},
#' \code{sysmem()$available}) report failure as \code{NA} rather than
#' \code{NULL}, and \code{\%||\%} passes those straight through into arithmetic
#' and then into \code{if ()}, which errors with "missing value where
#' TRUE/FALSE needed".
#'
#' @param a value to test
#' @param b fallback used when \code{a} is \code{NULL} or \code{NA}
#' @return \code{a} if usable, otherwise \code{b}
#' @keywords internal
#' @noRd
`%|NA|%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 1L && is.na(a)) return(b)
  a
}


## ---------------------------------------------------------------------------
## Output directory layout
## ---------------------------------------------------------------------------

## The v3 layout. Defined in exactly one place so that changing it does not
## require editing string literals scattered across pipeline.R and plots.R.
##
##   outdir/
##     METHODS.txt
##     logs/pipeline.log
##     data/matrices/{betas_all,detP_all,design_mask,snp_betas,sdfs_all}.rds
##     data/matrices/batches/          (extreme = TRUE only)
##     data/metadata/{sample_sheet,failed_samples,pc_scores,failed_probes}.csv
##     qc/{qccache.rds,qc_plots.pdf}

.MQC_LAYOUT <- list(
  root         = c(),
  logs         = c("logs"),
  matrices     = c("data", "matrices"),
  batches      = c("data", "matrices", "batches"),
  metadata     = c("data", "metadata"),
  qc           = c("qc"),

  methods      = c("METHODS.txt"),
  log          = c("logs", "pipeline.log"),
  betas        = c("data", "matrices", "betas_all.rds"),
  detp         = c("data", "matrices", "detP_all.rds"),
  designmask   = c("data", "matrices", "design_mask.rds"),
  snpbetas     = c("data", "matrices", "snp_betas.rds"),
  sdfs         = c("data", "matrices", "sdfs_all.rds"),
  sheet        = c("data", "metadata", "sample_sheet.csv"),
  failedsample = c("data", "metadata", "failed_samples.csv"),
  failedprobe  = c("data", "metadata", "failed_probes.csv"),
  pcscores     = c("data", "metadata", "pc_scores.csv"),
  snpconc      = c("data", "metadata", "snp_concordance.csv"),
  runinfo      = c("data", "metadata", "run_info.rds"),
  cache        = c("qc", "qccache.rds"),
  plots        = c("qc", "qc_plots.pdf")
)

#' Resolve a path inside a methylQC output directory
#'
#' All reads and writes in the package go through this function so that the
#' on-disk layout is defined in a single place.
#'
#' @param dir the output directory.
#' @param what a key from the layout table, e.g. \code{"betas"}, \code{"sheet"},
#'   \code{"matrices"}. Directory keys return the directory itself.
#' @param create if \code{TRUE}, create the containing directory.
#' @return an absolute or relative file path (character scalar).
#' @export
#' @examples
#' mqcpath(tempdir(), "betas")
mqcpath <- function(dir, what, create = FALSE) {
  stopifnot(is.character(dir), length(dir) == 1L)
  if (!what %in% names(.MQC_LAYOUT)) {
    stop("unknown path key '", what, "'. Known keys: ",
         paste(names(.MQC_LAYOUT), collapse = ", "), call. = FALSE)
  }
  parts <- .MQC_LAYOUT[[what]]
  p <- if (length(parts) == 0L) dir else do.call(file.path, c(list(dir), as.list(parts)))
  if (isTRUE(create)) {
    target <- if (grepl("\\.[A-Za-z0-9]+$", basename(p))) dirname(p) else p
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
  }
  p
}

#' Create the full methylQC output directory tree
#'
#' @param dir output directory
#' @param batches whether to also create the per-batch subdirectory
#' @return \code{dir}, invisibly
#' @keywords internal
#' @noRd
mqcmakedirs <- function(dir, batches = FALSE) {
  keys <- c("logs", "matrices", "metadata", "qc")
  if (isTRUE(batches)) keys <- c(keys, "batches")
  for (k in keys) dir.create(mqcpath(dir, k), recursive = TRUE, showWarnings = FALSE)
  invisible(dir)
}

#' Check that a directory looks like a v3 methylQC output directory
#'
#' Gives an explicit error for a v2-era flat layout rather than a confusing
#' "file not found" further down.
#'
#' @param dir directory to check
#' @return \code{TRUE}, invisibly
#' @keywords internal
#' @noRd
mqccheckdir <- function(dir) {
  if (!dir.exists(dir)) stop("directory does not exist: ", dir, call. = FALSE)
  if (file.exists(mqcpath(dir, "betas"))) return(invisible(TRUE))
  if (file.exists(file.path(dir, "betas_all.rds"))) {
    stop("'", dir, "' uses the flat methylQC v2 layout. v3 stores matrices ",
         "under data/matrices/. Re-run qc() to regenerate, or move the files ",
         "manually.", call. = FALSE)
  }
  stop("'", dir, "' does not contain methylQC Stage 1 output ",
       "(expected ", mqcpath(dir, "betas"), "). Run qc() first.", call. = FALSE)
}


## ---------------------------------------------------------------------------
## Atomic writes
## ---------------------------------------------------------------------------

#' Write a file atomically
#'
#' Writes to a temporary file in the same directory and renames on success, so
#' that an interrupted run leaves the previous version intact rather than a
#' truncated file. This is what keeps the single sample_sheet.csv safe when
#' Stage 2 updates it in place.
#'
#' @param x object to write
#' @param path destination path
#' @param writer a function of \code{(x, path)} performing the write
#' @return \code{path}, invisibly
#' @keywords internal
#' @noRd
write_atomic <- function(x, path, writer) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp", Sys.getpid())
  ok <- FALSE
  on.exit(if (!ok && file.exists(tmp)) unlink(tmp), add = TRUE)
  writer(x, tmp)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE))
      stop("could not write ", path, call. = FALSE)
    unlink(tmp)
  }
  ok <- TRUE
  invisible(path)
}

#' @keywords internal
#' @noRd
write_csv_atomic <- function(x, path) {
  write_atomic(x, path, function(o, p) utils::write.csv(o, p, row.names = FALSE))
}

#' @keywords internal
#' @noRd
save_rds_atomic <- function(x, path) {
  write_atomic(x, path, function(o, p) saveRDS(o, p))
}


## ---------------------------------------------------------------------------
## Memory probes  (macOS + Linux only; Windows is not supported)
## ---------------------------------------------------------------------------

#' Query system memory
#'
#' Reads total and available physical memory. On Linux this parses
#' \code{/proc/meminfo} and uses \code{MemAvailable} rather than
#' \code{MemFree}, because \code{MemFree} ignores reclaimable page cache and
#' would badly under-report. On macOS it combines \code{sysctl hw.memsize}
#' with free, inactive and speculative page counts from \code{vm_stat}.
#'
#' @return a list with \code{total} and \code{available} in bytes, plus
#'   \code{source}. Both are \code{NA_real_} when the platform is unsupported
#'   or the probe fails.
#' @export
#' @examples
#' str(sysmem())
sysmem <- function() {
  na <- list(total = NA_real_, available = NA_real_, source = "unsupported")
  os <- Sys.info()[["sysname"]]

  if (identical(os, "Linux")) {
    mi <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) NULL)
    if (is.null(mi)) return(na)
    grab <- function(field) {
      ln <- grep(paste0("^", field, ":"), mi, value = TRUE)
      if (!length(ln)) return(NA_real_)
      as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", ln[1])) * 1024
    }
    return(list(total = grab("MemTotal"), available = grab("MemAvailable"),
                source = "/proc/meminfo"))
  }

  if (identical(os, "Darwin")) {
    num <- function(cmd, args) {
      out <- tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE)),
                      error = function(e) character(0))
      if (!length(out)) return(NA_real_)
      suppressWarnings(as.numeric(out[1]))
    }
    total <- num("sysctl", c("-n", "hw.memsize"))
    psize <- num("sysctl", c("-n", "hw.pagesize"))
    vm <- tryCatch(suppressWarnings(system2("vm_stat", stdout = TRUE, stderr = FALSE)),
                   error = function(e) character(0))
    avail <- NA_real_
    if (length(vm) && !is.na(psize)) {
      pages <- function(label) {
        ln <- grep(label, vm, value = TRUE, fixed = TRUE)
        if (!length(ln)) return(0)
        suppressWarnings(as.numeric(gsub("[^0-9]", "", ln[1])))
      }
      free <- pages("Pages free")
      inact <- pages("Pages inactive")
      spec <- pages("Pages speculative")
      avail <- (free + inact + spec) * psize
    }
    return(list(total = total, available = avail, source = "sysctl/vm_stat"))
  }

  na
}

#' Free space on the filesystem holding a path
#'
#' Used to guard the \code{sdfs_all.rds} write, which is the largest single
#' artefact the pipeline produces. \code{df -Pk} is POSIX and behaves
#' identically on macOS and Linux, so no new dependency is needed.
#'
#' @param path a path on the filesystem to query; need not exist yet, the
#'   nearest existing parent is used.
#' @return free bytes, or \code{NA_real_} if the probe fails.
#' @export
#' @examples
#' diskfree(tempdir())
diskfree <- function(path) {
  while (nzchar(path) && !dir.exists(path)) {
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  if (!dir.exists(path)) return(NA_real_)
  out <- tryCatch(
    suppressWarnings(system2("df", c("-Pk", shQuote(path)),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  if (length(out) < 2L) return(NA_real_)
  f <- strsplit(trimws(out[2]), "[[:space:]]+")[[1]]
  if (length(f) < 4L) return(NA_real_)
  suppressWarnings(as.numeric(f[4])) * 1024
}

#' Resident memory of the current R process
#'
#' @return bytes, or \code{NA_real_} if unavailable.
#' @export
#' @examples
#' procmem()
procmem <- function() {
  os <- Sys.info()[["sysname"]]
  if (identical(os, "Linux")) {
    st <- tryCatch(readLines("/proc/self/status", warn = FALSE), error = function(e) NULL)
    if (is.null(st)) return(NA_real_)
    ln <- grep("^VmRSS:", st, value = TRUE)
    if (!length(ln)) return(NA_real_)
    return(as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", ln[1])) * 1024)
  }
  if (identical(os, "Darwin")) {
    out <- tryCatch(
      suppressWarnings(system2("ps", c("-o", "rss=", "-p", Sys.getpid()),
                               stdout = TRUE, stderr = FALSE)),
      error = function(e) character(0))
    if (!length(out)) return(NA_real_)
    return(suppressWarnings(as.numeric(trimws(out[1]))) * 1024)
  }
  NA_real_
}

#' Format a byte count for humans
#'
#' @param x bytes
#' @param digits decimal places
#' @return character vector
#' @keywords internal
#' @noRd
fmtbytes <- function(x, digits = 1) {
  ifelse(is.na(x), "unknown",
         ifelse(x >= 2^30, sprintf(paste0("%.", digits, "f GiB"), x / 2^30),
                ifelse(x >= 2^20, sprintf(paste0("%.", digits, "f MiB"), x / 2^20),
                       sprintf("%.0f KiB", x / 2^10))))
}

#' Format a duration in seconds for humans
#' @param s seconds
#' @return character
#' @keywords internal
#' @noRd
fmtdur <- function(s) {
  if (is.na(s)) return("unknown")
  if (s < 90) return(sprintf("%.0f s", s))
  if (s < 5400) return(sprintf("%.1f min", s / 60))
  sprintf("%.1f h", s / 3600)
}


## ---------------------------------------------------------------------------
## Alignment guard
## ---------------------------------------------------------------------------

#' Align a named vector to a reference set of probe IDs
#'
#' The v2 code assigned betas positionally (\code{betas[, i] <- b}) while
#' aligning the mask and detection matrices by name, so the two could in
#' principle drift apart within a single sample. This guard makes the
#' alignment explicit for every matrix. \code{identical()} on two character
#' vectors of EPIC length costs under a millisecond because R interns strings,
#' so running it on every sample is free; the reorder only happens when the
#' order actually differs.
#'
#' @param v a named vector
#' @param ids reference probe IDs
#' @param what label used in error messages
#' @param fill value used for IDs absent from \code{v}; when \code{NULL},
#'   missing IDs are an error
#' @return \code{v} reordered to match \code{ids}
#' @keywords internal
#' @noRd
align_to <- function(v, ids, what = "vector", fill = NULL) {
  if (is.null(names(v))) {
    if (length(v) != length(ids))
      stop(what, ": unnamed and length ", length(v), " != ", length(ids),
           call. = FALSE)
    return(stats::setNames(v, ids))
  }
  if (identical(names(v), ids)) return(v)
  idx <- match(ids, names(v))
  if (anyNA(idx)) {
    if (is.null(fill))
      stop(what, ": ", sum(is.na(idx)), " of ", length(ids),
           " probe IDs are missing (e.g. '",
           ids[which(is.na(idx))[1]], "'). Probe sets differ.", call. = FALSE)
    out <- stats::setNames(rep(fill, length(ids)), ids)
    ok <- !is.na(idx)
    out[ok] <- unname(v[idx[ok]])
    return(out)
  }
  out <- stats::setNames(unname(v[idx]), ids)
  stopifnot(identical(names(out), ids))
  out
}

## ---------------------------------------------------------------------------
## Row statistics
##
## matrixStats is used when available because it avoids materialising an
## intermediate copy, which matters at 100,000 x 1500. The base fallbacks keep
## the package usable (and testable) without it.
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.rowVars <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE))
    return(matrixStats::rowVars(x))
  n <- ncol(x)
  if (n < 2L) return(rep(NA_real_, nrow(x)))
  m <- rowMeans(x, na.rm = TRUE)
  (rowSums(x * x, na.rm = TRUE) - n * m * m) / (n - 1)
}

#' @keywords internal
#' @noRd
.rowMedians <- function(x, na.rm = TRUE) {
  if (requireNamespace("matrixStats", quietly = TRUE))
    return(matrixStats::rowMedians(x, na.rm = na.rm))
  apply(x, 1L, stats::median, na.rm = na.rm)
}

#' Assert that two matrices are conformable by label, reordering if needed
#'
#' v2 checked only \code{identical(dim(a), dim(b))}, which passes for two
#' matrices of the same shape but different row or column order. This compares
#' the labels.
#'
#' @param a matrix to align
#' @param ref reference matrix
#' @param what label used in messages
#' @return \code{a}, reordered to match \code{ref} if necessary
#' @keywords internal
#' @noRd
conform_to <- function(a, ref, what = "matrix") {
  if (is.null(a)) return(NULL)
  if (!is.matrix(a)) stop(what, " is not a matrix", call. = FALSE)
  rn <- rownames(ref); cn <- colnames(ref)
  if (is.null(rownames(a)) || is.null(colnames(a))) {
    if (!identical(dim(a), dim(ref)))
      stop(what, ": unlabelled and dimensions differ from the reference",
           call. = FALSE)
    dimnames(a) <- dimnames(ref)
    return(a)
  }
  if (identical(rownames(a), rn) && identical(colnames(a), cn)) return(a)
  if (!setequal(rownames(a), rn))
    stop(what, ": row (probe) sets differ from the reference matrix",
         call. = FALSE)
  if (!setequal(colnames(a), cn))
    stop(what, ": column (sample) sets differ from the reference matrix",
         call. = FALSE)
  a[match(rn, rownames(a)), match(cn, colnames(a)), drop = FALSE]
}
