###############################################################################
# methylQC-package.R - Package-level documentation
#
# This file provides the top-level help page for the methylQC package,
# accessible via ?methylQC or help(package = "methylQC"). It is parsed
# by roxygen2 to generate the package .Rd file.
###############################################################################
#' methylQC: Illumina Methylation Array QC and Preprocessing
#'
#' A reproducible pipeline for Illumina methylation array data using
#' SeSAMe for preprocessing, EpiDISH for cell-type deconvolution
#' (blood tissues), and inline QC diagnostics. Samples and probes are
#' never removed automatically: the package produces flags and a
#' diagnostic report, and the user applies exclusions via
#' \code{\link{cleanmat}}.
#'
#' @section Quick start:
#' \preformatted{
#' library(methylQC)
#'
#' # Stage 1: QC and preprocessing (auto-detects platform)
#' stage1 <- qc(
#'   indir  = "data/idats",
#'   outdir = "results"
#' )
#'
#' # Review the QC PDF in results/, then run Stage 2:
#'
#' # Stage 2: per-cell-type flagging, EpiDISH, consolidated sheet
#' stage2 <- prep(dir = "results")
#'
#' # To act on flags, point cleanmat() at the QC artifacts directly:
#' detP  <- readRDS("results/detP_all.rds")
#' mask  <- readRDS("results/mask_all.rds")
#' betas <- readRDS("results/betas_all.rds")
#' flagged <- flagsamples(detP, callrate = 0.95,
#'                        csv = "results/flagged_samples.csv")
#' betas_clean <- cleanmat(
#'   betas, mask = mask, detP = detP,
#'   dropprobes  = "results/failed_probes.csv",   # CSV path or character vec
#'   dropsamples = "results/flagged_samples.csv", # CSV path or character vec
#'   probes      = c("cg", "ch"),
#'   platform    = "EPIC",
#'   impute      = TRUE)
#' }
#'
#' @section Configuration:
#' All tunable parameters are stored as package options with the prefix
#' \code{methylQC.}. See \code{\link{mqcdefaults}} for the full list
#' and \code{\link{mqcopts}} / \code{\link{mqcset}} to inspect or
#' override them.
#'
#' @keywords internal
"_PACKAGE"
