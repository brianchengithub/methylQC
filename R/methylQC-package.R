###############################################################################
# methylQC-package.R — Package-level documentation
#
# This file provides the top-level help page for the methylQC package,
# accessible via ?methylQC or help(package = "methylQC"). It is parsed
# by roxygen2 to generate the package .Rd file.
###############################################################################
#' methylQC: Illumina Methylation Array QC and Preprocessing
#'
#' A reproducible pipeline for Illumina methylation array data using
#' SeSAMe for preprocessing and QC, EpiDISH for cell-type deconvolution,
#' and EpiDISH for cell-type deconvolution of blood tissues.
#'
#' @section Quick start:
#' \preformatted{
#' library(methylQC)
#'
#' # Stage 1: QC and preprocessing (auto-detects platform)
#' stage1 <- qc(
#'   in_dir  = "data/idats",
#'   out_dir = "results"
#' )
#'
#' # Review exclusion CSVs in `results/`, edit if needed, then:
#'
#' # Stage 2: prep — apply exclusions, impute, deconvolve
#' stage2 <- prep(stage1_dir = "results")
#' }
#'
#' @section Configuration:
#' All tunable parameters are stored as package options with the prefix
#' \code{methylQC.}. See \code{\link{methylQC_defaults}} for the full list
#' and \code{\link{methylQC_options}} / \code{\link{methylQC_set}} to
#' inspect or override them.
#'
#' @keywords internal
"_PACKAGE"