###############################################################################
# methylQC-package.R — Package-level documentation
#
# This file provides the top-level help page for the methylQC package,
# accessible via ?methylQC or help(package = "methylQC"). It is parsed
# by roxygen2 to generate the package .Rd file.
###############################################################################
#' methylQC: Illumina Methylation Array QC and Preprocessing
#'
#' A reproducible two-stage pipeline for Illumina DNA methylation array
#' data. Stage 1 streams IDATs through SeSAMe (openSesame, prep code
#' "QCDPB"), producing unmasked noob-corrected beta and M-value matrices
#' alongside separate quality-mask and detection p-value matrices, and
#' computes per-sample QC metrics. Stage 2 flags samples and probes,
#' generates a multi-page QC diagnostic report, and runs EpiDISH
#' cell-type deconvolution as a purity check on blood-derived tissues.
#' Optional EPICv2 replicate de-duplication and probe-ID harmonization,
#' and similarity-based k-NN imputation, are also provided.
#'
#' @section Quick start:
#' \preformatted{
#' library(methylQC)
#' check_dependencies()         # verify the environment once
#'
#' # Stage 1: QC and preprocessing (auto-detects platform)
#' stage1 <- qc(
#'   in_dir  = "data/idats",
#'   out_dir = "results"
#' )
#'
#' # Review exclusion CSVs in `results/`, edit if needed, then:
#'
#' # Stage 2: prep — flag samples/probes, QC report, deconvolution
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