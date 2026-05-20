# methylQC

**QC and preprocessing pipeline for Illumina DNA methylation arrays.**

`methylQC` wraps [SeSAMe](https://bioconductor.org/packages/sesame/),
[EpiDISH](https://bioconductor.org/packages/EpiDISH/), and standard
Bioconductor tools into a two-stage pipeline for Illumina Infinium
arrays (EPIC, EPICv2, 450K).

The pipeline outputs **raw (unmasked) matrices** alongside a quality
mask and detection p-values; nothing is excluded automatically.
`cleanmat()` is the single function that turns those matrices into a
filtered analysis matrix, on the user's terms.

## Features

- **Streaming preprocessing** through SeSAMe `openSesame` (one IDAT
  pair at a time; bounded memory).
- **Raw matrices + separate masks** — unmasked betas, quality mask
  matrix, and detection p-value matrix are saved side by side.
- **Configurable column names** for sample ID, donor, sex, age,
  batch, and cell type (with alias lists).
- **Sex check** with adaptive within-cluster regression bands.
- **Epigenetic age prediction** via the hardcoded Horvath (2013)
  clock (353 CpGs).
- **SNP identity verification** via MDS on rs probes.
- **Nine-page QC report PDF**: detection rate, probe-failure tail
  histogram, MDS, intensity, sample beta density, sex check, age
  check, scree, and a 2×3 panel of PC-vs-associated-variable plots.
- **Cell-type deconvolution** (EpiDISH RPC) for blood-derived
  tissues.
- **`cleanmat()`** — single primitive for applying QC decisions, with optional
  chunked k-NN imputation.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("sesame", "sesameData", "EpiDISH", "impute"))

devtools::install_github("brianchengithub/methylQC")
```

## Quick start

```r
library(methylQC)

# Tell methylQC how your sample sheets are labelled
mqcset(
  donorcol = "SubjectID",
  sexcol   = "gender"
)

# Stage 1 + Stage 2 in one call
pipeline(
  indir  = "/path/to/idats",
  outdir = "/path/to/output"
)
```

Stage 1 writes raw matrices and a sample sheet with per-sample QC
metrics. Stage 2 writes the QC report PDF, optional EpiDISH cell
proportions, and a consolidated sample sheet with flag columns.

To apply filters before downstream analysis:

```r
betas <- readRDS("output/betas_all.rds")
mask  <- readRDS("output/mask_all.rds")
detP  <- readRDS("output/detP_all.rds")
ss    <- read.csv("output/sample_sheet.csv")

# Standard autosomal-CpG EWAS matrix; feed cleanmat the QC CSVs directly
betas_ewas <- cleanmat(
  betas, mask = mask, detP = detP, pthresh = 0.05,
  dropprobes  = "output/failed_probes.csv",
  dropsamples = "output/flagged_samples.csv",
  probes      = "cg",
  platform    = "EPIC")
```

If you want an advisory list of samples below a custom call-rate
threshold:

```r
flagsamples(detP, callrate = 0.95,
            csv = "output/flagged_samples.csv")
# Then cleanmat() will consume that file directly:
betas_clean <- cleanmat(
  betas, mask = mask, detP = detP,
  dropsamples = "output/flagged_samples.csv",
  dropprobes  = "output/failed_probes.csv",
  probes      = c("cg", "ch"),
  platform    = "EPIC")
```

## Pipeline outputs

| File                         | Description                                                              |
|------------------------------|--------------------------------------------------------------------------|
| `betas_all.rds`              | Unmasked beta matrix (no NAs from masking).                              |
| `mvals_all.rds`              | Unmasked M-value matrix.                                                 |
| `mask_all.rds`               | Logical mask (`TRUE` = failed quality mask or pOOBAH).                   |
| `detP_all.rds`               | Detection p-values (pOOBAH).                                             |
| `snp_betas.rds`              | rs-probe betas (samples × rs probes) for `snpcheck()`.                   |
| `sample_sheet.csv`           | Sample sheet + QC metrics + flags + (optional) cell proportions.         |
| `metadata.rds`               | Run metadata (platform, sample count, timestamp, package version).       |
| `qc_plots.pdf`               | Nine-page QC diagnostic report (Stage 2).                                |
| `pc_scores.csv`              | All PC scores from the page-9 PCA (Stage 2).                             |
| `failed_probes.csv`          | Probes failing in ≥ `failmin` fraction of samples (Stage 2).             |
| `cell_proportions.csv`       | EpiDISH proportions, blood tissues only (Stage 2).                       |
| `pipeline_diagnostics.log`   | Per-step log (timestamps, parameter values, summary stats).              |

There is **no** `exclude_samples.csv`, **no** `exclude_probes.csv`,
**no** `probe_call_rates.csv` in v2.0.0. If you want a sample call-rate
flag list, call `flagsamples()`. Probe-level masking lives in
`mask_all.rds` and `detP_all.rds`.

See [METHODS.md](METHODS.md) for the full technical reference: exact
thresholds, the QCDP + B preprocessing order, the `frac_dt`
vs `colMeans(detP <= 0.05)` distinction, the QC-plot algorithms, the
`cleanmat()` order of operations, and design rationale.

## `cleanmat()` in one screen

```r
# Just quality mask
betas_masked <- cleanmat(betas, mask = mask, probes = NULL)

# Strict detection + k-NN imputation; consume the QC CSVs directly
betas_strict <- cleanmat(
  betas, mask = mask, detP = detP, pthresh = 0.01,
  dropprobes  = "output/failed_probes.csv",
  dropsamples = "output/flagged_samples.csv",
  probes      = "cg", platform = "EPIC",
  impute      = TRUE, knnk = 10)

# Sex chromosome probes only
betas_sex <- cleanmat(
  betas, mask = mask, detP = detP,
  probes = "sex", platform = "EPIC")

# SNP probes only
betas_snp <- cleanmat(
  betas, probes = "snp")
```

`probes` accepts any subset of `c("cg", "ch", "sex", "snp", "other")`,
or `NULL` to keep all probes. `"cg"` always means autosomal CpG; sex
chromosome cg probes live in `"sex"`.

## Design

- **Raw matrices + masks.** Signal and QC are decoupled. Different
  downstream analyses need different masks.
- **Flag, don't filter.** The pipeline produces flags and a report.
  Filtering is the user's responsibility, via `cleanmat()`.
- **MDS on all probes, PCA on top-variable.** The page-3 MDS uses
  every complete-case cg/ch non-sex probe; the page-9 PCA uses the
  top `ntop` (default 100,000) by variance.
- **No batch correction.** Analysis-specific. The page-9 panel
  surfaces sample-sheet variables that associate with the leading
  PCs.

## Citation

- **SeSAMe**: Zhou W, Triche TJ, Laird PW, Shen H.
  *Nucleic Acids Research* 2018;46(20):e123.
- **NOOB**: Triche TJ Jr, Weisenberger DJ, Van Den Berg D, Laird PW,
  Siegmund KD. *Nucleic Acids Research* 2013;41(7):e90.
- **Horvath clock**: Horvath S. *Genome Biology* 2013;14(10):R115.
- **EpiDISH**: Teschendorff AE, Breeze CE, Zheng SC, Beck S.
  *BMC Bioinformatics* 2017;18(1):105.
