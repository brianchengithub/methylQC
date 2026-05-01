# methylQC

**QC pipeline for Illumina DNA methylation arrays**

`methylQC` wraps [SeSAMe](https://bioconductor.org/packages/sesame/), [EpiDISH](https://bioconductor.org/packages/EpiDISH/), and standard Bioconductor tools into a two-stage pipeline for Illumina Infinium arrays (EPIC, EPICv2, 450K).

The pipeline outputs **raw (unmasked) matrices** alongside quality masks and detection p-values, giving you full control over filtering. A user-facing `apply_mask()` function handles masking, exclusions, and optional k-NN imputation.

## Features

- **Streaming preprocessing** via SeSAMe's `openSesame` (memory-efficient)
- **Raw matrices + separate masks**: unmasked betas, quality mask matrix, and detection p-value matrix
- **Configurable column names** for sample ID, donor, sex, age, batch, and cell type
- **Sex check** with optimized threshold clustering and regression-based confidence bands
- **Epigenetic age prediction** via hardcoded Horvath (2013) clock (353 CpGs)
- **SNP identity verification** via MDS on rs probes
- **Multi-page QC report PDF**: detection rate, MDS, intensity, probe failure, beta density, sex check, age check, PCA
- **Cell-type deconvolution** via EpiDISH for blood tissues
- **`apply_mask()`**: flexible filtering with optional chunked k-NN imputation

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("sesame", "sesameData", "EpiDISH", "impute"))

devtools::install_github("your-username/methylQC")
```

## Quick start

```r
library(methylQC)

# Configure column names to match your sample sheets
methylQC_set(
  donor_col = "SubjectID",
  reported_sex_col = "gender"
)

# Run the pipeline
pipeline(
  in_dir  = "/path/to/idats",
  out_dir = "/path/to/output"
)

# Apply masks and filter for downstream analyses
betas  <- readRDS("betas_all.rds")
mask   <- readRDS("mask_all.rds")
detP   <- readRDS("detP_all.rds")
excl_p <- read.csv("exclude_probes.csv")
ss     <- read.csv("sample_sheet.csv")

betas_clean <- apply_mask(
  betas,
  mask = mask,
  detP = detP, detP_thresh = 0.05,
  exclude_probes = excl_p$probe_id,
  exclude_samples = ss$sample_id[ss$flagged],
  probe_types = "cg",
  impute = TRUE
)
```

## Pipeline outputs

| File | Description |
|------|-------------|
| `betas_all.rds` | Unmasked beta matrix (no NAs from masking) |
| `mvals_all.rds` | Unmasked M-value matrix |
| `mask_all.rds` | Logical quality mask (TRUE = failed) |
| `detP_all.rds` | Detection p-values (pOOBAH) |
| `sample_sheet.csv` | All samples with QC metrics, flags, and cell proportions |
| `exclude_probes.csv` | Probes to remove (sex chrom, SNP, low call rate) |
| `qc_plots.pdf` | Multi-page QC diagnostic report |
| `snp_betas.rds` | SNP probe betas for identity checks |

See [METHODS.md](METHODS.md) for detailed technical documentation.

## `apply_mask()` — flexible filtering

Each argument is optional. Apply only what you need:

```r
# Just quality mask
betas_masked <- apply_mask(betas, mask = mask)

# Strict detection threshold + imputation
betas_strict <- apply_mask(betas, detP = detP, detP_thresh = 0.01,
                           impute = TRUE, knn_k = 5)

# Full filtering, no imputation (handle NAs yourself)
betas_filtered <- apply_mask(betas, mask = mask, detP = detP,
                             exclude_probes = excl_p$probe_id,
                             exclude_samples = ss$sample_id[ss$flagged])
```

## Design philosophy

- **Raw matrices + masks**: You choose what to trust. Different analyses may need different thresholds.
- **Flag, don't filter**: The pipeline flags; you exclude. `sample_sheet.csv` has all the flags.
- **Optional imputation**: `apply_mask(..., impute = TRUE)` runs chunked k-NN. Or skip it entirely.
- **No batch correction**: Analysis-specific. PCA plots help assess batch effects.

## Citation

- **SeSAMe**: Zhou W, et al. *Nucleic Acids Research*. 2018;46(20):e123.
- **Horvath clock**: Horvath S. *Genome Biology*. 2013;14(10):R115.
- **EpiDISH**: Teschendorff AE, et al. *BMC Bioinformatics*. 2017;18(1):105.
