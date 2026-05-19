# methylQC

**QC pipeline for Illumina DNA methylation arrays**

`methylQC` wraps [SeSAMe](https://bioconductor.org/packages/sesame/), [EpiDISH](https://bioconductor.org/packages/EpiDISH/), and standard Bioconductor tools into a two-stage pipeline for Illumina Infinium arrays (EPIC, EPICv2, 450K).

The pipeline outputs **unmasked, noob-corrected matrices** alongside quality masks and detection p-values, giving you full control over filtering. A user-facing `apply_mask()` function handles masking, exclusions, and optional similarity-based k-NN imputation.

## Features

- **Streaming preprocessing** via SeSAMe's `openSesame` (memory-efficient), prep code `QCDPB`
- **Unmasked, noob-corrected matrices + separate masks**: beta and M-value matrices, a quality-mask matrix, and a detection p-value matrix
- **Configurable column names** for sample ID, donor, sex, age, batch, and cell type
- **Two distinct QC layers**: sample *quality* vs. sample *identity / swaps* (see below)
- **Epigenetic age prediction** via hard-coded Horvath (2013) clock (353 CpGs), with zero-shot reference imputation for clock CpGs absent from a platform
- **Cell-type deconvolution** via EpiDISH (12 immune cell-type panel) as a purity check
- **Optional EPICv2 de-duplication + probe-ID harmonization** to EPIC or HM450
- **Multi-page QC report PDF**: detection rate, MDS, intensity, probe failure, beta density, sex check, age check, PCA
- **`apply_mask()`**: flexible filtering with optional similarity-based k-NN imputation

## The two QC layers: quality vs. identity

methylQC separates QC into two conceptually distinct sets of checks. **A sample can pass one and fail the other**, and conflating them hides real problems.

**Sample quality** — *Is the data technically sound?* This is a conservative screen for abnormally concerning samples. It relies primarily on:

- **MDS outliers** — samples that sit far from the bulk in multidimensional scaling space.
- **Total intensity outliers** — samples with abnormally low overall signal.
- **Detection p-value call rates** — the fraction of probes passing pOOBAH, examined per sample.

The quality layer flags samples; it does not silently delete data. Detection p-values are provided as a full matrix so each analysis can choose its own threshold.

**Sample identity / swaps** — *Is this sample who it claims to be?* A swapped or mislabeled sample can produce perfectly high-quality data and therefore would **never** be caught by the quality layer. Identity checks are independent:

- **Sex check** — predicted sex (from sex chromosome intensities) vs. reported sex.
- **Age check** — Horvath epigenetic age vs. reported age.
- **SNP genotyping** — concordance of the 59 rs probes (e.g., MDS clustering of genotypes; concordance against external genotype data where available).

Because the failure modes are different, the two layers are reported separately throughout the pipeline outputs.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("sesame", "sesameData", "EpiDISH"))

devtools::install_github("brianchengithub/methylQC")
```

After installing, verify your environment once:

```r
library(methylQC)
check_dependencies()
```

`check_dependencies()` reports the R version, the versions of required
packages, the state of the SeSAMe data cache, and whether the default
EpiDISH reference panel is available. It **does not install anything** —
for anything missing or out of date it prints the exact command to run.
If you have never initialized the SeSAMe cache, run `sesameData::sesameDataCacheAll()` once.

## Quick start

```r
library(methylQC)

# Configure column names to match your sample sheets
methylQC_set(
  donor_col = "SubjectID",
  reported_sex_col = "gender"
)

# Run the pipeline. in_dir is the single top-level directory; IDATs and
# sample sheets are discovered recursively within it.
pipeline(
  in_dir  = "/path/to/idats/CellType",
  out_dir = "/path/to/output/CellType"
)

# Apply masks and filter downstream
betas  <- readRDS("output/betas_all.rds")
mask   <- readRDS("output/mask_all.rds")
detP   <- readRDS("output/detP_all.rds")
excl_p <- read.csv("output/exclude_probes.csv")
ss     <- read.csv("output/sample_sheet.csv")

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

## Preprocessing: the QCDPB prep code

Stage 1 runs `openSesame` with prep code **`QCDPB`**:

| Code | Step |
|------|------|
| `Q` | Quality masking (design issues, cross-hybridization) |
| `C` | Channel (color) inference |
| `D` | Dye-bias correction |
| `P` | pOOBAH detection p-values |
| `B` | noob background correction |

The order matters. **pOOBAH (`P`) runs before noob (`B`)** because noob modifies the out-of-band signal that pOOBAH uses as its null distribution. Consequently:

- **Detection p-values** (`detP_all.rds`) are computed on **non-noob** signal — the statistically correct choice.
- **Beta and M-value matrices** (`betas_all.rds`, `mvals_all.rds`) **are noob-corrected**.
- The **MDS** plot uses the noob-corrected betas.
- **Sex chromosome intensities** and **mean intensity** are deliberately taken from the **non-noob** signal: noob only subtracts background and cannot recover weak signal, so it offers no benefit for chromosome-presence calls and would only push weak-signal samples toward zero.
- The **Horvath clock** is computed on the **unmasked, noob-corrected** betas.

## Pipeline outputs

| File | Description |
|------|-------------|
| `betas_all.rds` | Unmasked, noob-corrected beta matrix (probes x samples) |
| `mvals_all.rds` | Unmasked, noob-corrected M-value matrix |
| `mask_all.rds` | Logical quality mask (TRUE = failed) |
| `detP_all.rds` | Detection p-values (pOOBAH, non-noob) |
| `sdfs_all.rds` | List of final SigDF objects (saved by default) |
| `sample_sheet.csv` | All samples with QC metrics, flags, and cell proportions |
| `exclude_probes.csv` | Probes flagged (sex chrom, SNP, low call rate) |
| `exclude_samples.csv` | Samples flagged (low detection, low intensity) |
| `qc_plots.pdf` | Multi-page QC diagnostic report |
| `snp_betas.rds` | SNP probe betas for identity checks |
| `epicv2_dedup_log.csv` | Per-CpG kept/dropped log (only if EPICv2 harmonization is on) |

See [METHODS.md](METHODS.md) for detailed technical documentation.

## EpiDISH cell-type deconvolution

EpiDISH (RPC) is run as a **purity check** on blood-derived tissues — whole blood, PBMC, buffy coat, and isolated/sorted blood cell populations. For a cleanly sorted population a single cell type should dominate; high immune signal in a putatively non-immune sample indicates contamination.

EpiDISH (RPC) constrains estimates to sum to 1, so it is only meaningful for immune/blood tissue — running an immune-only reference on a non-immune tissue would force a meaningless simplex. The pipeline therefore restricts deconvolution to the blood/immune allowlist (`epidish_cell_types`).

### Reference panels

The reference is set via the `epidish_reference` option and may be a single name or a character vector (each panel is run and returned as a named list).

| Reference (`EpiDISH`) | Cell types | Notes |
|------------------------|-----------|-------|
| `cent12CT.m` *(default)* | 12 leukocyte subtypes — naive/memory CD4T, naive/memory CD8T, naive/memory B, Treg, NK, monocytes, neutrophils, eosinophils, basophils | Most granular; distinguishes naive vs. memory subsets. |
| `centDHSbloodDMC.m` | 7 types — B, NK, CD4T, CD8T, monocytes, neutrophils, eosinophils | The classic blood reference. |
| `centEpiFibIC.m` | Epithelial, fibroblast, total immune | For solid-tissue epithelial/fibroblast/immune fractions (not blood subtyping). |
| `centBloodSub.m` | Blood sub-types | Available in recent EpiDISH versions. |

```r
# Use multiple references at once
methylQC_set(epidish_reference = c("cent12CT.m", "centDHSbloodDMC.m"))
```

## EPICv2 de-duplication and probe-ID harmonization

EPICv2 targets many CpGs with multiple replicate probes (IDs sharing a `cg`/`ch` prefix but differing by suffix). For a cohort analysis every sample must use the **same physical probe** for a given CpG, or probe-design differences become technical noise.

This behavior is **off by default**. Enable it with a single toggle:

```r
methylQC_set(epicv2_harmonize = TRUE,
             epicv2_target    = "EPIC")   # or "HM450"
```

When enabled (and the platform is EPICv2), Stage 1:

1. **De-duplicates** replicate `cg`/`ch` probes cohort-consistently: for each replicated CpG, the probe with the **fewest cross-sample detection failures** is kept (ties → first probe by ID order). `rs` and control probes are left untouched. This is a *cohort-level* decision, unlike SeSAMe's per-sample `collapseToPfx(method = "minPvalue")`.
2. **Harmonizes** the surviving probe IDs to the target legacy platform (`EPIC` by default, or `HM450`) via SeSAMe's `mLiftOver`, which uses conversion mappings shipped in the `sesameData` cache.

A per-CpG kept/dropped log is written to `epicv2_dedup_log.csv` (no console dump).

### Standalone function

`harmonize_epicv2()` runs the same de-dup + harmonization on matrices you already have, independent of the pipeline:

```r
harmonized <- harmonize_epicv2(
  mat_path        = "results/betas_all.rds",   # beta OR M-value matrix
  detP_path       = "results/detP_all.rds",    # detection p-value matrix
  target_platform = "EPIC",                    # or "HM450"
  dedup_log_path  = "results/epicv2_dedup_log.csv"
)
```

**Input matrix format.** Both inputs are matrices with **probes in rows, samples in columns**, and EPICv2 probe IDs as the row names. The detection p-value matrix is required — the de-duplication criterion is the cross-sample detection failure rate, which cannot be derived from betas alone. Both matrices are exactly what Stage 1 writes (`betas_all.rds` / `mvals_all.rds` and `detP_all.rds`). To generate them yourself from SeSAMe:

```r
library(sesame)
sdfs  <- openSesame(idat_dir, prep = "QCDP", func = NULL)   # list of SigDFs
betas <- do.call(cbind, lapply(sdfs, getBetas, mask = FALSE))
detP  <- do.call(cbind, lapply(sdfs, pOOBAH, return.pval = TRUE))
```

## `apply_mask()` — flexible filtering

Each argument is optional. Apply only what you need:

```r
# Just quality mask
betas_masked <- apply_mask(betas, mask = mask)

# Strict detection threshold + imputation
betas_strict <- apply_mask(betas, detP = detP, detP_thresh = 0.01,
                           impute = TRUE)

# Full filtering, no imputation (handle NAs yourself)
betas_filtered <- apply_mask(betas, mask = mask, detP = detP,
                             exclude_probes = excl_p$probe_id,
                             exclude_samples = ss$sample_id[ss$flagged])
```

### Imputation

`apply_mask(..., impute = TRUE)` runs **similarity-based k-NN imputation**:

1. Sample-to-sample Euclidean distances are computed on the top `knn_var_probes` (default 30,000) most-variable probes. Restricting to variable probes is what makes "nearest" meaningful — most of the array is near-constant and would otherwise dominate the distance.
2. Each missing value is filled with the **inverse-distance-weighted mean** of the `knn_k` (default 50) nearest-neighbour samples' values at that probe.

There is no constant fallback: a value missing in a sample and in all of its usable neighbours is left as `NA`.

## Design philosophy

- **Unmasked matrices + masks**: You choose what to trust. Different analyses may need different thresholds.
- **Flag, don't filter**: The pipeline flags; you exclude. `sample_sheet.csv` carries all flags.
- **Quality vs. identity are separate**: a clean swapped sample fails identity checks but passes quality checks — both are reported.
- **Optional imputation**: `apply_mask(..., impute = TRUE)` runs similarity k-NN, or skip it entirely.
- **No silent package installs**: `check_dependencies()` reports problems and the exact fix; nothing is installed behind your back.

## Citation

- **SeSAMe**: Zhou W, et al. *Nucleic Acids Research*. 2018;46(20):e123.
- **mLiftOver**: Zhou W, et al. *Bioinformatics*. 2024;40(7):btae423.
- **Horvath clock**: Horvath S. *Genome Biology*. 2013;14(10):R115.
- **EpiDISH**: Teschendorff AE, et al. *BMC Bioinformatics*. 2017;18(1):105.
