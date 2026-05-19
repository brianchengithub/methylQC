# methylQC: Methods and Technical Documentation

**Version 1.2.0**

## 1. Overview

methylQC is a two-stage quality control pipeline for Illumina Infinium DNA methylation arrays (EPIC, EPICv2, 450K). It produces **unmasked, noob-corrected** beta and M-value matrices alongside separate quality mask and detection p-value matrices, giving users full control over which probe values to trust.

The pipeline **flags** samples and probes but does not automatically filter them. A user-facing `apply_mask()` function applies any combination of quality masking, detection p-value thresholds, sex chromosome exclusion, probe/sample exclusions, probe type filtering, and optional similarity-based k-NN imputation.

### 1.1 Two QC layers: quality vs. identity

methylQC deliberately separates QC into two conceptually distinct sets of checks, because they detect different failure modes and a sample can pass one while failing the other.

**Sample quality** asks *is the data technically sound?* It is a conservative screen for abnormally concerning samples and relies primarily on **MDS outliers**, **total intensity outliers**, and **detection p-value call rates**. The quality layer flags samples but does not delete data; detection p-values are exported as a full matrix so each downstream analysis sets its own threshold.

**Sample identity / swaps** asks *is this sample who it claims to be?* It relies on the **sex check**, the **age check**, and **SNP genotyping**. These are independent of quality: a swapped or mislabeled sample frequently generates perfectly high-quality data and therefore would never be flagged by the quality layer. For this reason identity checks (`sex_mismatch`, `age_outlier`, SNP concordance) never contribute to the `flagged` quality column — they are reported in their own columns.

## 2. Stage 1: Preprocessing (`qc()`)

### 2.1 IDAT discovery and sample sheet parsing

The user supplies a single top-level directory. `discover_idats()` confines all discovery to that tree:

1. Every IDAT pair (`_Grn.idat` / `_Red.idat`) is found recursively.
2. Every sample sheet matching the configurable pattern is found recursively and read (tab/comma auto-detected). Multiple sheets are **concatenated**, taking the union of their columns and filling absent columns with `NA`.
3. Sheet rows are reconciled against the IDATs actually on disk. IDATs present on disk but absent from every sheet have a minimal row **synthesized** and appended (not written back to disk). Rows whose IDATs are missing on disk trigger a **warning** and are dropped. A mismatch is always a warning, never an error.
4. Duplicate `sample_id` values (one sample appearing in more than one sheet row) are resolved by keeping the row with the **most non-missing values** (ties → first occurrence); a warning naming the `sample_id` and source sheet file(s) is emitted.

The array platform is auto-detected from IDAT headers. Sentrix IDs are decomposed into Chip, Row, and Col.

### 2.2 Signal processing (openSesame)

Each sample is processed through SeSAMe's `openSesame(func = NULL, prep = "QCDPB")`. The `func = NULL` argument returns the processed SigDF. The `prep = "QCDPB"` pipeline applies five steps:

1. **Q — Quality masking**: Probes with design issues are flagged in the SigDF's `mask` column (multi-mapping, SNP-overlapping CpG/SBE sites, cross-hybridizing sequences). This is a fixed, platform-level mask independent of sample signal.
2. **C — Channel inference**: Infinium-I probe color channel is inferred from the data.
3. **D — Dye-bias correction**: Non-linear Cy3/Cy5 channel efficiency correction.
4. **P — pOOBAH detection p-values**: Each probe's signal is compared to the out-of-band hybridization distribution. Probes with p > 0.05 are additionally flagged in the `mask` column. Sample-specific.
5. **B — noob background correction**: Normal-exponential out-of-band background subtraction.

**Order matters: `P` runs before `B`.** noob modifies the out-of-band signal that pOOBAH uses as its null distribution, so detection p-values must be computed *before* noob. To honour this, the pipeline derives two SigDF states per sample:

- A **pre-noob** SigDF (prep code up to and including `P`, i.e. `QCDP`): the **quality mask** and **detection p-values** are extracted from here.
- The **full QCDPB** SigDF: the **beta and M-value** matrices are extracted from here, via `getBetas(sdf, mask = FALSE)` — betas come from noob-corrected intensities but masked probes are **not** set to NA.

Consequently `detP_all.rds` is non-noob (correct), while `betas_all.rds` / `mvals_all.rds` are noob-corrected. The MDS plot uses the noob-corrected betas.

### 2.3 Output matrices

All matrices are probes (rows) × samples (columns) in `.rds` format. Row names are probe IDs; column names are sample IDs.

| File | Type | Description |
|------|------|-------------|
| `betas_all.rds` | numeric | Unmasked, noob-corrected betas. No NAs from masking. |
| `mvals_all.rds` | numeric | `log2((beta + 1e-6) / (1 - beta + 1e-6))`, noob-corrected. |
| `mask_all.rds` | logical | TRUE = failed quality mask or detection |
| `detP_all.rds` | numeric | pOOBAH detection p-values (non-noob) |
| `sdfs_all.rds` | list | Final QCDPB SigDF objects (saved by default) |
| `snp_betas.rds` | numeric | **samples × probes** (transposed). Continuous rs probe betas. |

### 2.4 Probe types

Probes are identified by ID prefix. Counts vary by platform:

| Prefix | Type | EPIC | EPICv2 | HM450 |
|--------|------|------|--------|-------|
| `cg` | CpG methylation | ~846,000 | ~930,000 | ~482,000 |
| `ch` | Non-CpG methylation | ~2,800 | ~2,900 | ~3,100 |
| `rs` | SNP genotyping | 59 | 59 | 65 |
| other | Non-standard | ~635 | varies | varies |

The rs probe count is auto-detected by prefix (`grepl("^rs", ...)`), not hardcoded.

The "other" probes are defined by exclusion: any probe whose ID does not start with `cg`, `ch`, or `rs`. On EPIC, these include `nv` (negative verification) probes and other platform-specific probes. Users can inspect them:

```r
excl_p <- read.csv("results/exclude_probes.csv")
other <- excl_p[grepl("non_cg_ch_probe", excl_p$reason), ]
table(substr(other$probe_id, 1, 2))
```

Sex chromosome probes (~19,642 on EPIC) are a **subset of the cg probes** — all sex probes have the `cg` prefix. They are handled by a dedicated `exclude_sex` parameter in `apply_mask()`, independent of `probe_types`.

### 2.5 Per-sample QC metrics

Computed from each SigDF during streaming:

**Detection rate (`frac_dt`)**: Fraction of probes passing pOOBAH.

**Mean intensity, probe decomposition**: `n_masked_apriori`, `n_failed_detection`, `n_usable`, `frac_usable`.

**Horvath epigenetic age (`horvath_age`)**: 353 CpG coefficients hard-coded from Horvath (2013). Computed on the **unmasked, noob-corrected** betas (`getBetas(sdf, mask = FALSE)` on the QCDPB SigDF). Per-CpG value precedence: (1) a clock CpG present on the platform contributes its actual beta — a present-but-failed probe still contributes its real value, since the age check is a coarse identity screen and is not re-imputed; (2) a clock CpG **absent from the platform manifest** contributes a hard-coded zero-shot reference median (`.horvath_zeroshot_betas`, 19 CpGs, derived from a large HM450 blood dataset). There is no 0.5 fallback. If fewer than 200 of the 353 clock CpGs can be resolved, the prediction is `NA`.

**Sex chromosome intensities (`sex_chrX_intensity`, `sex_chrY_intensity`)**: Median total intensity (MG + MR + UG + UR) over curated probe sets, extracted from the **pre-noob (`QCDP`)** SigDF. Sex inference and `mean_intensity` deliberately use **non-noob** signal: noob only subtracts background and cannot recover weak signal, so it offers no accuracy benefit for chromosome-presence calls and would only push weak-signal samples toward zero, making ambiguous samples worse. The curated probe sets are:

- **chrY.clean** (314 probes): Non-PAR Y-chromosome probes from which cross-hybridizing probes have been removed. "Cross-hybridizing" means the probe sequence has significant homology to **any non-Y genomic region** (autosomal or chrX), not just chrX. Of ~548 non-PAR chrY probes on EPIC, 234 are excluded, leaving 314 with clean sex-dimorphic signal.

- **chrX.xlinked** (3,433 probes): X-inactivation-specific probes, also restricted to **non-PAR** regions of chrX. These are the subset of chrX probes that are subject to X-chromosome inactivation in females, identified from X-inactivation status annotations in the SeSAMe `probeInfo` database (derived from Carrel & Willard 2005 and Cotton et al. 2015 escape-from-inactivation catalogs). Of ~19,093 non-PAR chrX probes on EPIC, 3,433 are selected as consistently inactivated, excluding probes that escape X-inactivation (which show similar methylation in both sexes and reduce discriminating power).

Both sets are hardcoded in `sex_probe_lists.R` from `sesameDataGet("EPIC.probeInfo")`.

## 3. Stage 2: Flagging and QC Report (`prep_single()`)

### 3.1 Sample flagging (quality layer)

Samples are **flagged** (not removed) based on detection rate and mean intensity — the sample *quality* layer. Sex mismatches and age outliers belong to the *identity* layer and do NOT contribute to `flagged`: a swapped sample is often high quality, so folding identity into the quality flag would hide it.

### 3.2 Probe flagging

Probes are flagged in `exclude_probes.csv`: low call rate, sex chromosome, SNP, and other non-cg/ch.

### 3.3 QC diagnostic report (PDF)

**Page 1 — Detection rate**: Bar chart of `frac_dt`, flagged in red.

**Page 2 — MDS**: Classical MDS on the **top N most variable** autosomal cg/ch probes (ranked by `matrixStats::rowVars()`, highest variance first; configurable via `n_top_variable`, default 100,000) with complete data across ALL samples.

**Page 3 — Intensity histogram**: Colored by MDS outlier status.

**Page 4 — Probe failure rate**: All probes, all samples.

**Page 5 — Beta density**: Per-sample density curves. Probes subsampled to 100,000 by **random sampling** (uniform without replacement). Lines colored by QC status (OK / high probe failure / low intensity).

**Page 6 — Sex check** *(identity layer)*: chrX.xlinked vs. chrY.clean intensity scatter. Optimized chrY threshold minimizes total absolute residuals from within-cluster regressions. Confidence bands are drawn **parallel to each cluster's regression line**, extending ±5 SD (from the within-cluster orthogonal distances) perpendicular to the line.

*Edge case*: If a sample falls below the chrY threshold but within 5 SD of the male regression line, the crude threshold is used as the tiebreaker — the sample is classified as female.

> **Visual inspection recommended.** The sex check algorithm is generally robust, but edge cases (XXY karyotypes, mosaic sex chromosome loss, contaminated samples) can produce ambiguous results. Always inspect the sex check plot. Samples flagged as mismatches or unclear should be investigated individually — examine their position relative to both bands and consider re-extracting metadata.

**Page 7 — Age check** *(identity layer)*: Reported vs. Horvath predicted age. Confidence band extends ±3 SD perpendicular to the best-fit regression line. The clock runs on unmasked, noob-corrected betas.

> **Visual inspection recommended.** The Horvath clock was trained primarily on blood and brain tissues and performs less well on other tissue types (e.g., muscle, r ≈ 0.70). The 3 SD threshold may flag samples that are biologically normal for the tissue. Review the plot to distinguish genuine annotation errors from expected tissue-specific clock behavior.

**Pages 8–9 — PCA**: Scree plot followed by 6 PC scatter plots condensed onto 2 pages (3 per page via `gridExtra::grid.arrange()`). PCA uses the **top N most variable** autosomal cg/ch probes (same `n_top_variable`, non-flagged samples only, `prcomp(scale. = TRUE)`).

For each of the first 6 PCs, the algorithm tests every non-QC sample sheet column for association with that PC's scores:

- **Continuous variables** (numeric columns): Pearson correlation coefficient is computed between the PC scores and the variable. The variable with the highest |r| is selected.
- **Categorical variables** (factor/character columns with 2–20 unique levels): A one-way ANOVA F-test is computed (`lm(PC_scores ~ factor(variable))`). The variable with the smallest p-value is selected (must be p < 0.05 to qualify).

Continuous and categorical candidates compete on a common scale: |r| is compared as `-abs(r)`, and ANOVA p-value is compared as `log10(p)`. The variable with the best score (highest |r| or lowest p-value) is shown as the plot title and used for point coloring (viridis scale for continuous, discrete colors for categorical). Columns excluded from testing: `sample_id`, `Basename`, `batch_folder`, `horvath_age`, `detected_platform`, SeSAMe QC statistics, and sex chromosome intensities.

### 3.4 EpiDISH cell-type deconvolution

EpiDISH is run as a **purity / contamination check**. For blood-derived tissues — whole blood, PBMC, buffy coat, and isolated/sorted blood cell populations — `run_epidish()` deconvolves the **unmasked** beta matrix. A cleanly sorted population should be dominated by a single cell type; an unexpectedly mixed profile, or immune signal in a putatively non-immune sample, indicates contamination.

The default reference is **`cent12CT.m`** — 12 leukocyte subtypes (naive/memory CD4T, naive/memory CD8T, naive/memory B, Treg, NK, monocytes, neutrophils, eosinophils, basophils). The reference is configurable via the `epidish_reference` option, which accepts either a single panel name or a character vector of names (each is run, and the results are returned as a named list).

EpiDISH's RPC method constrains estimates to **sum to 1**. This is only meaningful for immune/blood tissue: running an immune-only reference on a non-immune tissue would force the estimates into a meaningless simplex regardless of the true composition. Deconvolution is therefore restricted to the blood/immune allowlist (`epidish_cell_types`); it is **not** run on non-blood tissue.

#### Reference panels

| Reference | Cell types |
|-----------|-----------|
| `cent12CT.m` *(default)* | 12 leukocyte subtypes (naive/memory CD4T & CD8T & B, Treg, NK, Mono, Neutro, Eosino, Baso) |
| `centDHSbloodDMC.m` | 7 types: B, NK, CD4T, CD8T, Mono, Neutro, Eosino |
| `centEpiFibIC.m` | Epithelial, fibroblast, total immune (solid tissue) |
| `centBloodSub.m` | Blood sub-types (recent EpiDISH versions) |

#### Runtime robustness

methylQC does not silently install or update packages. Instead:

- `DESCRIPTION` declares minimum versions (`sesame >= 1.20.0`, `EpiDISH >= 2.18.0`, etc.).
- `check_dependencies()` verifies the R version, package versions, the SeSAMe data cache, and the presence of the default EpiDISH reference, printing the exact fix command for anything missing.
- At runtime, `resolve_epidish_reference()` confirms the requested reference panel exists in the installed EpiDISH and stops with an actionable message if not.

#### Non-blood tissue deconvolution

For non-blood tissues, run deconvolution manually with an appropriate reference:

```r
library(EpiDISH)
betas <- readRDS("results/betas_all.rds")

# Custom reference for any tissue:
custom_ref <- as.matrix(read.csv("my_reference.csv", row.names = 1))
overlap <- intersect(rownames(betas), rownames(custom_ref))
result <- epidish(beta.m = betas[overlap, ], ref.m = custom_ref[overlap, ],
                  method = "RPC")$estF
```

### 3.5 Consolidated sample sheet

A single `sample_sheet.csv` contains ALL samples with columns added by the pipeline:

| Column | Type | How computed |
|--------|------|--------------|
| `flagged` | logical | `frac_dt` < threshold or `mean_intensity` < threshold |
| `flag_reason` | character | e.g., `"low_detection(0.831)"` |
| `reported_sex_normalized` | character | M/F/NA via `normalize_sex()` |
| `inferred_sex_intensity` | character | F/M/Unclear from intensity-based clustering |
| `sex_mismatch` | logical | Reported ≠ inferred AND both are confident |
| `sex_unclear` | logical | Outside both 5 SD confidence bands |
| `age_outlier` | logical | > 3 SD from age regression line |
| `reported_age_parsed` | numeric | Robust parsing of "45y", "NA", etc. |
| `B`, `NK`, `CD4Tnv`, `CD4Tmem`, `Treg`, `CD8Tnv`, `CD8Tmem`, `Mono`, `Neu`, `Eos`, `Baso`, ... | numeric | EpiDISH proportions, one column per reference cell type (blood/immune tissues only; columns depend on the reference panel) |

## 4. Applying masks and filtering downstream (`apply_mask()`)

### Basic usage

```r
betas <- readRDS("results/betas_all.rds")
mask  <- readRDS("results/mask_all.rds")
detP  <- readRDS("results/detP_all.rds")
excl_p <- read.csv("results/exclude_probes.csv")
ss     <- read.csv("results/sample_sheet.csv")
```

### Common workflows

```r
# --- Standard EWAS: autosomal CpG, mask + detP, exclude flagged samples ---
betas_ewas <- apply_mask(
  betas,
  mask = mask,
  detP = detP, detP_thresh = 0.05,
  exclude_probes = excl_p$probe_id,
  exclude_samples = ss$sample_id[ss$flagged],
  exclude_sex = TRUE, platform = "EPIC",
  probe_types = "cg"
)

# --- Include non-CpG methylation probes ---
betas_all_meth <- apply_mask(
  betas,
  mask = mask, detP = detP,
  exclude_probes = excl_p$probe_id,
  exclude_samples = ss$sample_id[ss$flagged],
  exclude_sex = TRUE, platform = "EPIC",
  probe_types = c("cg", "ch")
)

# --- Keep sex chromosome probes (e.g., for X-inactivation analysis) ---
betas_with_sex <- apply_mask(
  betas,
  mask = mask, detP = detP,
  exclude_probes = excl_p$probe_id,
  exclude_samples = ss$sample_id[ss$flagged],
  exclude_sex = FALSE,              # <-- retain sex chrom probes
  probe_types = "cg"
)

# --- Sex chromosome probes ONLY ---
# First apply mask/detP, then keep only probes on sex chromosomes
betas_sex_only <- apply_mask(
  betas,
  mask = mask, detP = detP,
  exclude_sex = FALSE,
  probe_types = "cg"
)
sex_probes <- get_sex_probes("EPIC")
betas_sex_only <- betas_sex_only[rownames(betas_sex_only) %in% sex_probes, ]

# --- SNP probes only (for identity verification) ---
betas_snp <- apply_mask(betas, exclude_sex = FALSE, platform = "EPIC",
                         probe_types = "rs")

# --- Minimal: just quality mask, keep everything ---
betas_masked <- apply_mask(betas, mask = mask, exclude_sex = FALSE)

# --- Strict detection threshold + imputation ---
betas_strict <- apply_mask(
  betas,
  mask = mask,
  detP = detP, detP_thresh = 0.01,  # stricter than default 0.05
  exclude_probes = excl_p$probe_id,
  exclude_sex = TRUE, platform = "EPIC",
  probe_types = "cg",
  impute = TRUE, knn_k = 50, knn_var_probes = 30000
)

# --- Apply the same filtering to M-values ---
mvals <- readRDS("results/mvals_all.rds")
mvals_clean <- apply_mask(
  mvals,
  mask = mask, detP = detP,
  exclude_probes = excl_p$probe_id,
  exclude_samples = ss$sample_id[ss$flagged],
  exclude_sex = TRUE, platform = "EPIC",
  probe_types = "cg"
)
```

### Order of operations

1. Quality mask → set masked values to NA
2. Detection p-value → set values with p > threshold to NA
3. Remove excluded samples (columns)
4. Remove excluded probes (rows)
5. Remove sex chromosome probes (if `exclude_sex = TRUE`)
6. Filter by probe type prefix
7. Similarity-based k-NN imputation (if `impute = TRUE`)

### Imputation algorithm

`apply_mask(..., impute = TRUE)` runs similarity-based k-NN imputation, replacing the earlier `impute::impute.knn` approach (the `impute` package is no longer a dependency):

1. **Variable-probe restriction.** Sample-to-sample distances are computed on the top `knn_var_probes` (default 30,000) most-variable probes. The bulk of the array is near-constant and would dominate an all-probe distance; restricting to variable probes is what makes "nearest neighbour" biologically meaningful.
2. **Distance.** Pairwise **Euclidean** distance between samples, computed over the probes both samples have non-missing and rescaled to a full-length distance so partial overlap is handled.
3. **Imputation.** Each missing value is filled with the **inverse-distance-weighted mean** of the `knn_k` (default 50) nearest-neighbour samples' values at that probe. Neighbours that are themselves missing at that probe are skipped.
4. **No fabrication.** There is no constant fallback. A value missing in a sample and in all of its usable neighbours is left as `NA`.

### Probe retention reference

| Goal | `exclude_sex` | `probe_types` |
|------|---------------|---------------|
| Autosomal CpG (standard EWAS) | `TRUE` | `"cg"` |
| All methylation (CpG + non-CpG) | `TRUE` | `c("cg", "ch")` |
| Include sex chromosomes | `FALSE` | `"cg"` or `c("cg", "ch")` |
| Sex chromosome probes only | `FALSE` + post-filter | `"cg"` |
| SNP probes only | `FALSE` | `"rs"` |
| Everything | `FALSE` | `NULL` |

## 5. Configuration

### Column names

| Parameter | Default | Aliases |
|-----------|---------|---------|
| `sample_id_col` | sample_id | Sample_ID, SampleID, Sample_Name |
| `donor_col` | Donor | Subject, Participant, SubjectID, DonorID |
| `reported_sex_col` | Reported_Sex | Sex, gender, Gender, Female, Male |
| `reported_age_col` | Age | age_years, AgeYears, age_at_sample |
| `batch_col` | batch_folder | Batch, Plate, Slide, Chip |
| `cell_col` | Cell | Cell_Type, Tissue, Sample_Type, Source |

### QC thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `prep_code` | `"QCDPB"` | SeSAMe openSesame prep code |
| `sample_call_min` | 0.95 | Minimum frac_dt |
| `probe_call_min` | 0.95 | Minimum probe call rate |
| `intensity_min` | 1300 | Minimum mean intensity |
| `n_top_variable` | 100000 | Most-variable probes for MDS/PCA |
| `knn_k` | 50 | Nearest-neighbour samples for imputation |
| `knn_var_probes` | 30000 | Most-variable probes for imputation distances |
| `save_sdfs` | TRUE | Save the final SigDF list (`sdfs_all.rds`) |

### EpiDISH and EPICv2 options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `epidish_reference` | `"cent12CT.m"` | Reference panel name, or a vector of names |
| `epidish_method` | `"RPC"` | EpiDISH deconvolution method |
| `epicv2_harmonize` | `FALSE` | Enable EPICv2 de-duplication + ID harmonization |
| `epicv2_target` | `"EPIC"` | Harmonization target: `"EPIC"` or `"HM450"` |

### Changing thresholds

```r
# Before running the pipeline
methylQC_set(sample_call_min = 0.90, intensity_min = 1000)
pipeline(in_dir = "data/idats", out_dir = "results")

# Or pass as arguments (temporary)
pipeline(in_dir = "data/idats", out_dir = "results",
         sample_call_min = 0.90)

# Re-run Stage 2 only with different thresholds
methylQC_set(sample_call_min = 0.85)
prep(stage1_dir = "results")

# Change detection p-value threshold at apply_mask() time
betas_strict <- apply_mask(betas, detP = detP, detP_thresh = 0.01)

# Check current settings
str(methylQC_options())

# Reset to defaults
methylQC_reset()
```

## 6. EPICv2 de-duplication and probe ID harmonization

EPICv2 targets many CpGs with multiple replicate probes — probe IDs sharing a `cg`/`ch` prefix but differing by a design suffix (Infinium chemistry, strand, replicate index). For a cohort analysis every sample must use the **same physical probe** for a given CpG; otherwise probe-design differences become a source of technical variation that confounds downstream comparisons.

This behaviour is controlled by a single option, **off by default**:

```r
methylQC_set(epicv2_harmonize = TRUE,
             epicv2_target    = "EPIC")   # or "HM450"
```

### 6.1 Algorithm

When `epicv2_harmonize = TRUE` and the detected platform is EPICv2, Stage 1 performs two steps after preprocessing:

**Step 1 — Cohort-consistent de-duplication.** For each replicated CpG (grouped by the `cg`/`ch` prefix, stripping the suffix), the per-probe **cross-sample detection failure rate** is computed (`mean(detP > detP_thresh)` across all samples). The probe with the lowest failure rate is kept; ties are broken by taking the first probe in ID order. Only `cg`/`ch` probes are de-duplicated — `rs` (SNP) and control probes are left untouched.

This is a **cohort-level** decision: every sample ends up using the same physical probe for a given CpG. It differs from SeSAMe's `collapseToPfx(method = "minPvalue")`, which selects per-sample and can therefore assign different replicate probes to the same CpG in different samples.

**Step 2 — Probe ID harmonization.** Because de-duplication has already removed every replicate, the surviving EPICv2 probe IDs are lifted to the target legacy platform (`EPIC` by default, or `HM450`) with SeSAMe's `mLiftOver`. With no replicates left to collapse, this step only renames probe IDs. The conversion mappings are the ones shipped in the `sesameData` cache (no external file, no bundled data).

A per-CpG kept/dropped log is written to `epicv2_dedup_log.csv` (columns: `cpg`, `kept_probe`, `kept_fail_rate`, `dropped_probes`, `n_replicates`). There is no console dump.

When harmonization runs inside the pipeline, the beta, M-value, mask, and detection p-value matrices are all re-aligned to the harmonized probe set so they stay mutually consistent.

### 6.2 Standalone function

`harmonize_epicv2()` applies the same de-duplication + harmonization to matrices you already have:

```r
harmonized <- harmonize_epicv2(
  mat_path        = "results/betas_all.rds",   # beta OR M-value matrix
  detP_path       = "results/detP_all.rds",    # detection p-value matrix
  target_platform = "EPIC",                    # or "HM450"
  dedup_log_path  = "results/epicv2_dedup_log.csv")
```

**Input format.** `mat` and `detP` are matrices with **probes in rows, samples in columns**, EPICv2 probe IDs as row names. The detection p-value matrix is **required** — the de-duplication criterion is the cross-sample detection failure rate, which cannot be recovered from an unmasked beta matrix. Both inputs are exactly what Stage 1 writes. To generate them directly from SeSAMe:

```r
library(sesame)
sdfs  <- openSesame(idat_dir, prep = "QCDP", func = NULL)   # list of SigDFs
betas <- do.call(cbind, lapply(sdfs, getBetas, mask = FALSE))
detP  <- do.call(cbind, lapply(sdfs, pOOBAH, return.pval = TRUE))
```

### 6.3 Interaction with the Horvath clock

Harmonizing EPICv2 IDs to EPIC means the hard-coded 353-CpG Horvath clock applies directly, without the per-platform model variants that would otherwise be needed for EPICv2's renamed probes. Clock CpGs genuinely absent from the platform are still handled by zero-shot reference imputation (§2.5).

## 7. SNP identity verification (`check_snp_identity()`)

Uses rs probes (auto-detected; 59 on EPIC/EPICv2, 65 on HM450). MDS on Euclidean distances of continuous beta values (no genotype calling). Per-participant centroids → nearest-centroid flagging. This is part of the **sample identity** QC layer.

## 8. Design decisions

**Unmasked matrices + separate masks**: Signal and QC are decoupled. Users apply their own thresholds.

**Flag, don't filter**: Exclusion lists are advisory.

**Quality vs. identity QC are separate**: A swapped sample can be high quality; folding identity checks into the quality flag would hide real swaps. The two layers are reported independently.

**QCDPB with P before B**: detection p-values are computed pre-noob (correct null distribution), while betas/M-values are noob-corrected.

**Sex probes separate from probe_types**: `exclude_sex` (default TRUE) handles sex chromosomes independently. Sex probes are cg probes — without this parameter they'd be indistinguishable from autosomal cg probes.

**Intensity-based sex inference on non-noob signal**: Uses curated non-cross-hybridizing, non-PAR probe sets with a data-driven threshold. noob cannot recover weak signal, so non-noob intensities are used.

**Cohort-consistent EPICv2 de-duplication**: Replicate selection is a cohort-level decision (cross-sample detection failure), so every sample uses the same physical probe per CpG.

**Similarity-based k-NN imputation**: Neighbours are ranked on the most-variable probes, where "nearest" is meaningful. No constant fallback — unimputable values stay NA.

**No silent package installation**: `check_dependencies()` reports problems and the exact fix command; nothing is installed automatically (CRAN/Bioconductor policy, and safe on no-admin/HPC environments).

**No batch correction**: Analysis-specific. PCA plots help assess.

**Hardcoded Horvath clock and sex probe lists**: Avoids silent failures from missing SeSAMe cache data.

## 9. References

1. Zhou W, Laird PW, Snyder M. *Nucleic Acids Research*. 2018;46(20):e123.
2. Zhou W, et al. mLiftOver: harmonizing data across Infinium DNA methylation platforms. *Bioinformatics*. 2024;40(7):btae423.
3. Horvath S. *Genome Biology*. 2013;14(10):R115.
4. Teschendorff AE, et al. *BMC Bioinformatics*. 2017;18(1):105.
5. Carrel L, Willard HF. *Nature*. 2005;434(7031):400–404.
6. Cotton AM, et al. *Genome Research*. 2015;25(8):1091–1099.
7. Guintivano J, Aryee MJ, Kaminsky ZA. *Epigenetics*. 2013;8(3):290–302.
