# methylQC: Methods and Technical Documentation

**Version 1.2.0**

## 1. Overview

methylQC is a two-stage quality control pipeline for Illumina Infinium DNA methylation arrays (EPIC, EPICv2, 450K). It produces **raw (unmasked)** beta- and M-value matrices alongside separate quality mask and detection p-value matrices, giving users full control over which probe values to trust.

The pipeline **flags** samples and probes but does not automatically filter them. A user-facing `apply_mask()` function applies any combination of quality masking, detection p-value thresholds, sex chromosome exclusion, probe/sample exclusions, probe type filtering, and optional k-NN imputation.

## 2. Stage 1: Preprocessing (`qc()`)

### 2.1 IDAT discovery and sample sheet parsing

`discover_idats()` recursively scans the input directory for paired IDAT files. Sample sheets are auto-discovered and parsed with tab/comma auto-detection. The array platform is auto-detected from IDAT headers. Sentrix IDs are decomposed into Chip, Row, and Col.

### 2.2 Signal processing (openSesame)

Each sample is processed through SeSAMe's `openSesame(func = NULL, prep = "QCDP")`. The `func = NULL` argument returns the processed SigDF. The `prep = "QCDP"` pipeline applies four steps to the raw signal:

1. **Q — Quality masking**: Probes with design issues are flagged in the SigDF's `mask` column (multi-mapping, SNP-overlapping CpG/SBE sites, cross-hybridizing sequences). This is a fixed, platform-level mask independent of sample signal.
2. **C — Channel correction (NOOB)**: Normal-exponential out-of-band background subtraction.
3. **D — Dye-bias correction**: Non-linear Cy3/Cy5 channel efficiency correction.
4. **P — pOOBAH detection p-values**: Each probe's signal is compared to the out-of-band hybridization distribution. Probes with p > 0.05 are additionally flagged in the `mask` column. This is sample-specific.

After openSesame, beta values are extracted via `getBetas(sdf, mask = FALSE)`, which computes betas from the **corrected** intensities (noob + dye-bias already applied) but does **NOT** set masked probes to NA. The mask and detection p-values are extracted separately.

### 2.3 Output matrices

All matrices are probes (rows) × samples (columns) in `.rds` format. Row names are probe IDs; column names are sample IDs.

| File | Type | Description |
|------|------|-------------|
| `betas_all.rds` | numeric | Unmasked betas from corrected intensities. No NAs. |
| `mvals_all.rds` | numeric | `log2((beta + 1e-6) / (1 - beta + 1e-6))` |
| `mask_all.rds` | logical | TRUE = failed quality mask or detection |
| `detP_all.rds` | numeric | pOOBAH detection p-values |
| `snp_betas.rds` | numeric | **samples × probes** (transposed). Continuous `rs` probe betas. |

### 2.4 Probe types

Probes are identified by ID prefix. Counts vary by platform:

| Prefix | Type | EPIC | EPICv2 | HM450 |
|--------|------|------|--------|-------|
| `cg` | CpG methylation | ~846,000 | ~930,000 | ~482,000 |
| `ch` | Non-CpG methylation | ~2,800 | ~2,900 | ~3,100 |
| `rs` | SNP genotyping | 59 | 59 | 65 |
| other | Non-standard | ~635 | varies | varies |

The `rs` probe count is auto-detected by prefix (`grepl("^rs", ...)`), not hardcoded.

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

**Horvath epigenetic age (`horvath_age`)**: 353 CpG coefficients hardcoded from Horvath (2013). Uses **masked** beta values (`getBetas(sdf)`) for clock prediction to exclude unreliable values.

**Sex chromosome intensities (`sex_chrX_intensity`, `sex_chrY_intensity`)**: Median total intensity (MG + MR + UG + UR) over curated probe sets:

- **chrY.clean** (314 probes): Non-PAR Y-chromosome probes from which cross-hybridizing probes have been removed. "Cross-hybridizing" means the probe sequence has significant homology to **any non-Y genomic region** (autosomal or chrX). Of ~548 non-PAR chrY probes on EPIC, 234 are excluded, leaving 314 with clean sex-dimorphic signal.

- **chrX.xlinked** (3,433 probes): X-inactivation-specific probes, also restricted to **non-PAR** regions of chrX. These are the subset of chrX probes that are subject to X-chromosome inactivation in females, identified from X-inactivation status annotations in the SeSAMe `probeInfo` database (derived from Carrel & Willard 2005 and Cotton et al. 2015 escape-from-inactivation catalogs). Of ~19,093 non-PAR chrX probes on EPIC, 3,433 are selected as consistently inactivated, excluding probes that escape X-inactivation (which show similar methylation in both sexes and reduce discriminating power).

Both sets are hardcoded in `sex_probe_lists.R` from `sesameDataGet("EPIC.probeInfo")`.

## 3. Stage 2: Flagging and QC Report (`prep_single()`)

### 3.1 Sample flagging

Samples are **flagged** (not removed) based on detection rate and mean intensity. Sex mismatches and age outliers are reported but do NOT contribute to `flagged`.

### 3.2 Probe flagging

Probes are flagged in `exclude_probes.csv`: low call rate, sex chromosome, SNP, and other non-cg/ch.

### 3.3 QC diagnostic report (PDF)

**Page 1 — Detection rate**: Bar chart of `frac_dt`, flagged in red.

**Page 2 — MDS**: Classical MDS on the **top N most variable** autosomal cg/ch probes (ranked by `matrixStats::rowVars()`, highest variance first; configurable via `n_top_variable`, default 100,000) with complete data across ALL samples.

**Page 3 — Intensity histogram**: Colored by MDS outlier status.

**Page 4 — Probe failure rate**: All probes, all samples.

**Page 5 — Beta density**: Per-sample density curves. Probes subsampled to 100,000 by **random sampling** (uniform without replacement). Lines colored by QC status (OK / high probe failure / low intensity).

**Page 6 — Sex check**: chrX.xlinked vs. chrY.clean intensity scatter. Optimized chrY threshold minimizes total absolute residuals from within-cluster regressions. Confidence bands are drawn **parallel to each cluster's regression line**, extending ±5 SD (from the within-cluster orthogonal distances) perpendicular to the line.

*Edge case*: If a sample falls below the chrY threshold but within 5 SD of the male regression line, the crude threshold is used as the tiebreaker — the sample is classified as female.

> **Visual inspection recommended.** The sex check algorithm is generally robust, but edge cases (XXY karyotypes, mosaic sex chromosome loss, contaminated samples) can produce ambiguous results. Always inspect the sex check plot. Samples flagged as mismatches or unclear should be investigated individually — examine their position relative to both bands and consider re-extracting metadata.

**Page 7 — Age check**: Reported vs. Horvath predicted age. Confidence band extends ±3 SD perpendicular to the best-fit regression line.

> **Visual inspection recommended.** The Horvath clock was trained primarily on blood and brain tissues and performs less well on other tissue types (e.g., muscle, r ≈ 0.70). The 3 SD threshold may flag samples that are biologically normal for the tissue. Review the plot to distinguish genuine annotation errors from expected tissue-specific clock behavior.

**Pages 8–9 — PCA**: Scree plot followed by 6 PC scatter plots condensed onto 2 pages (3 per page via `gridExtra::grid.arrange()`). PCA uses the **top N most variable** autosomal cg/ch probes (same `n_top_variable`, non-flagged samples only, `prcomp(scale. = TRUE)`).

For each of the first 6 PCs, the algorithm tests every non-QC sample sheet column for association with that PC's scores:

- **Continuous variables** (numeric columns): Pearson correlation coefficient is computed between the PC scores and the variable. The variable with the highest |r| is selected.
- **Categorical variables** (factor/character columns with 2–20 unique levels): A one-way ANOVA F-test is computed (`lm(PC_scores ~ factor(variable))`). The variable with the smallest p-value is selected (must be p < 0.05 to qualify).

Continuous and categorical candidates compete on a common scale: |r| is compared as `-abs(r)`, and ANOVA p-value is compared as `log10(p)`. The variable with the best score (highest |r| or lowest p-value) is shown as the plot title and used for point coloring (viridis scale for continuous, discrete colors for categorical). Columns excluded from testing: `sample_id`, `Basename`, `batch_folder`, `horvath_age`, `detected_platform`, SeSAMe QC statistics, and sex chromosome intensities.

### 3.4 EpiDISH cell-type deconvolution

For blood-derived tissues, EpiDISH (RPC method, `centDHSbloodDMC.m` reference, 7 cell types) runs on the **raw unmasked** beta matrix.

#### Non-blood tissue deconvolution

For other tissues, run deconvolution manually with appropriate references:

```r
library(EpiDISH)
betas <- readRDS("results/betas_all.rds")

# Brain (neuronal vs. glial)
# Reference: Guintivano et al. 2013, Epigenetics 8(3):290-302
data(BrainDMC, package = "EpiDISH")
result <- epidish(beta.m = betas, ref.m = BrainDMC, method = "RPC")$estF

# Saliva / buccal
# Mixture of buccal epithelial + leukocytes. Use blood reference for
# leukocyte fraction; estimate epithelial as 1 - sum(leukocyte props).
# Dedicated reference: Middleton et al. 2022, Epigenetics 17(2)

# Skeletal muscle — no standard EpiDISH reference.
# See Segales et al. 2020, Epigenetics 15(1-2) for fiber-type markers.

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
| `B`, `NK`, `CD4T`, `CD8T`, `Mono`, `Neutro`, `Eosino` | numeric | EpiDISH proportions (blood tissues only) |

## 4. Applying masks and filtering downstream (`apply_mask()`)

### Basic usage

```r
betas <- readRDS("betas_all.rds")
mask  <- readRDS("mask_all.rds")
detP  <- readRDS("detP_all.rds")
excl_p <- read.csv("exclude_probes.csv")
ss     <- read.csv("sample_sheet.csv")
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
  impute = TRUE, knn_k = 10, chunk_size = 50000
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
7. k-NN imputation (if `impute = TRUE`)

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
| `sample_call_min` | 0.95 | Minimum frac_dt |
| `probe_call_min` | 0.95 | Minimum probe call rate |
| `intensity_min` | 1300 | Minimum mean intensity |
| `n_top_variable` | 100000 | Most-variable probes for MDS/PCA |

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

## 6. SNP identity verification (`check_snp_identity()`)

Uses rs probes (auto-detected; 59 on EPIC/EPICv2, 65 on HM450). MDS on Euclidean distances of continuous beta values (no genotype calling). Per-participant centroids → nearest-centroid flagging.

## 7. Design decisions

**Raw matrices + separate masks**: Signal and QC are decoupled. Users apply their own thresholds.

**Flag, don't filter**: Exclusion lists are advisory.

**Sex probes separate from probe_types**: `exclude_sex` (default TRUE) handles sex chromosomes independently. Sex probes are cg probes — without this parameter they'd be indistinguishable from autosomal cg probes.

**Intensity-based sex inference**: Uses curated non-cross-hybridizing, non-PAR probe sets with a data-driven threshold that adapts to each dataset. Transparent visual assessment on the sex check plot.

**Optional imputation**: Chunked k-NN via `apply_mask(..., impute = TRUE)`.

**No batch correction**: Analysis-specific. PCA plots help assess.

**Hardcoded Horvath clock and sex probe lists**: Avoids silent failures from missing SeSAMe cache data.

## 8. References

1. Zhou W, Laird PW, Snyder M. *Nucleic Acids Research*. 2018;46(20):e123.
2. Horvath S. *Genome Biology*. 2013;14(10):R115.
3. Teschendorff AE, et al. *BMC Bioinformatics*. 2017;18(1):105.
4. Carrel L, Willard HF. *Nature*. 2005;434(7031):400–404.
5. Cotton AM, et al. *Genome Research*. 2015;25(8):1091–1099.
6. Guintivano J, Aryee MJ, Kaminsky ZA. *Epigenetics*. 2013;8(3):290–302.
