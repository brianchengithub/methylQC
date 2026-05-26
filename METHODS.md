# methylQC: Methods and Technical Documentation

**Version 2.0.0**

---

## 0. What changed in v2.0.0

Breaking changes from v1.x. There are no deprecation shims.

- **No automatic exclusion.** v1.x wrote `exclude_samples.csv` and
  `exclude_probes.csv` and built a "filtered" beta matrix internally.
  v2.0.0 produces **flags only**. The user decides what to exclude and
  applies it via `cleanmat()`.
- **No probe-exclusion table at all.** `exclude_probes.csv`,
  `probe_call_rates.csv`, and the `low_call_rate` flag are gone.
  Probe-level masking is carried entirely by the `mask` and `detP`
  matrices.
- **`build_sample_exclusions()` is gone.** Its replacement, internal
  `qcflags()`, computes per-sample low-detection and low-intensity
  flags **only** — sex mismatch and age outliers are no longer rolled
  into the `flagged` column.
- **New exported utility `flagsamples()`** — given a detection p-value
  matrix and a call-rate threshold, writes a CSV of samples below the
  threshold. Advisory only.
- **QC report restructured to 9 pages** with new colour rules
  (see §3.3).
- **Function and option names shortened** to single tokens
  (no underscores) for the user-facing API. Internal helpers retain
  underscores.

API rename map for the user-facing surface:

| v1.x                              | v2.0.0          |
|-----------------------------------|-----------------|
| `apply_mask`                      | `cleanmat`      |
| `build_sample_exclusions`         | (internal `qcflags`) |
| `build_probe_exclusions`          | (removed)       |
| `check_sample_metadata`           | `checkmeta`     |
| `check_snp_identity`              | `snpcheck`      |
| `compute_qc_metrics`              | `qcmetrics`     |
| `compute_qc_metrics_streaming`    | `qcstream`      |
| `discover_idats`                  | `discover`      |
| `edit_exclusions`                 | (removed)       |
| `extract_snp_betas`               | `snpbetas`      |
| `log_capture`                     | `logcapture`    |
| `make_logger`                     | `makelog`       |
| `methylQC_defaults`               | `mqcdefaults`   |
| `methylQC_options`                | `mqcopts`       |
| `methylQC_reset`                  | `mqcreset`      |
| `methylQC_set`                    | `mqcset`        |
| `prep_single`                     | `prepcell`      |
| `run_epidish`                     | `rundish`       |
| `run_opensesame`                  | `runsesame`     |
| `write_qc_report`                 | `qcreport`      |
| (new)                             | `flagsamples`   |
| (new)                             | `checkdeps`     |
| (new)                             | `liftv2`        |

---

## 1. Overview

methylQC is a two-stage QC pipeline for Illumina Infinium DNA
methylation arrays (EPIC, EPICv2, 450K).

- **Stage 1 (`qc()`)** streams IDATs through SeSAMe's `openSesame`
  pipeline (QCDP), computes per-sample QC metrics, extracts SNP probe
  betas, and writes raw (unmasked) beta- and M-value matrices alongside
  separate quality-mask and detection p-value matrices.
- **Stage 2 (`prep()`/`prepcell()`)** flags samples for plot colouring,
  produces a multi-page QC PDF, runs EpiDISH cell-type deconvolution
  for blood-derived tissues, and writes a consolidated sample sheet.

**Nothing is removed automatically.** The pipeline flags; the user
applies exclusions via `cleanmat()`. The optional `flagsamples()`
utility produces an advisory CSV of samples whose call rate falls
below a user-chosen threshold.

---

## 2. Stage 1: Preprocessing (`qc()`)

### 2.1 IDAT discovery and sample sheet parsing

`discover(dir, sheetpattern, basecol, platform, logger)` recursively
scans the input directory for paired `_Grn.idat` / `_Red.idat` files.
For each folder containing IDATs it locates a sample sheet (regex
`sheetpattern`, default `"sample.*sheet.*\\.(csv|txt|tsv)$"`) with
tab-vs-comma auto-detection, or synthesises one from filenames. The
array platform is auto-detected from IDAT headers and checked for
consistency across folders. Sentrix-format IDs are decomposed into
`Chip`, `Row`, and `Col`.

### 2.2 Signal processing: openSesame and the QCDP pipeline

Each sample is processed through `sesame::openSesame(func = NULL,
prep = "QCDP")`. `func = NULL` returns the processed `SigDF`. The
`prep = "QCDP"` string is read left-to-right by SeSAMe: each letter is
a preprocessing step applied in order.

1. **Q — Quality masking.** SeSAMe sets the `mask` column to `TRUE`
   for probes with known design problems (multi-mapping, SNP overlap
   at the CpG or single-base-extension site, cross-hybridisation).
   This mask is platform-level and signal-independent.
2. **C — Channel correction (NOOB).** Normal-exponential out-of-band
   background subtraction. Sample-specific.
3. **D — Dye-bias correction.** Non-linear correction of Cy3/Cy5
   channel efficiency.
4. **P — pOOBAH detection p-values.** Each probe's signal is compared
   to the out-of-band hybridisation distribution; probes with
   `p > 0.05` are additionally flagged in the `mask` column.

Note that the `prep = "QCDP"` step order above is what SeSAMe runs.
methylQC adds a fifth step purely at extraction time — call it B for
"betas":

5. **B — Beta extraction with `mask = FALSE`.** After QCDP, methylQC
   calls `getBetas(sdf, mask = FALSE)` to produce **unmasked** betas
   from the corrected intensities. The `mask` column is preserved and
   returned as a separate matrix; the beta matrix itself contains no
   NAs from masking. This decoupling is the central design choice:
   downstream code applies masks (or doesn't) at its own discretion.

### 2.3 Output matrices

All matrices are probes (rows) × samples (columns), in `.rds`, with
probe IDs as row names and sample IDs as column names. The single
exception is `snp_betas.rds`, which is transposed (samples × rs probes).

| File              | Type    | Description                                                            |
|-------------------|---------|------------------------------------------------------------------------|
| `betas_all.rds`   | numeric | Unmasked betas from corrected intensities. No NAs from masking.        |
| `mvals_all.rds`   | numeric | `log2((beta + 1e-6) / (1 - beta + 1e-6))`.                             |
| `mask_all.rds`    | logical | `TRUE` = probe failed the quality mask **or** pOOBAH detection.        |
| `detP_all.rds`    | numeric | pOOBAH detection p-values.                                             |
| `snp_betas.rds`   | numeric | **Samples × rs probes**. Continuous rs-probe betas for identity check. |
| `sample_sheet.csv`| table   | Merged sample sheet + per-sample QC metrics.                           |
| `metadata.rds`    | list    | Platform, input dir, sample count, batch count, timestamp, pkg version.|

### 2.4 Probe types

Probes are classified by ID prefix (no hard-coded counts):

| Prefix | Type                  | Approx counts (EPIC / EPICv2 / 450K) |
|--------|-----------------------|--------------------------------------|
| `cg`   | CpG methylation       | ~846k / ~930k / ~482k                |
| `ch`   | Non-CpG methylation   | ~2.8k / ~2.9k / ~3.1k                |
| `rs`   | SNP genotyping        | 59 / 59 / 65                         |
| other  | Platform-specific     | ~635 / varies / varies               |

"Other" is defined by exclusion: any probe whose ID does not start
with `cg`, `ch`, or `rs`. On EPIC this includes `nv` (negative
verification) probes among others.

**Sex chromosome probes** are a subset of `cg` probes — all sex
probes have the `cg` prefix. In `cleanmat()` they are partitioned
out of `"cg"` into their own category `"sex"` so the user can include
or exclude them independently.

### 2.5 Per-sample QC metrics

Computed by `qcmetrics()` (in-memory list of SigDFs) or `qcstream()`
(streaming, one IDAT pair at a time). Both produce the same columns.

From `sesameQC_calcStats(funs = c("detection", "intensity", "numProbes",
"channel", "dyeBias"))`:

- **`frac_dt`** — fraction of probes passing detection (see §2.6 for
  how this differs from a per-sample call rate computed from the
  saved `detP_all.rds`).
- **`mean_intensity`** — mean total intensity across probes.
- Plus channel diagnostics, dye-bias diagnostics, and per-probe-type
  counts.

Additional columns computed by methylQC:

- **Probe decomposition:** `n_probes_total`, `n_masked_apriori`,
  `n_failed_detection`, `n_detected`, `n_usable`,
  `frac_apriori_mask`, `frac_detection_among_unmasked`, `frac_usable`.
  These split SeSAMe's combined mask into the quality-mask and
  detection-mask components.
- **`horvath_age`** — Horvath (2013) epigenetic age. Uses
  `getBetas(sdf)` (i.e. **masked** betas) so masked probes don't
  enter the clock; requires ≥200 clock probes present.
- **`sex_chrX_intensity` / `sex_chrY_intensity`** — median total
  intensity (`MG + MR + UG + UR`) over the curated probe sets
  `.sesame_chrX_xlinked` (3,433 non-PAR X-inactivated probes) and
  `.sesame_chrY_clean` (314 non-PAR Y-specific probes with
  cross-hybridising probes removed). Both lists are baked into
  `sex_probe_lists.R` to avoid silent failures when the SeSAMe data
  cache is missing.

`checkmeta(ss, qcdf, logger)` cross-checks reported sex and reported
age against the inferred values **for the log only**; it does not
produce flags.

### 2.6 `frac_dt` is not the same as `colMeans(detP <= 0.05)`

This distinction matters in v2.0.0 because `flagsamples()` computes
sample call rate from `detP_all.rds` directly, whereas the QC report
and sample sheet use SeSAMe's `frac_dt`.

- **`frac_dt`** is reported by `sesameQC_calcStats(funs = "detection")`
  and stored on the sample sheet. SeSAMe computes it from the SigDF's
  internal state; the exact denominator depends on the SeSAMe version
  but is generally the fraction of **non-quality-masked** probes that
  pass pOOBAH.
- **`colMeans(detP <= 0.05)`** computed on `detP_all.rds` uses **all
  probes**, regardless of quality mask. It is also strictly bounded by
  whatever pOOBAH actually wrote into `detP`.

These two values can differ by a few percent. `flagsamples()` returns
the latter; it is reproducible, independent of SeSAMe internals, and
the natural quantity to threshold for advisory exclusion. Compare to
`frac_dt` only if you want to understand a single sample, not to
build an exclusion list.

---

## 3. Stage 2: Flagging and QC Report (`prepcell()`)

### 3.1 Sample flagging — `qcflags()` (internal)

Per-sample flags used **only to colour QC plots**. Writes nothing to
disk, produces no exclusion list. Two flag types:

- **Low detection** — `frac_dt < samplemin` (default 0.95, strict `<`).
- **Low intensity** — `mean_intensity < intmin` (default 1300,
  strict `<`).

Sex mismatch and age outliers are deliberately **not** rolled into
`flagged`. They have their own dedicated columns in the consolidated
sample sheet and their own QC-report panels.

### 3.2 `flagsamples()` — advisory, opt-in CSV of low-call-rate samples

Standalone exported utility. Signature:

```r
flagsamples(detP,
            callrate = NULL,          # default = mqcopts()$samplemin (0.95)
            pthresh  = NULL,          # default = mqcopts()$detp      (0.05)
            csv      = "flagged_samples.csv",
            logger   = NULL)
```

For each sample, computes `colMeans(detP <= pthresh, na.rm = TRUE)`
(see §2.6) and writes the rows with call rate below `callrate` to
`csv`. Returns the full table invisibly. **Does not modify any
matrix and does not accept hand-picked IDs** — those go into
`cleanmat(dropsamples = …)`.

### 3.3 QC diagnostic report (PDF)

Nine pages, in this order. Colour rules per panel are listed
explicitly; legend visibility follows the same rules.

| Page | Panel                                | Colour by                                   |
|------|--------------------------------------|---------------------------------------------|
| 1    | Detection rate per sample (bar)      | Low-detection flag (`frac_dt < samplemin`)  |
| 2    | Per-probe sample-failure histogram (full range) | Single colour; Y axis capped to expose the tail |
| 3    | MDS on all non-sex cg/ch probes      | MDS outlier (4 SD from geometric median)    |
| 4    | Mean intensity per sample (histogram)| MDS outlier (precedence) > low-intensity > OK |
| 5    | Sample beta-value density            | Low-intensity flag                          |
| 6    | Sex check (chrX vs chrY intensity)   | Reported sex (F / M / n/a)                  |
| 7    | Reported vs Horvath age              | Sex-mismatch flag                           |
| 8    | Scree plot                           | (no colour-coding)                          |
| 9    | PC vs associated variable, 2×3 panel | The variable's levels/values                |

**Page 1 — Detection rate.** Horizontal bar per sample, ordered by
`frac_dt`. Dashed line at `samplemin`. Bars below the threshold are
red.

**Page 2 — Per-probe sample-failure histogram (full range, Y capped).**
Replaces v1.x's per-probe call-rate distribution. The histogram spans
the full 0–1 X axis in 40 equal-width bins so the shape of the
distribution — and in particular *where it starts to drop* — is
visible. The Y axis is capped at 1.5 × the tallest bin centred at
`fail_rate ≥ 0.05`; this deliberately clips the dominant near-zero
spike (probes that fail in almost no samples) while leaving every
other bar at full height. The clipped bin's count is reported in the
subtitle so the spike's magnitude is still known.

A dashed vertical line at `failmin` (default 0.10, configurable via
`mqcset(failmin = …)`) marks the CSV cutoff. Probes with
`fail_rate ≥ failmin` are written to `failed_probes.csv`. The CSV
behaviour is unchanged from earlier v2 drafts; only the plot is
broader now.

Failure rate is built from `detP` by default
(`rowMeans(detP > pthresh)`); if `detP` is unavailable, falls back to
`mask`, then to `is.na(betas)`. Quality-masked probes are excluded
from the histogram by default (`inclqual = FALSE`); set
`inclqual = TRUE` to include them.

**Page 3 — MDS.** Classical multidimensional scaling on **all**
complete-case cg/ch non-sex probes (v1.x took a top-variable subset
of 100,000 — v2.0.0 uses all). Distance matrix is Euclidean over
samples (`dist(t(betas))`); embedding via `cmdscale(d, k = 2)`.
Outliers are samples whose distance from the geometric median (via
Weiszfeld iteration) exceeds 4 standard deviations of the all-samples
distance distribution. A dashed circle marks the 4-SD ring.

**Page 4 — Mean intensity.** Histogram of `mean_intensity` with
samples grouped by colour. The colour is three-category with a strict
precedence rule:

```
MDS outlier  >  Low intensity  >  OK
```

A sample that is both an MDS outlier *and* below `intmin` is shown as
"MDS outlier" — MDS-outlier colour wins. Dashed vertical line at
`intmin`.

**Page 5 — Beta density.** One density curve per sample, overlaid on
a single panel. Density is computed with `density(x, from = 0, to = 1,
n = 512)` on each sample's cg/ch beta values (NAs excluded; values
clamped to `[0, 1]`). Lines are coloured by the low-intensity flag
(red for low-intensity, blue for OK). No legend.

**Page 6 — Sex check.** Scatter of `sex_chrX_intensity` (X)
vs `sex_chrY_intensity` (Y) coloured by reported sex (F / M / n/a).
Algorithm:

1. Search `chrY` quantiles (`probs = seq(0.15, 0.85, 0.01)`) for the
   threshold that minimises within-cluster sum-of-absolute-residuals
   from two regressions, one for `chrY <= thr` (female) and one for
   `chrY > thr` (male).
2. Fit `lm(chrY ~ chrX)` separately within each cluster.
3. For each sample, compute orthogonal distance to each line; classify
   as F or M if within 5 SD of one band and >5 SD of the other.
   Within both bands → use the crude threshold as tiebreaker; outside
   both → "Unclear".
4. Mismatch = reported ≠ inferred and inferred is not "Unclear".

Confidence bands (5 SD perpendicular distance) are drawn shaded around
each regression line. Visual inspection is recommended — XXY,
mosaic-loss-of-Y, and contaminated samples can land in edge regions.

**Page 7 — Age check.** Reported vs Horvath predicted age. Robust age
parser strips `"y"`, `"yrs"`, etc. Fits `lm(horvath ~ reported)`.
Outliers are samples whose orthogonal distance to the regression line
exceeds 3 SD of the residual distribution. Points are coloured by
the sex-mismatch flag (red for mismatch, blue for OK). Visual
inspection is recommended — the Horvath clock is trained on blood
and brain; expect lower fit in other tissues.

**Page 8 — Scree.** Bar plot of `% variance explained` for the first
`min(20, npc)` PCs. Title includes the actual number of probes used
in the PCA so the user can sanity-check the input size. This is the
only PCA panel that displays percent-variance.

**Page 9 — PC vs associated variable, 2×3 single page.** For each of
the first 6 PCs, the algorithm scans every non-QC sample-sheet column
and selects the one most associated with the PC's scores:

- **Numeric columns:** Pearson `|r|` with the PC. Variable maximising
  `|r|` wins.
- **Categorical columns (2–20 unique levels):** one-way ANOVA F-test
  `lm(PC ~ factor(var))`; variable minimising p-value wins, must be
  p < 0.05 to qualify.
- Continuous and categorical compete on a common scale: continuous
  uses `-abs(r)`, categorical uses `log10(p)`; lower wins.

Each panel renders the chosen variable on the x-axis (scatter for
numeric, jittered points + boxplot for categorical), with the
selection metric printed as subtitle (e.g. `|r|=0.812` or
`ANOVA p=3.1e-04`). No `% variance` is printed — that lives only on
the scree plot (page 8). Excluded from association testing:
`sample_id`, `Basename`, `batch_folder`, `sheet_path`, `horvath_age`,
`detected_platform`, raw IDAT columns, methylQC-internal QC columns,
and sex chromosome intensities (regex
`"^frac_|^n_|^num_|^mean_|^na_|^med[RG]$|^top[RG]$|^InfI_|^RG|^sex_chr"`).

PCA itself uses `prcomp(t(betas), scale. = TRUE)` on the top
`ntop` (default 100,000) most-variable probes among cg/ch non-sex
probes, complete cases only. PC scores for all PCs are written to
`pc_scores.csv`.

### 3.4 EpiDISH cell-type deconvolution

`rundish()` runs EpiDISH (`method = "RPC"` by default) using the
reference resolved from `mqcopts()$dishref`. Default reference is the
Salas et al. 2022 12-cell immune panel
([Nat Commun 13:761, 2022](https://doi.org/10.1038/s41467-021-27864-7)),
split by platform:

| Platform | Reference | Cells | Cell-type labels |
|----------|-----------|-------|------------------|
| EPIC     | `cent12CT.m`    | 12 | CD4Tnv, CD4Tmem, CD8Tnv, CD8Tmem, Treg, Bnv, Bmem, NK, Mono, Neu, Eos, Bas |
| 450k     | `cent12CT450k.m`| 12 | same 12 cell types |

To revert to the original 7-cell Houseman/DHS reference, set
`mqcset(dishref = "centDHSbloodDMC.m")` — that matrix is dual-platform
and produces columns `B, NK, CD4T, CD8T, Mono, Neutro, Eosino`.

**No built-in tissue allowlist.** Earlier methylQC drafts hard-coded
a list of blood-tissue tokens (`PBMC, WB, WBC, …`) and silently
skipped EpiDISH for anything else. v2.0.0 dropped that: cell-type
labels are user-defined and not reliable ("BM" can mean bone marrow
or B-memory; "Mono" can mean monocyte or whatever the lab called it).
The only gate is data-level: `rundish()` errors if the reference's
CpGs don't sufficiently overlap (< 50) the input matrix.

**Two calling forms** of the same function:

```r
# Form A (in-memory): returns the proportion matrix, no I/O
props <- rundish(betas, platform = "EPIC")

# Form B (on-disk): writes cell_proportions.csv AND merges columns
# into sample_sheet.csv on disk; returns proportions invisibly
rundish("results/PBMC")
rundish("results", celltype = "BM")          # filter by sample-sheet label
rundish("results", samples  = c("S1","S2"))  # filter by ID
rundish("results/PBMC", ref = "centDHSbloodDMC.m")   # override ref
rundish("results/PBMC", method = "CBS")              # override method
```

The on-disk form is dispatched when the first argument is a length-1
string that points to an existing directory. It loads
`betas_all.rds`, optionally subsets, runs the in-memory primitive,
writes `cell_proportions.csv`, and left-joins the proportion columns
into `sample_sheet.csv` (removing any stale prop columns from a prior
run before merging).

**In-pipeline vs standalone — same end state.** When EpiDISH runs as
part of Stage 2 (`prep(..., dish = TRUE)`, the default), the
proportions get merged into the per-cell-type `sample_sheet.csv` via
`build_consolidated_sample_sheet()`. When EpiDISH runs standalone
(`rundish("dir/")`) the proportions get merged via the same
left-join logic. Both paths produce equivalent merged sheets.

**Custom references.** Pass any EpiDISH reference name as a string
via `ref`, or set it as the default:

```r
# Brain (Guintivano et al. 2013)
rundish("results/brain", ref = "BrainDMC")

# Run a non-EpiDISH custom reference outside methylQC
library(EpiDISH)
betas <- readRDS("results/betas_all.rds")
custom_ref <- as.matrix(read.csv("my_reference.csv", row.names = 1))
overlap <- intersect(rownames(betas), rownames(custom_ref))
result <- epidish(beta.m = betas[overlap, ],
                  ref.m  = custom_ref[overlap, ],
                  method = "RPC")$estF
```

**Opting out of in-pipeline EpiDISH.** Pass `dish = FALSE` to
`prep()` or `pipeline()`; then run `rundish()` standalone afterwards
on whichever directories you want, with whichever reference.

### 3.5 Consolidated sample sheet

A single `sample_sheet.csv` per Stage 2 slice (one per cell type if
`bycell = TRUE`) contains all input rows with these added columns:

| Column                     | Type       | Source                                              |
|----------------------------|------------|-----------------------------------------------------|
| `flagged`                  | logical    | `frac_dt < samplemin` OR `mean_intensity < intmin`  |
| `flag_reason`              | character  | e.g. `"low_detection(0.831);low_intensity(1102)"`   |
| `reported_sex_normalized`  | character  | `normalize_sex()` of `sexcol`                       |
| `inferred_sex_intensity`   | character  | F / M / Unclear from page-6 sex check               |
| `sex_mismatch`             | logical    | reported ≠ inferred, both confident                 |
| `sex_unclear`              | logical    | outside both 5 SD bands                             |
| `age_outlier`              | logical    | >3 SD from age regression                           |
| `reported_age_parsed`      | numeric    | robust parser of `agecol`                           |
| `B, NK, CD4T, CD8T, Mono, Neutro, Eosino` | numeric | EpiDISH (blood only)                  |

---

## 4. Applying QC decisions: `cleanmat()`

`cleanmat()` is the single primitive for turning the raw Stage 1
matrices into an analysis-ready matrix. It does five things in a
fixed order; the user controls each by passing the relevant
argument, but methylQC never picks IDs on the user's behalf.

The function was named `applymask()` in earlier v2 drafts. It was
renamed when probe and sample exclusion were folded in, because
"masking" is no longer the only thing it does.

### 4.1 Signature

```r
cleanmat(mat,
         mask        = NULL,          # logical matrix; TRUE = mask out
         detP        = NULL,          # numeric matrix; mask where detP > pthresh
         pthresh     = NULL,          # default mqcopts()$detp (0.05)
         dropprobes  = NULL,          # NULL | char vec of probe IDs | CSV path
         dropsamples = NULL,          # NULL | char vec of sample IDs | CSV path
         probes      = c("cg", "ch"), # categories to keep
         platform    = NULL,          # required if "cg" or "sex" in probes
         impute      = FALSE,
         knnk        = 10L,
         chunk       = 50000L)
```

`dropprobes` and `dropsamples` each accept:

- **`NULL`** — no exclusion of that kind.
- **A character vector of IDs** — explicit list. `dropprobes` matches
  against `rownames(mat)`; `dropsamples` against `colnames(mat)`.
- **A length-1 string that is a path to a CSV file** — the function
  reads it and pulls IDs from the `probe_id` (or `sample_id`) column;
  if no such column is present, it falls back to the first column.
  This is what consumes methylQC's own `failed_probes.csv` and
  `flagged_samples.csv` directly.

Disambiguation: if the argument is length 1 AND `file.exists()`
returns TRUE, it is treated as a path. Otherwise it is treated as
an ID vector. A length-1 ID like `"cg00000029"` is therefore
interpreted as an ID, not a path, as long as no file by that name
exists in the working directory.

`probes` is a character vector of probe categories to **keep**:

| Token   | Meaning                                                    |
|---------|------------------------------------------------------------|
| `"cg"`  | Autosomal CpG (cg probes NOT on chrX/chrY).                |
| `"ch"`  | Non-CpG methylation (ch probes).                           |
| `"sex"` | Sex chromosome probes (chrX, chrY) of any prefix.          |
| `"snp"` | rs probes.                                                 |
| `"other"`| Everything else (e.g. `nv`).                              |

Pass `NULL` to keep every probe (no category filter). The five
categories form a strict partition — `"cg"` always excludes
sex-chromosome cg probes, and `"sex"` always excludes from the
others.

`platform` is required iff `"cg"` or `"sex"` is in `probes`, because
the sex-vs-autosomal split is read from the SeSAMe manifest.

### 4.2 Order of operations

1. **Quality mask** → set positions where `mask == TRUE` to `NA`.
2. **Detection** → set positions where `detP > pthresh` to `NA`.
3. **Drop samples** → resolve `dropsamples` (vector or CSV) and drop
   matching columns.
4. **Drop probes** → resolve `dropprobes` (vector or CSV) and drop
   matching rows.
5. **Probe category filter** → keep only probes in any of the
   `probes` categories.
6. **k-NN imputation** (optional) → chunked `impute::impute.knn`
   with `k = knnk`, chunk size `chunk`.

Probe and sample exclusions happen **before** the category filter so
the user's explicit drops take precedence and the category filter
operates on the post-exclusion matrix.

### 4.3 Common workflows

```r
betas <- readRDS("results/betas_all.rds")
mask  <- readRDS("results/mask_all.rds")
detP  <- readRDS("results/detP_all.rds")

# The "feed cleanmat the QC artifacts directly" workflow. This is the
# happy path: cleanmat reads probe IDs from failed_probes.csv and
# sample IDs from flagged_samples.csv on disk, applies masks, drops
# them, keeps autosomal CpG, and imputes.
betas_ewas <- cleanmat(
  betas, mask = mask, detP = detP,
  dropprobes  = "results/failed_probes.csv",
  dropsamples = "results/flagged_samples.csv",
  probes      = "cg",
  platform    = "EPIC",
  impute      = TRUE)

# Mix forms: probe IDs from a CSV, sample IDs from an in-memory vector
ss <- read.csv("results/sample_sheet.csv")
betas_clean <- cleanmat(
  betas, mask = mask, detP = detP,
  dropprobes  = "results/failed_probes.csv",
  dropsamples = ss$sample_id[ss$flagged],
  probes      = "cg", platform = "EPIC")

# CpG + non-CpG methylation, autosomal
betas_all_meth <- cleanmat(
  betas, mask = mask, detP = detP,
  probes = c("cg", "ch"), platform = "EPIC")

# Keep sex chromosome probes alongside autosomes
betas_with_sex <- cleanmat(
  betas, mask = mask, detP = detP,
  probes = c("cg", "sex"), platform = "EPIC")

# Sex chromosome probes only
betas_sex_only <- cleanmat(
  betas, mask = mask, detP = detP,
  probes = "sex", platform = "EPIC")

# SNP probes only (for identity work) — no platform needed
betas_snp <- cleanmat(
  betas, probes = "snp")

# Minimal: just the quality mask, no category filter, no exclusions
betas_masked <- cleanmat(
  betas, mask = mask, probes = NULL)

# Stricter detection threshold + imputation
betas_strict <- cleanmat(
  betas, mask = mask, detP = detP, pthresh = 0.01,
  dropprobes  = "results/failed_probes.csv",
  dropsamples = "results/flagged_samples.csv",
  probes      = "cg", platform = "EPIC",
  impute      = TRUE, knnk = 10, chunk = 50000)

# Same filter chain on M-values
mvals <- readRDS("results/mvals_all.rds")
mvals_clean <- cleanmat(
  mvals, mask = mask, detP = detP,
  dropprobes  = "results/failed_probes.csv",
  dropsamples = "results/flagged_samples.csv",
  probes      = "cg", platform = "EPIC")
```

---

## 5. Configuration

All tunable parameters are stored as package options under the
`methylQC.` prefix. Inspect with `mqcopts()`, override with
`mqcset()`, reset with `mqcreset()`.

### 5.1 Column-name options

| Option         | Default          | Aliases                                                  |
|----------------|------------------|----------------------------------------------------------|
| `idcol`        | `sample_id`      | `Sample_ID`, `SampleID`, `Sample_Name`                   |
| `basecol`      | `Basename`       | (single value)                                           |
| `donorcol`     | `Donor`          | `Subject`, `Participant`, `SubjectID`, `DonorID`         |
| `sexcol`       | `Reported_Sex`   | `Sex`, `gender`, `Gender`, `Female`, `Male`              |
| `agecol`       | `Age`            | `age_years`, `AgeYears`, `age_at_sample`                 |
| `batchcol`     | `batch_folder`   | `Batch`, `Plate`, `Slide`, `Chip`                        |
| `cellcol`      | `Cell`           | `Cell_Type`, `Tissue`, `Sample_Type`, `Source`           |
| `sheetpattern` | `"sample.*sheet.*\\.(csv|txt|tsv)$"` | (regex)                                |

Each `*col` has a matching `*aliases` option for the search list when
the primary name is absent (e.g. `donoraliases`, `sexaliases`).

### 5.2 Threshold options

| Option       | Default | Meaning                                                                                       |
|--------------|---------|-----------------------------------------------------------------------------------------------|
| `samplemin`  | 0.95    | Sample call-rate threshold. Sample is flagged low-detection if `frac_dt < samplemin` (strict).|
| `probemin`   | 0.95    | Probe pass-rate threshold (used in the inline probe-breakdown console summary).               |
| `intmin`     | 1300    | Sample is flagged low-intensity if `mean_intensity < intmin` (strict).                        |
| `detp`       | 0.05    | A probe passes detection at `detP <= detp`; fails at `detP > detp`.                           |
| `ntop`       | 100000  | Number of top-variance probes used by PCA (page 9). MDS (page 3) uses **all** cg/ch non-sex.  |
| `failmin`    | 0.05    | Page-2 CSV cutoff: probes with `fail_rate >= failmin` written to `failed_probes.csv`; dashed vertical line on the plot. Matches meffil's `detectionp.cpgs.threshold`. |
| `inclqual`   | `FALSE` | Page-2 setting: include quality-masked probes in the failure histogram if `TRUE`.             |
| `cores`      | 1       | Cores for SeSAMe (kept low because streaming is per-sample).                                  |
| `savesdf`    | `FALSE` | Persist the full SigDF list to disk (memory-expensive).                                       |
| `collapse`   | `FALSE` | Collapse EPICv2 replicates to prefix.                                                         |
| `collapsemethod` | `"mean"` | `"mean"` or `"minPvalue"`.                                                                |
| `dishref`    | `list(EPIC = "cent12CT.m", "450k" = "cent12CT450k.m")` | EpiDISH reference. Single string OR named-by-platform list. Default = 12-cell Salas 2022. Set `mqcset(dishref = "centDHSbloodDMC.m")` for the legacy 7-cell reference. |
| `dishmethod` | `"RPC"` | EpiDISH method. |

### 5.3 Changing thresholds

```r
# Set globally before any methylQC call
mqcset(samplemin = 0.90, intmin = 1000)
pipeline(indir = "data/idats", outdir = "results")

# Or pass as ... — temporary override for one call
pipeline(indir = "data/idats", outdir = "results", samplemin = 0.90)

# Re-run Stage 2 with looser thresholds, same Stage 1 outputs
mqcset(samplemin = 0.85)
prep(dir = "results")

# Per-call detection threshold
betas_strict <- cleanmat(betas, mask = mask, detP = detP, pthresh = 0.01)

# Inspect / reset
str(mqcopts())
mqcreset()
```

---

## 6. SNP identity verification: `snpcheck()`

Operates on `snp_betas.rds` (samples × rs probes; ~59 on EPIC/EPICv2,
~65 on 450K).

1. NAs replaced with per-probe median (rs probes cluster at 0/0.5/1,
   so median is a sensible imputation).
2. Euclidean distance matrix → classical MDS to 2 dimensions.
3. If a donor column is available, per-donor centroids are computed;
   each sample is assigned to its nearest centroid. Samples whose
   nearest centroid is not their reported donor are flagged.
4. Outputs `snp_mds_coordinates.csv`, `snp_identity_flags.csv` (if
   any), and a two-page PDF (all samples; zoomed view of flagged
   donors).

Signature:

```r
snpcheck(rsbetas,
         ss       = NULL,
         donorcol = NULL,
         outdir   = ".",
         width    = 12, height = 10,
         logger   = NULL)
```

`rsbetas` may be a matrix or a path to `snp_betas.rds`.

---

## 7. EPICv2 replicate collapse and probe lift-over: `liftv2()`

EPICv2 manifests target many CpGs with multiple replicate probes
(`cg`/`ch` probe IDs sharing a stem but differing by suffix). For
cohort analyses every sample must use the **same physical probe** for
a given CpG, otherwise replicate-choice variation enters the analysis
as technical noise. `liftv2()` resolves this cohort-consistently and
then lifts the surviving probe IDs to a legacy platform (EPIC or
HM450) via SeSAMe's `mLiftOver`.

### 7.1 Algorithm

For each CpG with replicate probes:

1. Strip the EPICv2 suffix to build a CpG-level key
   (`sub("_.*$", "", probe_id)`).
2. Compute each probe's cross-sample detection failure rate
   (`rowMeans(detP > pthresh)`); probes with all-NA detection are
   treated as worst.
3. Keep the probe with the lowest failure rate. Ties resolve to the
   first probe by ID order. This is a **cohort-level** decision; it
   selects the same physical probe for every sample, unlike SeSAMe's
   `collapseToPfx(method = "minPvalue")` which picks per-sample.
4. After de-duplication, call `sesame::mLiftOver(mat, platform,
   impute = FALSE)` to lift IDs to the target platform.

`rs` (SNP) probes and control probes are not touched.

### 7.2 Signature

```r
liftv2(mat      = NULL,         # numeric matrix, EPICv2 probes x samples
       detP     = NULL,         # detection p-value matrix (required)
       matpath  = NULL,         # alternative to mat: path to .rds
       detppath = NULL,         # alternative to detP: path to .rds
       platform = c("EPIC", "HM450"),
       pthresh  = NULL,         # default mqcopts()$detp (0.05)
       dedupcsv = NULL,         # optional per-CpG kept/dropped log
       details  = FALSE,
       logger   = NULL)
```

If `details = FALSE` (default), returns the de-duplicated, lifted
matrix. If `details = TRUE`, returns a list with `mat`, `kept_ids`
(surviving EPICv2 IDs), `id_map` (named character vector: EPICv2 ID
→ harmonised ID), and `platform`. The detail list is what you need
to re-align companion matrices (mask, detP) with the internal helper
`apply_epicv2_map(companion, kept_ids, id_map)` so all matrices stay
consistent.

### 7.3 Typical workflow

```r
# Stage 1 has already produced betas_all.rds and detP_all.rds on EPICv2.
lifted <- liftv2(
  matpath  = "results/betas_all.rds",
  detppath = "results/detP_all.rds",
  platform = "EPIC",
  dedupcsv = "results/epicv2_dedup_log.csv",
  details  = TRUE)

# Lift the companion matrices through the same kept_ids + id_map
mask    <- readRDS("results/mask_all.rds")
detP    <- readRDS("results/detP_all.rds")
mask_lifted <- methylQC:::apply_epicv2_map(mask, lifted$kept_ids,
                                           lifted$id_map, fill = FALSE)
detP_lifted <- methylQC:::apply_epicv2_map(detP, lifted$kept_ids,
                                           lifted$id_map, fill = NA_real_)

# Now hand to cleanmat() exactly like any non-v2 cohort
betas_clean <- cleanmat(
  lifted$mat, mask = mask_lifted, detP = detP_lifted,
  probes = "cg", platform = "EPIC")
```

### 7.4 When NOT to use `liftv2()`

- Single-platform analyses where every sample is on EPICv2 and you
  don't need cross-platform comparability. Use SeSAMe's
  `collapseToPfx` (set `mqcset(collapse = TRUE)`) or skip collapse
  entirely.
- Cohorts already collapsed by SeSAMe's per-sample method. `liftv2()`
  expects the matrix to still contain replicate suffixes.

---

## 8. Environment verification: `checkdeps()`

`checkdeps(quiet = FALSE)` is a standalone diagnostic. It does **not**
install or update anything; it verifies that the runtime is correctly
provisioned and prints the exact install command for anything missing
or stale.

What it checks:

- **R version.** Must be ≥ 4.1.0.
- **Required Bioconductor / CRAN packages and minimum versions** —
  baked into the internal `.methylqc_min_versions` list:

  | Package      | Minimum version |
  |--------------|-----------------|
  | sesame       | 1.20.0          |
  | sesameData   | 1.20.0          |
  | EpiDISH      | 2.18.0          |
  | ggplot2      | 3.4.0           |
  | matrixStats  | 0.62.0          |

- **SeSAMe data cache** — calls
  `sesameData::sesameDataGet("genomeInfo.hg38")` as a tracer; failure
  means the cache is not initialised. Fix:
  `sesameData::sesameDataCacheAll()`.
- **EpiDISH reference panel(s)** — checks that every reference name
  in `mqcopts()$dishref` (default: `cent12CT.m`, `cent12CT450k.m`) is
  loadable from the installed
  EpiDISH.

Return value: invisibly, `list(ok = <logical>, problems = <character>)`.
Run on a fresh box before the first pipeline run; the printed report
plus the install commands give you a deterministic provisioning
checklist.

---

## 9. Design rationale

**Raw matrices + separate masks.** Signal and QC are decoupled.
Different downstream analyses need different masks; baking them into
`betas_all.rds` would force one choice on everyone.

**Flag, don't filter.** v1.x wrote exclusion CSVs that looked
authoritative but were advisory; users either ignored them or treated
them as gospel. v2.0.0 produces no exclusion lists by default. The
only sample-level exclusion utility (`flagsamples()`) is opt-in and
explicit in what it does. The only place exclusion actually happens
is `cleanmat()`, which takes IDs the user provides.

**MDS on all probes (change from v1.x top-N).** Top-variable selection
before MDS biases the embedding toward whatever drove that variance
— often a batch effect. Using all complete-case probes gives a more
faithful low-dimensional summary at modest extra cost
(`dist(t(betas))` is the same shape regardless of probe count).

**PCA still uses top-N.** PCA is computationally heavier than MDS,
and top-variable selection is the standard approach for visualising
methylation PCA; we retain it (`ntop = 100000`).

**PC vs associated variable, not PC-vs-PC scatter.** Six PC-vs-PC
scatter plots tell you the shape of variation. PC-vs-variable plots
tell you what each axis means. Six panels on one page also fit a
single sheet (2×3 via `gridExtra::grid.arrange`).

**Full-range probe-failure histogram with a capped Y axis.** Plotting
the per-probe failure distribution at its natural scale is dominated
by the spike near 0% and tells you nothing. Plotting *only* the tail
(an earlier v2 draft) tells you what's already exported to CSV but
hides the shape of the rest of the distribution — specifically,
where the distribution starts to drop. v2.0.0 plots the full 0-1 range
in 40 bins and caps the Y axis at 1.5 × the tallest bin centred at
`fail_rate ≥ 0.05`, deliberately clipping the near-zero spike. The
clipped spike's count is reported in the subtitle. A dashed line at
`failmin` marks the CSV cutoff; the tail probes themselves continue
to be exported to `failed_probes.csv` unchanged.

**`failmin = 0.05` default.** The defensible range across mainstream
EWAS packages is 0% (ChAMP, minfi traditional) through 10% (meffil's
`detectionp.cpgs.threshold = 0.1`), with intermediate choices at 5%
(DNAmArray workflow). 5% is the moderate end of that range and is
appropriate when k-NN imputation is available downstream: at typical
EWAS cohort sizes (n ≳ 100), `cleanmat(..., impute = TRUE)` recovers
the per-probe NAs from the `mask` / `detP` matrices reliably, so the
cost of keeping a probe that fails in 5–9% of samples is modest. The
0% / "drop if any sample fails" rule used historically by ChAMP and
minfi is increasingly seen as too aggressive at EPIC scale (10–20%
of cg probes typically discarded). The Lehne 2015 and Heiss & Just
2019 recommendations are framed in terms of the *per-sample* detection
p-value rather than the cohort fraction; methylQC's per-sample default
remains pOOBAH `p < 0.05` (SeSAMe default), so the cohort fraction is
the only knob protecting against systematically noisy probes — picking
it deliberately matters.

**Intensity plot precedence (MDS > low-intensity > OK).** MDS
outliers warrant the strongest signal regardless of intensity.
Showing a sample as "low intensity" when it is also an MDS outlier
would understate the problem.

**Beta density built from scratch.** v1.x's README claimed this plot
existed but no code produced it. v2.0.0 actually computes it
(`density(x, from = 0, to = 1, n = 512)` per sample on cg/ch probes,
overlaid lines coloured by low-intensity flag).

**Hardcoded Horvath clock and sex probe lists.** Avoids silent
failures when the SeSAMe data cache is unavailable.

**Intensity-based sex inference with regression bands.** Robust to
batch-to-batch shifts in absolute intensity because the bands adapt
to the data; transparently visualised.

**Optional imputation.** Chunked k-NN via
`cleanmat(..., impute = TRUE)`. Chunking is essential at EPIC scale.

**No batch correction.** Analysis-specific. PCA panel (page 9) helps
identify whether batch variables associate with the leading PCs.

**Short user-facing names; internals keep underscores.** A user typing
`cleanmat` repeatedly benefits from one token. Internal helpers
(`resolve_column`, `normalize_sex`, `predict_horvath_age`,
`parse_age_robust`, `geometric_median`, …) keep underscores because
they are referenced by other internals and clarity inside the codebase
matters more than typing speed.

---

## 10. References

1. Zhou W, Triche TJ, Laird PW, Shen H. *SeSAMe: reducing artifactual
   detection of DNA methylation by Infinium BeadChips in genomic
   deletions*. Nucleic Acids Research. 2018;46(20):e123.
2. Triche TJ, Weisenberger DJ, Van Den Berg D, Laird PW, Siegmund KD.
   *Low-level processing of Illumina Infinium DNA Methylation
   BeadArrays*. Nucleic Acids Research. 2013;41(7):e90. (NOOB.)
3. Horvath S. *DNA methylation age of human tissues and cell types*.
   Genome Biology. 2013;14(10):R115.
4. Teschendorff AE, Breeze CE, Zheng SC, Beck S. *A comparison of
   reference-based algorithms for correcting cell-type heterogeneity
   in Epigenome-Wide Association Studies*. BMC Bioinformatics.
   2017;18(1):105. (EpiDISH RPC.)
5. Carrel L, Willard HF. *X-inactivation profile reveals extensive
   variability in X-linked gene expression in females*. Nature.
   2005;434(7031):400–404.
6. Cotton AM, Price EM, Jones MJ, Balaton BP, Kobor MS, Brown CJ.
   *Landscape of DNA methylation on the X chromosome reflects CpG
   density, functional chromatin state and X-chromosome inactivation*.
   Genome Research. 2015;25(8):1091–1099.
7. Guintivano J, Aryee MJ, Kaminsky ZA. *A cell epigenotype specific
   model for the correction of brain cellular heterogeneity bias and
   its application to age, brain region and major depression*.
   Epigenetics. 2013;8(3):290–302.
