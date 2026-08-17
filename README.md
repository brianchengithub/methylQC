# methylQC

Quality control and preprocessing for Illumina Infinium DNA methylation
BeadChips (HM450, EPIC v1, EPIC v2, MSA, mouse), built on
[sesame](https://bioconductor.org/packages/sesame).

```r
library(methylQC)

pipeline("~/idats", "~/qcout")
```

That is the whole thing. One call takes a directory of IDATs to a finished QC
report: it sizes the run against the memory on your machine, reads and
preprocesses every array in parallel, computes the metrics, calls sex,
estimates cell composition and epigenetic age, writes a multi-page PDF and a
consolidated sample sheet, and prints a summary of what it found.

```
methylQC complete - ~/qcout

  Platform                  : EPIC, 866,553 probes
  Samples                   : 240 processed, 1 failed
  Stage 1 time              : 18.4 min

  Flagged samples           : 7
      low call rate         : 3
      low intensity         : 2
      sex mismatch          : 2

  Report    : ~/qcout/qc/qc_plots.pdf
  Sheet     : ~/qcout/data/metadata/sample_sheet.csv
  Methods   : ~/qcout/METHODS.txt
  Log       : ~/qcout/logs/pipeline.log
```

---

## What it does

**Preprocessing.** Each IDAT pair is read and passed through channel inference,
non-linear dye bias correction and noob background correction. Detection
p-values come from ELBAR, computed once, after noob — ELBAR locates background
empirically from the corrected signal, and running it earlier leaves the two
colour channels carrying different additive offsets, which makes the
low-intensity beta distribution bimodal and causes ELBAR to fall back to a
ten-probe background that disables detection entirely.

**Nothing on disk is masked.** `betas_all.rds` holds every probe for every
sample. The probe design mask and the detection p-values are stored separately,
and masking is applied as a view when you want it. This matters concretely:
several epigenetic clocks and published probe sets include probes that sit in
sesame's design mask, and keeping the two apart lets you recover such a probe
*and still see whether it detected in that sample*.

**Sample QC.** Call rate against every probe on the array, so a probe whose
p-value could not be computed counts against the sample rather than vanishing
from both sides of the fraction. Intensity outliers are found relative to the
cohort using the median absolute deviation of log2 intensity, with an absolute
floor kept as a whole-cohort warning, because a relative rule cannot notice
that every array in the run is bad.

**Sex calling.** From sex-chromosome intensity, using the package's own
cohort-relative caller rather than an external one. It splits the cohort on
chrY intensity at the point of greatest separation and accepts the split only
when the two groups are genuinely separated. When they are not — a single-sex
cohort, a tiny one, a noisy one — it reports no sex and no mismatches instead
of inventing a split. Samples that fit neither cluster are marked unclear
rather than forced into the nearer one.

**Epigenetic age.** Horvath (2013) 353-probe clock, with the coefficients
bundled in the package so results do not move when an annotation package
updates. Missing probes are filled from the cohort mean for that probe, and
per-sample coverage is recorded.

**Cell composition.** EpiDISH deconvolution against a choice of reference
panels. One failed array does not take the rest of the cohort with it.

**Identity.** rs-probe genotype concordance for every pair, to catch sample
swaps and confirm repeated measurements come from the same donor.

**EPICv2 replicates.** Collapsed automatically when replicate suffixes are
detected, keeping the probe recommended by Peters et al. (2024). Betas,
p-values and the design mask collapse together in one operation, so they cannot
drift apart.

**It tells you what it did.** `METHODS.txt` is written into every output
directory in the register of a journal Methods section, with the run's real
numbers substituted.

---

## Installation

```r
install.packages("pak")
pak::pkg_install("brianchengithub/methylQC")
```

`pak` resolves the Bioconductor dependencies on its own. On first use, cache the
sesame annotation once:

```r
sesameData::sesameDataCache()
```

**macOS and Linux only.** Parallelism uses forking via
`BiocParallel::MulticoreParam`.

---

## Output

```
outdir/
  METHODS.txt                         what was done, with this run's numbers
  logs/pipeline.log
  data/
    matrices/
      betas_all.rds                   beta values, UNMASKED, all probes
      detP_all.rds                    ELBAR detection p-values, all probes
      design_mask.rds                 named logical vector, platform-level
      snp_betas.rds                   rs probe betas
      sdfs_all.rds                    when SigDFs are retained
      batches/                        when extreme = TRUE
    metadata/
      sample_sheet.csv                one row per sample, every metric and flag
      failed_samples.csv
      failed_probes.csv
      pc_scores.csv
      snp_concordance.csv
      run_info.rds
  qc/
    qc_plots.pdf
    qccache.rds                       derived summaries for qcplots()
```

### Sample sheets

The sheet is found in tiers, strictest first: `SampleSheet.csv` next to the
IDATs is matched immediately and nothing else is opened. If that finds nothing
usable the search widens to `SampleSheet*.csv`, then to `samples.*`,
`targets.*`, `pheno.*` and similar, then to any delimited file whose name
contains one of those words — including `.txt`. A candidate only qualifies once
it parses as a table *and* has a column that looks like a sample identifier, so
casting a wide net does not drag in unrelated text files. The delimiter is
inferred from the header rather than assumed, and an Illumina `[Data]` preamble
is skipped.

Columns are then resolved by name and **confirmed against their values**: a
column called `Sex` whose entries are not sexes is rejected in favour of one
whose entries are, and a sheet whose headers no alias list anticipates still
resolves from content alone. The log names every column it chose and how.

Identifiers are matched flexibly too — your sheet may key on the Sentrix
barcode, the IDAT file name, a full path, or `Sentrix_ID` plus
`Sentrix_Position` in separate columns; each is tried and the pairing that
matches the most samples wins.

There is exactly **one** `sample_sheet.csv`. Stage 1 creates it, Stage 2 adds
columns to it, and `qcplots()` keeps its flags in step with the report. It is
never subset and never duplicated, and re-running is idempotent.

No samples and no probes are removed. Everything is flagged, nothing is
deleted; you apply your own exclusion criteria downstream:

```r
betas <- readRDS("~/qcout/data/matrices/betas_all.rds")
detp  <- readRDS("~/qcout/data/matrices/detP_all.rds")
dm    <- readRDS("~/qcout/data/matrices/design_mask.rds")

clean <- cleanmat(betas, detp, dm, maskuse = "both")   # or "detection", "design", "none"
m     <- mvals(clean)
```

M-values are not stored, because they are a one-line transform of a matrix
already on disk.

---

## Retuning a threshold

Changing a threshold does not reprocess anything.

```r
qcplots("~/qcout", detp = 0.01, dry = TRUE)      # what would this cost?
qcplots("~/qcout", detp = 0.01)                  # do it
qcplots("~/qcout", detp = 0.01, suffix = "d01")  # side by side, nothing overwritten
```

PCA, MDS and the beta density panels are frozen — their inputs do not depend on
any tunable threshold — so a retune is typically sub-second. Call rates come
from a cached per-sample p-value histogram and probe failure rates from cached
grid counts; only an off-grid threshold reopens a matrix file. Without
`suffix`, the sample sheet is rewritten to match so the PDF and the CSV cannot
disagree; with `suffix`, nothing canonical is touched.

---

## Memory

`pipeline()` preflights every run: it measures the real per-sample cost on your
machine with your array, reads available RAM, and picks a worker count and
batch size. You can see the plan without committing to it:

```r
qcplan("~/idats")
```

Retaining the preprocessed `SigDF` objects (`savesdf`) is on by default, but a
full set is about 34 MB per EPIC sample, so the preflight withdraws it when the
cohort is large enough that keeping them would not fit — and says so, with the
numbers and what to do about it. Nothing else depends on them, so the run is
otherwise unaffected.

For cohorts that do not fit at all:

```r
pipeline("~/idats", "~/qcout", extreme = TRUE)   # never holds a full matrix
b <- makemat("~/qcout", "betas", write = TRUE)   # reassemble afterwards
prep("~/qcout")                                  # then Stage 2
```

If a run crosses its memory cap it writes partial matrices to disk under
`*_partial.rds` before stopping, then reports measured bytes per sample, fixed
overhead, the projected requirement and a suggested batch size. Partial output
is deliberately not named like complete output, so no later stage mistakes it
for a finished run.

---

## The individual stages

`pipeline()` is `qc()` followed by `prep()`. Each is callable on its own when
you want to re-run one without the others.

| Function | Does |
|---|---|
| `pipeline()` | everything, in one call |
| `qcplan()` | preflight only: memory, workers, batch size, SigDF retention |
| `qc()` | Stage 1: read, preprocess, write matrices (+ Stage 2 unless `stage2 = FALSE`) |
| `prep()` | Stage 2: metrics, sex, age, cell composition, report |
| `qcplots()` | regenerate the threshold-dependent panels from cache |
| `makemat()` | reassemble per-batch matrices from an `extreme = TRUE` run |
| `cleanmat()` | apply a masking policy to a beta matrix |
| `sexcall()` | the cohort sex caller, usable on its own |

---

## Options

```r
mqcopts()                    # current
mqcdefaults()                # factory
mqcset(detp = 0.01)          # change
mqcreset()                   # revert
```

| Option | Default | Meaning |
|---|---|---|
| `detp` | `0.05` | ELBAR detection p-value threshold |
| `maskuse` | `"both"` | which masks `cleanmat()` applies |
| `samplemin` | `0.95` | minimum call rate |
| `intmad` | `3` | MAD multiplier for intensity outliers |
| `intfloor` | `1300` | absolute floor, cohort-level **warning only**; `NA` disables |
| `failmin` | `0.10` | probe failure rate above which a probe is listed |
| `sexsep` | `1.0` | chrY separation floor below which no sex is called |
| `sexband` | `5.0` | distance from a cluster axis, in SDs, beyond which a sample is unclear |
| `sexmin` | `8` | cohort size below which no sex is called |
| `workers` / `batch` | `NULL` | inherited from `qcplan()` |
| `memfrac` | `0.80` | fraction of available RAM the run may use |
| `savesdf` | `TRUE` | retain SigDFs; withdrawn automatically if too large |
| `extreme` | `FALSE` | extreme memory conservation |
| `collapse` | `NA` | EPICv2 replicate handling; `NA` decides automatically |
| `collapsemethod` | `"peters"` | replicate selection rule |
| `dishref` / `dishmethod` | `"blood"` / `"RPC"` | EpiDISH settings |
| `snpmin` | `0.70` | concordance below which an unrelated pair is not reported |
| `sheetpatterns` | 4 patterns | sample sheet filename search, strictest first |
| `idaliases` | see `mqcdefaults()` | accepted sample identifier column names |

---

## EPICv2 replicate table

The default `collapsemethod = "peters"` keeps the probe recommended as superior
by Peters et al. (2024), from a table frozen at build time. Build it once:

```r
pak::pkg_install("bioc::EPICv2manifest")
methylQC::build_epicv2_table()      # from the package source directory
```

Until then it falls back to `"minpval"` with a warning.

---

## Development

`man/` and `NAMESPACE` are generated, not hand-edited:

```r
roxygen2::roxygenise()
```

```sh
R CMD build methylQC && R CMD check methylQC_3.0.2.tar.gz
```

See `METHODS.md` for the methodology in full, and `METHODS.txt` in any output
directory for the run-specific version.

MIT licensed. Issues and pull requests welcome.
