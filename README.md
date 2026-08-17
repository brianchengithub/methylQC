# methylQC

Quality control and preprocessing for Illumina Infinium DNA methylation
BeadChips (HM450, EPIC v1, EPIC v2, MSA, mouse), built on
[sesame](https://bioconductor.org/packages/sesame).

Reads IDATs, preprocesses every array in parallel, flags bad samples, calls sex,
estimates cell composition and epigenetic age, and writes a QC report. Nothing
is deleted — samples and probes are flagged, and you apply your own exclusions
downstream.

Everything a run produces goes into one output directory, organised into
subdirectories: the matrices under `data/matrices/`, the per-sample table and
sidecar CSVs under `data/metadata/`, the report under `qc/`, and the log under
`logs/`. You rarely need to know those paths — the functions below take the
output directory and find what they need.

See [METHODS.md](METHODS.md) for what the pipeline does and why.

---

## Install

```r
install.packages("pak")
pak::pkg_install("brianchengithub/methylQC")
```

That is all — `pak` resolves sesame and the other Bioconductor dependencies,
and the annotation data is downloaded automatically on the first run.

**macOS and Linux only** — parallelism uses forking.

---

## 1. Run the pipeline

```r
library(methylQC)

pipeline(
  indir  = "/path/to/idat/directory",
  outdir = "/path/to/output/directory"
)
```

One call does everything: finds the IDATs and your sample sheet, sizes the run
against available memory, processes every array, and writes the report. It
prints a summary at the end, and names every flagged sample and its reason in
the console and the log.

**What you get:**

```
/path/to/output/directory/
  METHODS.txt                          what was done, with this run's numbers
  logs/
    pipeline.log                       every decision the run made
  data/
    matrices/
      betas_all.rds                    beta values, UNMASKED, all probes
      detP_all.rds                     detection p-values, all probes
      design_mask.rds                  named logical vector, platform-level
      snp_betas.rds                    rs probe betas
    metadata/
      sample_sheet.csv                 one row per sample, every metric and flag
      failed_samples.csv               arrays that could not be processed
      failed_probes.csv                probes failing in many samples
      pc_scores.csv
      snp_concordance.csv
  qc/
    qc_plots.pdf                       the report
    qccache.rds                        lets qcplots() redraw without reprocessing
```

Useful arguments to `pipeline()`:

| Argument | Default | Meaning |
|---|---|---|
| `platform` | inferred | force the array type, e.g. `"EPIC"`, `"EPICv2"`, `"HM450"` |
| `sheet` | auto-detected | path to a sample sheet, if it is not found automatically |
| `workers` | from preflight | forked worker processes |
| `batch` | from preflight | samples per task |
| `extreme` | `FALSE` | never hold a full matrix, for cohorts that will not fit in memory |

---

## 2. Redraw the report at a different threshold

No reprocessing. This reads `qccache.rds` only and takes about a second.

```r
qcplots("/path/to/output/directory", detp = 0.01)
```

Point it at a **parent folder** to do every project beneath it at once:

```r
qcplots("~/projects", detp = 0.01)
```

Everything you can change:

| Argument | Default | Accepts | Meaning |
|---|---|---|---|
| `detp` | `0.05` | 0–1 | detection p-value above which a probe is called failed |
| `samplemin` | `0.95` | 0–1 | call rate below which a sample is flagged |
| `failmin` | `0.05` | 0–1 | probe failure rate above which a probe is listed and marked |
| `intmad` | `3` | 0–20 | MADs below the cohort median before a sample is a low-intensity outlier |
| `intfloor` | `1300` | number or `NA` | absolute intensity floor; a **whole-cohort warning only**, `NA` disables |
| `inclqual` | `FALSE` | `TRUE`/`FALSE` | include design-masked probes in the probe-failure panel |

The sample sheet is rewritten to match, so the report and the CSV cannot
disagree about how many samples failed.

---

## 3. Load masked betas for downstream work

Nothing on disk is masked, so you pick the policy when you use the data. Give
the output directory and say which masks you want — the matrices are found,
read and combined for you:

```r
b <- loadbetas("/path/to/output/directory", maskuse = "both")
```

| `maskuse` | Applies |
|---|---|
| `"both"` (default) | design mask **and** detection p-values |
| `"detection"` | detection p-values only |
| `"design"` | design mask only |
| `"none"` | nothing; the stored matrix as-is |

Masked cells become `NA`; no probe or sample is removed.

Two further arguments, each doing exactly one thing and independent of the
other:

| Argument | Default | Accepts | Does |
|---|---|---|---|
| `values` | `"beta"` | `"beta"`, `"M"` | `"M"` returns `log2(beta / (1 - beta))`; a per-value transform, nothing else changes |
| `samples` | `"all"` | `"all"`, `"passing"` | `"passing"` drops the columns whose `flagged` is `TRUE` in `sample_sheet.csv` |

`flagged` is `TRUE` when any of `low_callrate`, `low_intensity`,
`sex_mismatch` or `mds_outlier` is set. `age_outlier` does not count — it
describes your reported metadata, not the array. For anything finer than
"passing", read the sheet and subset on the individual flag columns:

```r
ss <- read.csv("/path/to/output/directory/data/metadata/sample_sheet.csv")
keep <- ss$sample_id[!ss$low_intensity]        # tolerate a low call rate
b <- loadbetas("/path/to/output/directory")[, keep]
```

So, an EWAS-ready M-value matrix with QC failures removed:

```r
m <- loadbetas("/path/to/output/directory", values = "M", samples = "passing")
```

A **parent folder** returns one matrix per project, named by relative path:

```r
all <- loadbetas("~/projects")
names(all)     # "cohortA/out"  "cohortB/out"  "pilot/out"
```

---

## 4. Harmonise EPIC v2 to EPIC v1

EPIC v2 measures some CpGs with several probes, distinguished by a suffix
(`cg00004963_TC21`). Epigenetic clocks and cell-type reference panels key on
bare `cg` identifiers and match nothing until those are resolved.

`collapsev2()` is the single function that does this. It reduces each replicate
group to one probe — keeping the one Peters et al. (2024) found performs best
against matched EPIC v1 and whole-genome bisulphite data — and renames the rows
to bare `cg` identifiers, so the result is EPIC v1 compatible in both content
and naming. Betas, detection p-values and the design mask are collapsed
together, so they cannot drift apart.

**`pipeline()` does not do this to your stored data.** `betas_all.rds` keeps
EPIC v2's native probe names, because collapsing discards one probe of every
replicate pair and which one to keep is a judgement worth making deliberately.
Stage 2 collapses transiently in memory where the clock and the EpiDISH panels
require bare identifiers, and throws that copy away; the files on disk are
untouched.

Harmonise when you need it:

```r
b  <- loadbetas("/path/to/output/directory", maskuse = "none")
dp <- readRDS("/path/to/output/directory/data/matrices/detP_all.rds")
dm <- readRDS("/path/to/output/directory/data/matrices/design_mask.rds")

cv <- collapsev2(b, dp, dm)
betas <- cv$betas      # one row per CpG, bare cg identifiers
detp  <- cv$detp       # collapsed in the same operation, so they cannot disagree
```

The Peters table is derived from the `EPICv2manifest` annotation package, which
methylQC does not redistribute. Build it once, then reinstall:

```r
pak::pkg_install("bioc::EPICv2manifest")
methylQC::build_epicv2_table()      # from the package source directory
```

Until then `collapsev2()` falls back to keeping the replicate with the lowest
median detection p-value, and says so.

---

MIT licensed. Issues and pull requests welcome.
