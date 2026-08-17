# methylQC

Quality control and preprocessing for Illumina Infinium DNA methylation
BeadChips (HM450, EPIC v1, EPIC v2, MSA, mouse), built on
[sesame](https://bioconductor.org/packages/sesame).

Reads IDATs, preprocesses every array in parallel, flags bad samples, calls sex,
estimates cell composition and epigenetic age, and writes a QC report. Nothing
is deleted — samples and probes are flagged, and you apply your own exclusions
downstream.

See [METHODS.md](METHODS.md) for what the pipeline does and why.

---

## Install

```r
install.packages("pak")
pak::pkg_install("brianchengithub/methylQC")
```

On first use, cache the sesame annotation once:

```r
sesameData::sesameDataCache()
```

**macOS and Linux only** — parallelism uses forking.

---

## 1. Run the pipeline

```r
library(methylQC)

pipeline("~/idats", "~/qcout")
```

One call does everything. It finds the IDATs and your sample sheet, sizes the
run against available memory, processes every array, and writes the report.

```
methylQC complete - ~/qcout

  Platform                  : EPIC, 866,553 probes
  Samples                   : 240 processed, 1 failed
  Stage 1 time              : 18.4 min

  Flagged samples           : 7
      low call rate         : 3
      low intensity         : 2
      sex mismatch          : 2
      MDS outlier           : 0

  Report    : ~/qcout/qc/qc_plots.pdf
  Sheet     : ~/qcout/data/metadata/sample_sheet.csv
  Methods   : ~/qcout/METHODS.txt
  Log       : ~/qcout/logs/pipeline.log
```

Flagged sample IDs and their reasons are printed to the console and the log.

To see the memory plan without committing to a run:

```r
qcplan("~/idats")
```

**Output:**

```
~/qcout/
  METHODS.txt                          what was done, with this run's numbers
  logs/pipeline.log
  data/matrices/
    betas_all.rds                      beta values, UNMASKED, all probes
    detP_all.rds                       detection p-values, all probes
    design_mask.rds                    named logical vector, platform-level
    snp_betas.rds                      rs probe betas
  data/metadata/
    sample_sheet.csv                   one row per sample, every metric and flag
    failed_samples.csv, failed_probes.csv, pc_scores.csv, snp_concordance.csv
  qc/
    qc_plots.pdf                       the report
    qccache.rds                        lets qcplots() redraw without reprocessing
```

---

## 2. Redraw the report at a different threshold

No reprocessing — this reads `qccache.rds` only and takes about a second.

```r
qcplots("~/qcout", detp = 0.01)                  # redraw at a new threshold
qcplots("~/qcout", detp = 0.01, dry = TRUE)      # what would this cost?
qcplots("~/qcout", detp = 0.01, suffix = "d01")  # side by side, nothing overwritten
```

Tunable here: `detp`, `samplemin`, `failmin`, `intmad`, `intfloor`, `inclqual`.
Without `suffix`, `sample_sheet.csv` is updated to match so the report and the
sheet cannot disagree.

Point it at a **parent folder** to do every project beneath it at once:

```r
qcplots("~/projects", detp = 0.01)
```

---

## 3. Apply masks for downstream work

Nothing on disk is masked. The design mask and the detection p-values are
stored separately, so you choose the policy at the point of use:

```r
betas <- readRDS("~/qcout/data/matrices/betas_all.rds")
detp  <- readRDS("~/qcout/data/matrices/detP_all.rds")
dm    <- readRDS("~/qcout/data/matrices/design_mask.rds")

clean <- cleanmat(betas, detp, dm, maskuse = "both")   # design + detection
m     <- mvals(clean)                                   # M-values for an EWAS
```

`maskuse` is `"both"`, `"detection"`, `"design"` or `"none"`. Masked cells
become `NA`; no probe or sample is removed.

Drop the samples that failed QC:

```r
ss   <- read.csv("~/qcout/data/metadata/sample_sheet.csv")
keep <- ss$sample_id[!ss$flagged]
m    <- m[, keep]
```

---

## 4. Harmonise EPIC v2 for clocks and reference panels

EPIC v2 measures some CpGs with several probes, distinguished by a suffix
(`cg00004963_TC21`). Epigenetic clocks and cell-type reference panels key on
bare `cg` identifiers and match nothing until those are resolved.

`pipeline()` does this automatically when it sees the suffixes, so
`betas_all.rds` from an EPIC v2 run is already harmonised. To do it yourself to
a matrix from elsewhere:

```r
cv    <- collapsev2(betas, detp, dm)   # one row per CpG, bare cg identifiers
betas <- cv$betas                       # betas, detp and the mask collapse together
```

The default keeps the probe recommended by Peters et al. (2024). That table is
not shipped; build it once, then reinstall:

```r
pak::pkg_install("bioc::EPICv2manifest")
methylQC::build_epicv2_table()      # from the package source directory
```

Until then it falls back to `"minpval"` with a warning.

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
| `detp` | `0.05` | detection p-value threshold |
| `samplemin` | `0.95` | minimum call rate |
| `intmad` | `3` | MAD multiplier for intensity outliers |
| `intfloor` | `1300` | absolute floor, cohort-level **warning only** |
| `failmin` | `0.05` | probe failure rate above which a probe is listed |
| `mdssd` | `4` | SDs from the MDS centroid before a sample is an outlier |
| `maskuse` | `"both"` | which masks `cleanmat()` applies |
| `savesdf` | `TRUE` | retain SigDFs; withdrawn automatically if too large |
| `extreme` | `FALSE` | never hold a full matrix, for cohorts that will not fit |
| `dishref` / `dishmethod` | `"blood"` / `"RPC"` | EpiDISH settings |

`mqcdefaults()` lists the rest, including the sex-calling and sample-sheet
options.

---

MIT licensed. Issues and pull requests welcome.
