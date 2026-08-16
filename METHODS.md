# methylQC methods

Version 3.0.0. This document explains what the pipeline does and why. Every
claim about sesame's behaviour below was traced to a specific location in the
sesame source (`github.com/zwdzwd/sesame`, master at the time of writing) and
those locations are cited inline so that they can be re-checked.

`METHODS.txt`, written into every output directory, is the run-specific
version of this document with the real numbers substituted.

---

## 1. Overview

Two stages.

**Stage 1** reads each IDAT pair, applies channel inference, dye bias
correction and noob background correction, computes ELBAR detection
p-values, and writes the beta matrix, the detection p-value matrix and the
probe design mask.

**Stage 2** derives per-sample and per-probe quality metrics, infers sex,
estimates cell composition and epigenetic age, runs the dimensionality
reduction, and writes the QC report and the Methods file.

Preprocessing is entirely within-sample. No step borrows information from
other samples, so a cohort split across several output directories produces
bit-identical preprocessed values to the same samples processed together.
Cohort-relative *judgements* do differ: age outliers, PCA, probe call rates
and the intensity threshold are all computed within whichever set of samples
is present.

---

## 2. The prep chain

sesame's default prep string is `QCDPB`:

| Code | Function | Effect |
|---|---|---|
| `Q` | `qualityMask` | flag probes with design problems |
| `C` | `inferInfiniumIChannel` | re-infer type I probe colour channels |
| `D` | `dyeBiasNL` | non-linear dye bias correction |
| `P` | `pOOBAH` | detection p-values, written to the mask |
| `B` | `noob` | background correction |

methylQC v3 runs, per sample:

```
read IDAT
  -> prepSesame "C"   -> record instrument diagnostics
  -> prepSesame "D"   -> ELBAR(return.pval = TRUE)        [one call]
  -> prepSesame "B"   -> getBetas(mask = FALSE)
```

Four deliberate departures.

### 2.1 ELBAR replaces pOOBAH

Both derive a detection p-value by ranking a probe's reading against a
background population, but they build that population differently.

pOOBAH assumes out-of-band signal *is* background. Infinium type I probes
measure both alleles in one colour channel, so the other channel for those
probes should contain only noise; pOOBAH pools those readings, plus negative
controls, into two channel-specific empirical distributions and takes the
smaller of the two tail fractions (`R/detection.R:163-165`).

ELBAR does not take that on faith. It pools in-band and out-of-band readings,
sorts by total intensity, scans 500 intensity windows for the point at which
the spread of beta values begins to widen — real signal appearing — and
treats everything below that crossover as background, capped at 3000
(`R/detection.R:22-48`). The background population is derived from the data's
own intensity profile.

ELBAR carries two internal guards, both of which appear in
`logs/pipeline.log` when they fire:

- *"Background signal is dichotomous"* (`R/detection.R:37`). The dimmest
  intensity window already shows a beta spread above 0.5, so no background
  population can be located. ELBAR falls back to roughly the ten dimmest
  readings, which makes the resulting p-values very coarsely quantised.
- *"Background signal lacks variation"* (`R/detection.R:51`). The background
  pile's 10th-to-90th percentile spread is under 10 intensity units, so the
  empirical distribution is nearly a step function and almost nothing is
  masked.

ELBAR is substantially more expensive. Benchmarked on synthetic EPIC-scale
data (866,000 probes, single core), pOOBAH's hot path took 0.53 s and ELBAR's
4.37 s, a factor of 8.3. The cost is one line: the window scan runs
`df$beta[df$MU > t1][seq_len(500)]` inside a 500-iteration loop over a
million-row frame, so each iteration performs a full million-element
comparison and subset to keep 500 values. The frame is already sorted by
intensity, so a binary search would make this nearly free, but the
implementation scans linearly.

Reported behaviour: ELBAR masks fewer probes than pOOBAH. In the paper
introducing it (Nucleic Acids Research 2024;52(7):e38) the probe success rate
differed significantly (P = 2.7e-8) while the Spearman correlation difference
did not (P = 0.71), and that comparison was made in low-input datasets
spanning single cell to 250 ng. For standard-input cohorts the advantage may
be smaller. The gain is coverage, not demonstrated accuracy.

### 2.2 Detection runs once, before noob

noob rewrites all four signal columns, `MG`, `UG`, `MR`, `UR`, replacing each
intensity with a background-subtracted estimate plus a flat offset of 15
(`R/background.R:112-118`). Those columns are where the out-of-band
background lives — the out-of-band signal is a view into them, not a separate
object. So noob background-subtracts the background reference itself.

A detection calculation performed after noob therefore compares corrected
foreground against corrected background, which is a different comparison from
the one that produced the mask. sesame states the constraint itself in
`R/open.R:34`: pOOBAH must precede noob because noob modifies the out-of-band
signal.

methylQC v2 ran detection twice: once inside `openSesame` (setting the mask)
and again in `extract_detP()` on the finished, post-noob object. The stored
p-values consequently did not correspond to the stored mask. v3 captures the
single pre-noob call, which fixes the correspondence and removes one full
detection pass per sample.

Because the p-values are retained rather than thresholded, `detP <= detp`
reproduces the detection mask exactly, at whatever threshold is chosen later.

### 2.3 `Q` is dropped

`qualityMask()` writes only to `sdf$mask` (`R/mask.R:176-186`, via `addMask`
at `:15-21`) and touches no intensity. Nothing downstream in this chain reads
that column:

- `inferInfiniumIChannel` never reads the mask; its `mask_failed` argument
  defaults `FALSE` and only writes.
- `dyeBiasNL` defaults to `mask = TRUE`, which — despite the name — means
  *use all Infinium-I probes including masked ones*, with the source comment
  explaining they want the entire support range (`R/dye_bias.R:123-124`). Only
  `mask = FALSE` consults the mask, and `prepSesame` passes no extra
  arguments.
- Its early-exit branch depends on the dye bias statistic, which routes
  through `totalIntensities(sdf)` defaulting to `mask = FALSE`
  (`R/sesame.R:111-113`), so `Q` cannot flip that branch either.
- `ELBAR` calls `signalMU(mask = FALSE)` in both places.
- `noob` filters by `nonuniqMask(platform)` directly (`R/background.R:86`),
  not by the mask column.
- `getBetas(mask = FALSE)` skips masking entirely (`R/sesame.R:209`).

`CDB` and `QCDB` therefore produce identical intensities, identical ELBAR
p-values and identical betas. `Q`'s only effect would be to pre-load the mask
column so that the final mask is a union of design and detection failures with
no record of which — which is exactly what section 3 is designed to avoid.

The multi-mapping protection that sesame's ordering comment attributes to
`Q` is in fact hard-coded inside pOOBAH and noob via `nonuniqMask()`, and is
therefore unaffected.

### 2.4 `I` is not in the string

ELBAR is called directly between `D` and `B` so that its p-values can be
kept. Including `I` as well would compute the same quantity a second time and
discard the result.

---

## 3. Masking

**Nothing on disk is masked.** `betas_all.rds` contains every probe for every
sample, with no `NA` introduced by masking. Two mask types are stored
separately.

**Design mask.** Probes with known design problems: multi-mapping, overlap
with common SNPs, cross-hybridisation. This is a platform lookup, so it is
byte-identical for every sample on that platform, and is stored once as a
named logical vector. For EPIC that is about 3.5 MB, against 4.8 GB for the
probes-by-samples logical matrix that v2 stored at 1500 samples.

**Detection.** Stored as the full matrix of ELBAR p-values, not as a
thresholded mask, so the threshold is a decision made downstream rather than
baked into Stage 1.

Keeping them apart is not tidiness. sesame's `addMask()` only ever sets
`TRUE` and never records a reason, so a fused mask says a probe was masked
but not why. v2's plotting code then had to reverse-engineer the split as
`mask & (detP <= pthresh)`, which misclassifies whenever the mask and the
p-values disagree — and in v2 they always did, for the reason in 2.2.

The practical consequence: several epigenetic clocks and published probe sets
include probes that sit in the design mask. With the two stored separately,
such a probe can be recovered together with its detection status. With them
fused, subtracting the design mask from the union would silently reclassify
any probe that failed *both* checks as design-masked-but-detection-clean,
which is precisely the case that matters.

Masking is applied by `cleanmat()` under `maskuse`, defaulting to `"both"`.

`NA` is treated as failure everywhere. In v2, `mat[detp > pthresh] <- NA`
carried `NA` in the logical index, and R leaves `NA`-indexed elements
untouched, so a probe whose detection could not be computed was silently
treated as having passed.

---

## 4. Quality metrics

### 4.1 Where each statistic is measured

Instrument diagnostics are recorded after `C` and before `D` and `B`, because
they are background and instrument measurements:

- Out-of-band intensity **measures background**. Running it after noob
  measures the background remaining after background was subtracted.
- Dye bias is by definition what dye bias correction removes. Measured after
  `D` it is near-corrected and nearly uninformative.
- Channel switch counts must come after `C`, which produces them.

sesame's own detection statistic group is **not** used, because
`sesameQC_calcStats_detection` calls `pOOBAH()` directly at `R/QC.R:207`.
Under this design that would be a third detection pass per sample and would
contradict ELBAR everywhere else. Detection metrics are computed from the
ELBAR p-values already in hand.

Consequence for existing thresholds: `intmin = 1300` in v2 was calibrated
against intensity measured after the full `QCDPB` chain. The v3 quantity is
different, so both `mean_intensity_raw` and the corrected value are recorded
and the threshold applies to the raw one.

### 4.2 Call rate

```
call_rate = mean(!is.na(detP) & detP <= detp)
```

The denominator is every probe on the array. v2 used
`colMeans(detP <= pthresh, na.rm = TRUE)`, which drops missing probes from
both terms, so a sample with half its p-values missing scored a perfect 1.00 —
and those missing probes are exactly the ones that could not be assessed.

A sample whose call rate cannot be computed is flagged, not passed. In v2 the
resulting `NA` propagated into `out[out$flagged, ]`, which emitted a literal
`<NA>` row into the CSV while never flagging the real sample.

### 4.3 Intensity outliers

Flagged relative to the cohort, not against a fixed cutoff, because absolute
intensity depends on scanner, chemistry, DNA input and array version:

```
flag  when  log2(intensity) < median(log2(intensity)) - intmad * mad(log2(intensity))
```

Log first because intensity is right-skewed; one-sided because only low
intensity indicates failure; three MADs is the conventional line for
"extreme".

A relative rule **cannot detect a uniformly degraded cohort**, because the
median moves with it. Simulated at n = 200:

| Scenario | fixed 1300 | MAD, 3 | Tukey, 3×IQR |
|---|---|---|---|
| healthy cohort, 3 bad arrays | 3, all correct | 3, all correct | 3, all correct |
| healthy cohort, nothing wrong | 0 | 0 | 0 |
| **entire cohort degraded** | **176** | **0** | **0** |

So the absolute floor is kept, but as a **cohort-level warning** rather than a
per-sample flag. Flagging 176 of 200 samples as individual outliers is the
wrong framing; the cohort is the outlier. The floor value should be
calibrated from your own arrays: the first v3 run prints the observed
intensity distribution.

### 4.4 Failures fail closed

Any sample that cannot be processed is recorded in `failed_samples.csv`, kept
in the sample sheet with missing values, given an all-failed detection
profile, and flagged. A sample that could not be assessed should look
maximally suspect, never perfect.

---

## 5. Dimensionality reduction

PCA input is defined **without reference to any user-configurable
threshold**:

- cg and ch probes only
- sex chromosomes excluded
- design-masked probes excluded
- detection failure in at most 1% of samples, at a fixed internal threshold of
  0.05 — deliberately not `detp`
- residual `NA` filled with the probe mean
- top 100,000 by variance

Because nothing tunable reaches it, the PCA is computed once and frozen.
`qcplots()` re-renders its panels from the stored scores and never recomputes
it. That removes roughly 98% of the replot cost and, more importantly, means
the PCA panel and the rest of the report can never be built under different
thresholds without saying so.

v2 required a probe to have passed detection in *every* sample
(`complete.cases`), which makes the probe set depend both on `detp` and on
cohort size. Simulating a realistic per-probe failure distribution, the
fraction of probes with zero failures fell from 90.5% at n = 50 to 74.7% at
n = 1500, whereas a 1% tolerance held near 90% at every cohort size. Since
cell types are processed into separate output directories of unequal size,
the v2 behaviour meant their PCAs were built on systematically different
probe sets.

Storage note: `prcomp`'s `rotation` component is probes by PCs, which at
100,000 probes and 1500 samples is 1.2 GB. The plots need only the scores,
the standard deviations and the probe count, so only those are kept.

---

## 6. EPICv2 replicate probes

EPICv2 measures some CpG sites with more than one probe design, distinguished
by a suffix (`cg00004963_TC21`). Earlier arrays have no suffix.

Collapsing is automatic when suffixes are detected, because clock and
cell-type references key on bare `cg` identifiers and match nothing
otherwise.

**Selection, not averaging.** Probe design affects hybridisation, so
averaging two designs produces a value corresponding to neither EPICv1 nor
whole-genome bisulphite sequencing. Choosing by detection p-value picks the
*brightest* probe, which is a different question from which is most accurate.

The default keeps the probe recommended as superior by Peters et al. (2024,
BMC Genomics 25:251), whose recommendations derive from cross-platform
comparison against matched EPICv1 and WGBS data, measuring sensitivity to
methylation change and measurement precision against a platform consensus.
The table is frozen at package build time so that an annotation update cannot
silently change results for an existing analysis.

Calibration of probe values across array versions is **not** attempted.
Matched-cohort work reports array-level Spearman correlation of 0.968 to
0.981 between EPICv1 and EPICv2, but at the probe level only about 25% of
probes exceed 0.70 and the remainder average below 0.30 — and this holds
regardless of background correction or normalisation, including on raw data.
For most probes a fitted calibration would be regressing noise on noise.

Betas, p-values and the design mask are collapsed together in one operation
keyed on the same group membership, so they cannot disagree afterwards. v2
collapsed only the beta matrix, leaving the mask and detection lookups keyed
on suffixed identifiers that no longer matched; both returned entirely `NA`
and the run applied no quality control at all while producing files of the
correct shape.

---

## 7. Epigenetic age

No mask is applied. The published clocks were fitted on unmasked beta values,
so applying one would feed the model something different from what it was
trained on.

methylQC ships the Horvath (2013) 353-probe coefficients itself, from
Additional File 3 of the original paper. Earlier versions delegated to
`sesame::predictAge()` with a model fetched from sesameData under the keys
`Clock_Horvath353` and `Anno/HM450/Clock_Horvath353.rds`. Neither key exists —
the only clock in sesameData is `MM285.clock347`, for the mouse array — so the
lookup always failed and epigenetic age was never computed on any human array.
Bundling the coefficients also fixes results against annotation updates, which
matters when a run is repeated months later.

Probe coverage is reported per sample. EPICv2 in particular drops probes
several clocks depend on, and an age computed from a reduced probe set should
not look identical to one computed from the full set. A probe present on the
array but missing in a given sample is filled from **that probe's cohort mean**,
not from a constant: substituting a fixed value while retaining the full model
intercept biases every prediction in the same direction.

---

## 7a. Sex calling

Sex is called from sex-chromosome intensity by methylQC's own cohort-relative
caller, not by `sesame::inferSex()`.

Per sample, Stage 1 records the median total (M+U, both channels) intensity
over two curated probe sets: 314 chrY probes excluding pseudo-autosomal and
cross-hybridising probes, and 3,433 X-inactivation-specific probes excluding
pseudo-autosomal. Both are measured **inside the worker**, while the `SigDF` is
still in hand, so the check costs two numbers per sample and does not depend on
whether the preprocessed objects are retained.

Stage 2 then splits the cohort. Candidate cuts are the midpoints between
consecutive sorted log2 chrY values leaving at least three samples either side;
the chosen cut maximises the empty gap between the two groups divided by the
larger within-group standard deviation. Log2 first, because intensity is
multiplicative in scanner gain and DNA input.

The split is accepted only if **both** of the following hold:

- the separation score is at least `sexsep` (default 1.0), and
- the two groups' median chrY intensities differ by at least 1.5×.

Either condition alone is foolable. A gap test alone accepts a chance tail
split in a small single-sex cohort, because the widest spacing between order
statistics of a normal sample can reach a full within-group standard deviation
at small *n*. A ratio test alone accepts a smooth intensity gradient with no
real grouping. Together they separate a genuine male/female split — which
scores about 10 on the gap statistic and 3–4× on the ratio — from a chance one,
by an order of magnitude on at least one axis.

When the split is rejected, **no sex is called for any sample** and no
mismatches are reported. This is the important behaviour. The rule this
replaced chose its cut by minimising within-half regression residuals, a
quantity minimised by cutting anywhere through a homogeneous cloud, so it could
never discover that there was nothing to split. Simulated over 200 replicates
per cohort shape:

| cohort | spurious mismatches, old rule | new rule |
|---|---|---|
| 20F + 0M | 50.0% of the cohort, in 100% of runs | 0 |
| 0F + 20M | 49.8% of the cohort, in 100% of runs | 0 |
| 38F + 2M | 13.2% of the cohort, in 100% of runs | 0 |
| 10F + 10M | 5.4% | 0 |
| 50F + 50M | 0.2% | 0 |

Single-sex cohorts are entirely ordinary — a female breast series, a male
prostate series — so the old rule failed hardest exactly where a genuine sample
swap matters most.

Within an accepted split, each group is a line rather than a blob, because
chrX and chrY intensity both scale with overall array brightness. Orthogonal
distance from the fitted line is used to mark samples that fit neither
cluster; those are reported as **unclear** rather than forced into the nearer
group, which is where aneuploidy, contamination and mixed samples land.

---

## 8. Parallelism

Parallelism is by forking (`BiocParallel::MulticoreParam`), so macOS and
Linux only.

methylQC does **not** call `openSesame()` with a `BPPARAM`. That path has a
specific defect. `readIDATpair` contains (`R/sesame.R:346-350`):

```r
if (is.null(manifest)) {
    df_address <- sesameDataGet(paste0(attr(dm,'platform'), '.address'))
    manifest <- df_address$ordering
    controls <- df_address$controls
}
```

`openSesame()` defaults `manifest = NULL` and passes it through to every
worker, so each forked process independently calls `sesameDataGet()`. That
routes into ExperimentHub and BiocFileCache, which is SQLite, and N forked
processes contending for one SQLite file deadlock — a freeze with no error and
no CPU activity.

methylQC resolves the address object **once in the parent** and passes both
`manifest` and `controls` to every worker. They must go together: `controls`
is assigned only inside that same conditional, so supplying the manifest alone
leaves `attr(sdf, "controls")` unset and noob fails at `negControls()`.

This is not a modification to sesame. It is the manual composition sesame's
own vignette publishes for refined control, using the same exported functions
in the same order.

Other hardening: every worker wraps its own work in `tryCatch`, so a corrupt
IDAT cannot abort a cohort regardless of BiocParallel's error semantics; and
matrix columns are assigned by name with an `identical()` guard rather than
positionally.

---

## 9. Memory

Parent-side requirement is arithmetic, known before any IDAT is read:
`n_probes × n_samples × 8 bytes` per retained double matrix. At 866,553
probes and 1500 samples that is 9.7 GiB each for betas and detection
p-values.

Only the per-sample worker cost needs measuring, and `qcplan()` measures it
by processing one sample — which the run must do anyway — bracketed by
`gc(reset = TRUE)`.

Two v3 changes cut the v2 footprint substantially at no cost: the design mask
became a vector rather than a matrix (4.8 GB saved at 1500 samples), and
M-values are no longer materialised (a further 9.7 GB), since they are a
one-line transform of a matrix already on disk.

Batch size is a memory dial that costs no time. Worker memory is
`workers × batch × 13.2 MiB` per EPIC sample; per-task overhead is
milliseconds against seconds of per-sample work. Batches are kept to at least
four per worker so one slow sample cannot strand the tail.

If resident memory crosses the cap, the partial matrices are written to disk
**before** stopping. A linear model is fitted across completed batches —
resident memory grows linearly in samples because the matrices do — giving
measured bytes per sample, fixed overhead, the projected total, and a
suggested batch size.

`extreme = TRUE` never holds a full matrix: each batch writes its own file,
and `makemat()` reassembles them later.

---

## 10. Known limits

- Parallel execution is unavailable on Windows.
- `collapsemethod = "peters"` requires `build_epicv2_table()` to have been run
  once against the `EPICv2manifest` annotation package; until then it falls
  back to `"minpval"` with a warning.
- Cross-array-version calibration of probe values is not attempted (section 6).
- The absolute intensity floor is a placeholder until calibrated against your
  own arrays (section 4.3).
- Cohort-relative QC judgements are not comparable across output directories
  of different composition, even though the preprocessed values are (section 1).
- Sex calling is cohort-relative and therefore declines on single-sex, very
  small or very noisy cohorts (section 7a). It reports nothing rather than
  guessing; there is no single-sample sex call.
- Sex-chromosome probes are excluded from the PCA input only when
  `GenomeInfoDb` is installed, since it is needed to read chromosome
  assignments from the manifest. It is a `Suggests`, and its absence is
  reported as a warning rather than silently ignored.
- Epigenetic age uses the Horvath (2013) clock only. Numeric output has not
  been validated against real IDATs in this release; the implementation is
  covered by unit tests against the published coefficients.

---

## References

Zhou W, Triche TJ, Laird PW, Shen H. SeSAMe: reducing artifactual detection
of DNA methylation by Infinium BeadChips in genomically heterogeneous
contexts. *Nucleic Acids Research* 2018;46(20):e123.

Zhou W, et al. Low-input and single-cell methods for Infinium DNA methylation
BeadChips. *Nucleic Acids Research* 2024;52(7):e38. [ELBAR]

Triche TJ, Weisenberger DJ, Van Den Berg D, Laird PW, Siegmund KD. Low-level
processing of Illumina Infinium DNA Methylation BeadArrays. *Nucleic Acids
Research* 2013;41(7):e90. [noob]

Peters TJ, et al. Characterisation and reproducibility of the
HumanMethylationEPIC v2.0 BeadChip for DNA methylation profiling.
*BMC Genomics* 2024;25:251.

Lussier AA, et al. Technical variability across the 450K, EPICv1, and EPICv2
DNA methylation arrays. *Clinical Epigenetics* 2024;16:166.

Horvath S. DNA methylation age of human tissues and cell types.
*Genome Biology* 2013;14:R115.

Teschendorff AE, Breeze CE, Zheng SC, Beck S. A comparison of reference-based
algorithms for correcting cell-type heterogeneity in Epigenome-Wide
Association Studies. *BMC Bioinformatics* 2017;18:105. [EpiDISH]
