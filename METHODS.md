# methylQC methods

Version 3.0.2. This document explains what the pipeline does and why. Every
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
  -> prepSesame "DB"  -> dye bias correction, then noob
  -> ELBAR(return.pval = TRUE)                            [one call]
  -> getBetas(mask = FALSE)
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

### 2.2 Detection runs once, after noob

The ordering constraint is real, but it points the opposite way for the two
detection methods, and methylQC 3.0.1 got it wrong by carrying a
pOOBAH-specific argument across to ELBAR.

**For pOOBAH, detection must precede noob.** noob rewrites all four signal
columns, `MG`, `UG`, `MR`, `UR`, replacing each intensity with a
background-subtracted estimate plus a flat offset of 15
(`R/background.R:112-118`). Those columns are where the out-of-band background
lives — the out-of-band signal is a view into them, not a separate object — so
noob background-subtracts the background reference itself. pOOBAH *defines*
background as that out-of-band signal, so running it afterwards compares
corrected foreground against corrected background. sesame states the
constraint in `R/open.R:34`.

**For ELBAR, detection must follow noob.** ELBAR does not treat out-of-band
signal as background by definition. It pools in-band and out-of-band readings,
sorts by total intensity, and locates background empirically as the
low-intensity population whose beta values are undifferentiated around 0.5 —
then uses that population as the reference distribution for the p-value. It
re-derives that reference from whatever data it is handed, so running it on
noob-corrected signal is self-consistent.

Running it *before* noob is not. Type II probes measure M in green and U in
red — different channels — while type I probes measure both alleles in the
same channel. Before noob the two channels still carry different additive
background offsets, so dark type II probes are pushed toward beta 0 and 1
while dark type I probes stay near 0.5. The pooled low-intensity distribution
is then bimodal, there is no knee to find, and ELBAR emits `Background signal
is dichotomous`. `dyeBiasNL` rescales the channels against each other but does
not remove the additive offset; noob does.

Measured on `EPIC.1.SigDF`, the reference EPIC v1 array in sesameData:

| | beta spread, 500 dimmest probes | background probes used | distinct p-values | call rate |
|---|---|---|---|---|
| `C -> D -> ELBAR` (3.0.1) | 0.758 | 10 | 8 | 0.9999 |
| `C -> D -> B -> ELBAR` | 0.083 | ~500 | 1108 | 0.9920 |

When the check trips, ELBAR abandons the search and takes its background from
`df$MU[10]` — the ten dimmest probes. An ECDF over ten values leaves eight
distinct p-values for the entire array, every probe at p = 0, and a call rate
of 0.9999 by construction: detection is silently switched off. Comparing the
two orderings, 0.78% of probes (~6,800 on EPIC) are called detected under the
pre-noob ordering that the post-noob ordering correctly masks, and none the
other way.

methylQC v3.0.2 therefore runs `C` → (statistics) → `DB` → ELBAR →
`getBetas`. Per-sample statistics are still measured after `C` and before `D`
and `B`, so the reasoning in section 4.1 is unaffected.

methylQC v2 ran detection twice: once inside `openSesame` (setting the mask)
and again in `extract_detP()` on the finished object, so the stored p-values
did not correspond to the stored mask. v3 captures a single call, which fixes
the correspondence and removes one full detection pass per sample.

Because the p-values are retained rather than thresholded, `detP <= detp`
reproduces the detection mask exactly, at whatever threshold is chosen later.

### 2.3 `Q` is dropped

`Q` (`qualityMask`) writes only to `sdf$mask`. Nothing else in this chain reads
that column, so prepending it changes no number methylQC produces. Rather than
rest that on a reading of sesame's source — whose line numbers have already
drifted once — it is measured, on `EPIC.1.SigDF` under sesame 1.28.1, and
enforced by a regression test:

| | `CDB` vs `QCDB` |
|---|---|
| `MG`, `MR`, `UG`, `UR` after the chain | identical, max abs difference 0 |
| ELBAR detection p-values | identical |
| call rate at p ≤ 0.05 | 0.99204 vs 0.99204 |
| `getBetas(mask = FALSE)` | identical, no `NA` either way |
| `getBetas(mask = TRUE)` | 0 `NA` vs 105,454 `NA` |

`Q` masks 105,454 of 866,553 probes on EPIC, 12.2%. Its entire effect is that
last row: it decides what `getBetas(mask = TRUE)` deletes. methylQC calls
`getBetas(mask = FALSE)` and stores the mask separately as a named logical
vector in `design_mask.rds`, so the information is kept without punching
105,454 holes in the beta matrix.

Applying `Q` would therefore not "clean up" the preprocessing. `dyeBiasNL`
fits on all Infinium-I probes including masked ones, ELBAR calls
`signalMU(mask = FALSE)`, and `noob` filters by `nonuniqMask(platform)`
directly rather than by `sdf$mask` — all three ignore it by construction.

One statistic does move: `sesameQC_calcStats(sdf, "intensity")` respects the
mask, so mean intensity reads 3,153.7 after `QC` against 3,172.9 after `C`, 0.60%
lower. That cannot change which samples are flagged. The design mask is a
platform lookup — `qualityMask(resetMask(sdf))` depends only on the probe set
and the platform, so every sample in a cohort is masked at exactly the same
105,454 probes — and the intensity rule is cohort-relative, comparing each
sample against the cohort's own median absolute deviation. A shift applied
identically to every sample leaves that comparison untouched.

**Where the design mask *is* applied.** Dropping `Q` from the prep string does
not mean the mask is ignored; it means it is applied where it matters, by
methylQC rather than by sesame:

- **PCA and MDS input** exclude design-masked probes outright (section 5).
- **The probe-failure panel** excludes them unless `inclqual = TRUE`.
- **`cleanmat()`** applies them on demand, under the policy the caller picks.

Two places deliberately do *not* mask, because the published methods they
implement were not derived on masked data: **epigenetic age** (section 7) and
**cell composition**. For reference, 8 of the 334 Horvath clock probes present
on EPIC are design-masked (2.4%), as are 15 of the 315 EpiDISH blood reference
probes (4.8%). Dropping those would silently shrink a published model rather
than improve it.

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


### 3.1 Where the design mask is and is not applied

The design mask is stored, never baked in. Section 2.3 shows that prepending
`Q` to the prep string changes nothing about preprocessing, so nothing is lost
by omitting it — but the mask itself is used, in the places where it belongs.
Collected in one table so the accounting is visible:

| step | design mask | why |
|---|---|---|
| dye bias (`D`), noob (`B`), ELBAR | **not applied** | all three ignore `sdf$mask` by construction; applying it changes no number (section 2.3) |
| stored beta matrix | **not applied** | `getBetas(mask = FALSE)`; the mask ships alongside as `design_mask.rds` so no value is destroyed |
| stored detection matrix | **not applied** | p-values are kept for every probe so any threshold can be applied later |
| PCA and MDS input | **applied** | design-problem probes are structured noise and would drive components (section 5) |
| probe-failure panel | **applied**, unless `inclqual = TRUE` | a probe that is masked by design is not evidence about this cohort's detection |
| mean intensity statistic | **not applied** | the shift is uniform across samples and the rule is cohort-relative (section 2.3) |
| sex calling | **not applied** | uses curated chrX/chrY probe lists that were selected independently |
| epigenetic age | **not applied** | 8 of the 334 Horvath probes on EPIC are design-masked; dropping them shrinks a published model (section 7) |
| cell composition | **not applied** | 15 of the 315 EpiDISH blood reference probes are design-masked; same reasoning |
| `cleanmat()` / `loadbetas()` | **caller's choice** | `maskuse` is `"both"`, `"detection"`, `"design"` or `"none"` |

The through-line: the mask is applied wherever methylQC is making a judgement
of its own, and withheld wherever a published model or the stored data is at
stake. Anything withheld is recoverable, because the mask is a separate
object.

## 4. Quality metrics

### 4.1 Where each statistic is measured

Each group is taken at the last point before the step that would destroy it.
The table below was measured on `EPIC.1.SigDF`, the reference EPIC v1 array in
sesameData, rather than reasoned about.

| statistic | measured | why not later | evidence |
|---|---|---|---|
| Infinium-I channel switches | **before `C`** | `C` resolves the disagreement the statistic counts | 65 R→G and 904 G→R before `C`; 0 and 0 after |
| dye bias (`RGratio`, `RGdistort`) | after `C`, before `D` | `D` corrects dye bias by definition | `RGratio` 1.512 → 0.999 across `D` |
| out-of-band intensity | after `C`, before `B` | `noob` subtracts background; these *are* background | `mean_oob_red` 545 → 278 across `B` |
| mean intensity | after `C`, before `D`/`B` | absolute scale shifts under `noob` | see below |
| probe counts by type | after `C` | they should reflect the corrected channel assignment | `num_probes_IR` 92,192 → 93,031 |

The channel-switch placement is a correction. methylQC 3.0.1 measured the whole
group after `C`, on the stated reasoning that "`C` produces them". It does not
— it *consumes* them. `sesameQC_calcStats(sdf, "channel")` calls
`inferInfiniumIChannel(sdf, summary = TRUE)`, which re-infers the channel and
counts probes disagreeing with the declared one. After `C` the declared channel
has already been set to the inferred one, so the off-diagonal counts are zero
for every sample and the diagonal counts reproduce `num_probes_IR` and
`num_probes_IG` exactly. Four of the eighteen recorded columns were therefore
constant or duplicated in every 3.0.1 run.

Mean intensity is the weakest of the four constraints, and the claim made for
it in 3.0.1 was overstated. Scaling an array to 75, 50, 25 and 10 per cent of
its signal is recovered as 0.750, 0.500, 0.250 and 0.100 before `noob` and
0.751, 0.503, 0.254 and 0.105 after, so the cohort-relative MAD rule in section
4.3 would detect a degraded array from either measurement point. What `noob`
changes is the absolute level, by about 10 per cent — and `intfloor` is an
absolute threshold, so the uncorrected scale is the one it can be calibrated
against. That is the reason to measure early, not an appeal to the metric being
otherwise meaningless.

Nothing here interacts with the detection ordering in section 2.2. ELBAR
neither consumes nor produces any of these statistics, and none of them feed
ELBAR, so moving detection after `noob` leaves every measurement point above
unchanged.

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

**"Failed" here means one thing only: the array could not be processed at
all.** `process_one()` wraps the read-and-preprocess of one IDAT pair in a
handler that never throws, so a corrupt file, a truncated scan, an unreadable
manifest or an error inside sesame is caught, recorded, and the cohort
continues. The sample appears in `failed_samples.csv` with the error message,
and it keeps its row in the sample sheet with missing values rather than
vanishing.

This is a different and much rarer thing than *failing QC*. A sample with a
call rate of 0.80 processed perfectly well; it is flagged, not failed.

**What "an all-failed detection profile" means.** A failed sample produced no
p-values at all — there is nothing to put in its column of the detection
matrix. The pipeline fills that gap with the worst possible values instead of
leaving it empty: every probe is recorded as having failed detection, so the
sample's call rate computes to 0 and its column of the beta matrix is `NA`
throughout. Concretely, its p-value histogram is set to all-probes-in-the-
top-bin and its contribution to the per-probe failure counts is incremented
everywhere.

The alternative — leaving the column absent, or `NA`, and letting downstream
code use `na.rm = TRUE` — is how a sample that was never measured ends up with
a call rate of 1.00. That is the defect this replaces: v2 computed
`colMeans(detP <= pthresh, na.rm = TRUE)`, which removes missing probes from
the denominator as well as the numerator, so a sample with half its p-values
missing scored a perfect 1.00 and a sample with all of them missing scored
`NA`, which then propagated into the flagged-samples CSV as a literal `<NA>`
row while never flagging the real sample.

The principle applies wherever a value could not be computed: an unmeasurable
intensity is a flag, an `NA` detection p-value is a failed probe rather than a
passing one, and a `NA` call rate is a failed sample. A sample that could not
be assessed must look maximally suspect, never perfect.

### 4.5 A flag is not an exclusion, and flags are independent

Nothing is removed. `flagged` is a summary column for convenience; every
constituent flag stays in the sample sheet as its own column, and the analyst
decides what to do with each.

This matters because **the flags are computed independently and can legitimately
disagree**. A sample flagged only for call rate has still been fully
preprocessed: its betas are real numbers, its sex-chromosome intensities were
measured, and EpiDISH ran on it like any other. If its sex call matches the
reported sex and its cell composition looks ordinary, that is genuine
information, not an artefact — a modest excess of undetected probes does not
invalidate the ~95% that did detect, and both the sex call and the
deconvolution rest on small, well-behaved probe sets that may be entirely
intact.

So a low call rate is a reason to look, not a verdict. Common readings:

| pattern | likely meaning |
|---|---|
| low call rate, everything else normal | some probe-level dropout; usable for most purposes, and the sex and composition estimates are informative |
| low call rate **and** low intensity | genuinely weak array; treat the whole sample with suspicion |
| low call rate **and** sex mismatch | look for a swap or a mix before blaming the array |
| MDS outlier alone | not necessarily bad quality — check whether it tracks a batch or a tissue on the PC pages before excluding |

The one flag deliberately excluded from `flagged` is `age_outlier`, which
describes the reported metadata distribution rather than array quality;
excluding on it would discard sound arrays from an age-skewed cohort. It is
still reported in `flag_reason` and listed separately in the console output.

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

### 5.3 The principal component pages

PC1 through PC6 are tested against a fixed set of variables, one variable per
page, six panels to a page: chip slide, chip row, chip column, plate, well row,
well column, age and sex. They are drawn whether or not they associate,
because "this component is not the plate" is as useful to a reader as the
converse, and a page that was never drawn is indistinguishable from a variable
that has no effect. Anything not drawn is named in the log with the reason.

Chip position is parsed from the Sentrix barcode by `discover()`
(`200607130026_R06C01` → slide `200607130026`, row `06`, column `01`); well
position is parsed from a well column if the sheet has one, taking the letter
as the row and the digits as the column. A cohort can be laid out on chips, on
plates, or on both, and either layout can carry a batch effect.

Association is eta squared from one-way analysis of variance for a categorical
variable, and Pearson r for a continuous one; the panel type follows, a boxplot
against a scatter. "Continuous" means numeric with more than twelve distinct
values, so a numeric plate number is still treated as a factor.

methylQC's own outputs — `flagged`, `mds_outlier`, `low_callrate` and the rest
— are excluded from the candidate set. Associating a component with a flag
derived from the same beta matrix the PCA is built on is circular, and it
crowds out the batch or tissue variable the page exists to surface. Any
remaining sheet variable reaching eta squared 0.10 against some PC gets its own
page too, capped at six such pages.

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

### 6.1 How the collapse is performed

`collapsev2()` is the single operation. Betas, detection p-values and the
design mask are collapsed together, keyed on the same group membership, so
they cannot disagree afterwards. v2 collapsed only the beta matrix, leaving
the mask and detection lookups keyed on suffixed identifiers that no longer
matched; both returned entirely `NA` and the run applied no quality control at
all while producing files of the correct shape.

Step by step:

1. **Group.** `stripv2()` removes the replicate suffix with the regular
   expression `_[A-Z]{2}[0-9]{2}$`, so `cg00004963_TC21` becomes
   `cg00004963`. The pattern is deliberately narrow: `rs9363764` and
   `ch.2.1234F` are returned unchanged, and the function is therefore safe to
   apply to any platform. Probe identifiers are split on the stripped key.
2. **Select one row per group.** Groups of size one pass through untouched.
   For the rest, `method` decides:
   - `"peters"` (default) looks the group up in the frozen table and keeps the
     recommended probe. If the group is absent from the table, it falls back
     to `"minpval"` for that group alone and records the fallback.
   - `"minpval"` keeps the replicate with the lowest median detection p-value
     across samples — the brightest probe, which answers a different question
     from which is most accurate, hence its status as a fallback.
   - `"mean"` averages the replicates. Available, not recommended: probe
     design affects hybridisation, so an average of two designs corresponds to
     neither.
3. **Apply the same selection to everything.** The chosen row index is used to
   subset the beta matrix, the detection matrix and the design mask alike.
   Under `"mean"`, betas and p-values are averaged and a collapsed probe is
   design-masked only if *every* replicate was.
4. **Rename.** Output row names are the stripped keys, so the matrix leaves
   the function with bare `cg` identifiers — the EPICv1 naming that clocks and
   reference panels expect.
5. **Report.** A count of groups taking each path (`single`, `peters`,
   `minpval_fallback`, `mean`) is returned and written to the log.

### 6.2 The Peters table

The recommendations come from Peters et al. (2024, BMC Genomics 25:251), who
evaluated every EPICv2 replicate group against matched EPICv1 and whole-genome
bisulphite sequencing data and published, per group, which probe performs
best.

methylQC does **not** ship that table, because it is derived from the
`EPICv2manifest` annotation package and redistributing it would fork it.
`build_epicv2_table()` derives it once and writes
`inst/extdata/epicv2_replicates.rds`:

- It restricts the manifest to probes that actually have replicates, using the
  `namerep` column when present and falling back to duplicated `Name` values.
- It locates the recommendation column by pattern rather than by a hard-coded
  name — `sensitivity`, `precision`, `rep_result`, `posrep`, `superior` — since
  the exact name has moved between annotation releases.
- It picks, per group, the probe maximising that column when it is numeric, or
  the first flagged `Y`/`TRUE`/`1` when it is categorical, and falls back to
  manifest order otherwise.

The stored table has one row per replicate group with columns `group` (the
stripped identifier), `chosen` (the full EPICv2 `IlmnID` to keep),
`criterion` (which manifest column was used) and `built` (the date). Freezing
it at build time is the point: an annotation update cannot silently change
results for an analysis that has already been run, and `criterion` records
which release's column the choice came from.

Until the table is built, `collapsemethod = "peters"` falls back to
`"minpval"` with a warning, for the whole matrix rather than silently.

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

## 8. Sex calling

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
scores about 10 on the gap statistic and 3–4× on the ratio — from a chance one.

Once accepted, the **cut-off is placed halfway between the highest
high-confidence female and the next reading above her**, on the log2 scale.
Both terms come from this dataset, so the boundary is dataset-specific and
assumes nothing about the shape of either distribution. A cut-off derived from
the female cluster's MAD does assume normality, and MAD is 0.674 σ for a
normal, so a *k*-MAD rule is really a 0.674*k*-σ rule that a genuine female
crosses from time to time. Measured over 300 replicates, the gap rule gives
0.003 false mismatches per 100-sample cohort against 0.08 for a 3-MAD rule.

"High confidence" trims the female cluster at `sexcutsd` MADs (default 4)
before taking its maximum, and only when the cluster has at least ten members,
below which its own MAD is too noisy to trim on. The trim matters when a female
carries elevated chrY — contamination, a mixed sample, XXY — because she would
otherwise drag the boundary up behind her and let true males through
underneath it. With one such sample planted at chrY 1600 against females at 800
and males at 3000, she is isolated in 99–100% of replicates.

The cut-off is **hard**: every sample with usable sex-chromosome intensity gets
a call, and there is no ambiguous band. Distance from the line, in
female-cluster MADs, is reported as `sex_confidence` so a borderline call is
visible in the sample sheet without being withheld.

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
| 10F + 10M | 5.4% | 0.06 per cohort |
| 50F + 50M | 0.2% | 0.003 per cohort |

Single-sex cohorts are entirely ordinary — a female breast series, a male
prostate series — so the old rule failed hardest exactly where a genuine sample
swap matters most.

---

## 9. The QC report

Eight sample-level panels, then the principal component pages. Every panel is
redrawn from `qc/qccache.rds` by `qcplots()`; only the PCA and MDS
*coordinates* are frozen.

### 9.1 Colour

One palette, fixed a priori, used by every panel that colours samples:

| level | colour | meaning |
|---|---|---|
| OK | steelblue | no sample-level flag |
| Low detection | red3 | call rate below `samplemin` |
| Low intensity | goldenrod3 | more than `intmad` MADs below the cohort median |
| MDS outlier | purple3 | further than `mdssd` SDs from the MDS centroid |

A sample can trip several, so the categories are ordered by precedence: MDS
outlier beats low intensity, which beats low detection. Every scale is drawn
with `drop = FALSE`, so all four levels appear in the legend of every panel
whether or not the cohort contains any — otherwise a clean cohort and a
category the code had forgotten to compute would look identical.

### 9.2 The panels

1. **Detection rate per sample.** Horizontal bars, one per sample, ordered by
   call rate, with a dashed line at `samplemin`. Sample labels are dropped
   above 120 samples, where they stop being legible.
2. **Mean intensity per sample.** Histogram, filled by flag, dashed line at the
   MAD-derived cut-off.
3. **Per-probe sample-failure rate.** Histogram of the fraction of samples in
   which each probe failed, over the full 0–1 range. The y axis is capped at
   1.5× the tallest bar at `fail_rate >= 0.05`, because nearly every probe
   passes in nearly every sample and that leftmost spike is one to two orders
   of magnitude taller than everything else — uncapped it compresses the whole
   informative tail onto the axis. The subtitle states the clipped count and
   the cap when clipping occurs. That 0.05 marks where the spike ends and is
   independent of `failmin`; the two coinciding at the default is a
   coincidence. Design-masked probes are excluded unless `inclqual = TRUE`, and
   probes at or above `failmin` are written to `failed_probes.csv`.
4. **Beta value density.** One curve per sample, flagged samples opaque and in
   their flag colour over a faded cohort.
5. **EpiDISH cell-type proportions.** Stacked composition per sample, labelled
   by sample identifier and ordered by compositional similarity, so alike
   samples sit together and a departure reads as a break in the stack. The
   ordering is the leaf order of an average-linkage hierarchical clustering on
   Euclidean distance; the dendrogram itself is not drawn, since it would take
   a third of the page to say what the ordering already says. Row sums are
   reported in the subtitle rather than forced to 1: RPC does not constrain
   them, and a row summing well below 1 means the reference panel did not
   explain that sample.
6. **Sex chromosome intensity.** chrX against chrY, coloured by *reported* sex,
   with the cut-off drawn. Neither mismatches nor QC failures are marked here:
   the point of the panel is to show the two clusters the call was made from,
   and everything else on it competes with that. Mismatches are in the sample
   sheet, in the log and in the console listing.
7. **Reported versus epigenetic age.** With a transparent ±3 SD band about the
   fitted line.
8. **MDS of retained probes.** With a dashed circle at `mdssd` SDs from the
   cohort centroid.

Then a scree plot, and the principal component pages (section 5.3).

## 10. Parallelism

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

## 11. Memory

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

## 12. Options

Every option, its default and its effect. `mqcopts()` returns the current
values, `mqcset()` changes them, `mqcreset()` restores the defaults. Passing
`NULL` to `mqcset()` restores that one option.

**Detection and masking**

| Option | Default | Effect |
|---|---|---|
| `detp` | `0.05` | ELBAR p-value above which a probe is called failed |
| `maskuse` | `"both"` | which masks `cleanmat()` and `loadbetas()` apply |
| `inclqual` | `FALSE` | include design-masked probes in the probe-failure panel |

**Sample-level thresholds**

| Option | Default | Effect |
|---|---|---|
| `samplemin` | `0.95` | call rate below which a sample is flagged |
| `intmad` | `3` | MADs below the cohort median of log2 intensity before a sample is flagged |
| `intfloor` | `1300` | absolute intensity floor; whole-cohort warning only, `NA` disables |
| `failmin` | `0.05` | probe failure rate above which a probe is listed |
| `mdssd` | `4` | SDs from the MDS centroid before a sample is an outlier |

**Sex calling** (section 8)

| Option | Default | Effect |
|---|---|---|
| `sexsep` | `1.0` | separation floor below which the cohort is declared unimodal and no sex is called |
| `sexcutsd` | `4` | trim width for the high-confidence female set, in MADs; selects who defines the boundary, not the boundary |
| `sexmin` | `8` | cohort size below which no sex is called |

**EPICv2** (section 6)

| Option | Default | Effect |
|---|---|---|
| `collapse` | `NA` | `NA` decides from the probe identifiers; `TRUE`/`FALSE` force it |
| `collapsemethod` | `"peters"` | `"peters"`, `"minpval"` or `"mean"` |

**Cell composition and identity**

| Option | Default | Effect |
|---|---|---|
| `dishref` | `"blood"` | `"blood"`, `"bloodsub"`, `"epithelial"`, `"epifibfat"`, `"12ct"` |
| `dishmethod` | `"RPC"` | `"RPC"`, `"CBS"` or `"CP"` |
| `snpmin` | `0.70` | concordance below which an unrelated sample pair is not reported |

**Compute** (sections 10 and 11)

| Option | Default | Effect |
|---|---|---|
| `workers` | `NULL` | forked processes; `NULL` inherits from the preflight |
| `batch` | `NULL` | samples per task; `NULL` inherits from the preflight |
| `memfrac` | `0.80` | fraction of available RAM the run may use |
| `savesdf` | `TRUE` | retain SigDFs; the preflight withdraws it when they will not fit |
| `extreme` | `FALSE` | never hold a full matrix |

**Sample sheet resolution** (section 4)

`sheetpatterns` is the tiered filename search, strictest first. `idcol`,
`sexcol`, `agecol`, `batchcol`, `donorcol` and `cellcol` name the preferred
column for each role, and the matching `*aliases` vectors list the alternatives
accepted. Each is confirmed against the column's own values before being used,
so an unexpected header still resolves and a misleading one is rejected.

---

## 13. Known limits

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
