## The acceptance tests for the sex caller. The rule these enforce is that a
## cohort-relative caller must REFUSE when the cohort gives it nothing to
## separate, rather than splitting a homogeneous cloud and reporting the
## resulting halves as sex mismatches.

## Intensity generators matched to real chrX/chrY behaviour: females carry
## roughly twice the X signal and background-level Y, males the reverse.
femaleX <- function(n) rnorm(n, 6000, 300)
femaleY <- function(n) rnorm(n, 800, 80)
maleX   <- function(n) rnorm(n, 3200, 300)
maleY   <- function(n) rnorm(n, 3000, 250)

test_that("a bimodal cohort is called, and called correctly", {
  set.seed(1)
  n <- 15
  x <- c(femaleX(n), maleX(n)); y <- c(femaleY(n), maleY(n))
  truth <- c(rep("F", n), rep("M", n))
  r <- sexcall(x, y, paste0("s", seq_len(2 * n)))

  expect_true(attr(r, "called"))
  expect_gt(attr(r, "separation"), 1)
  expect_equal(r$inferred_sex, truth)
  expect_false(any(r$sex_unclear))
})

test_that("a single-sex cohort is REFUSED, not split in half", {
  ## The v2 rule scored candidate cuts by within-half regression residuals,
  ## which is minimised by cutting anywhere through one cloud. Simulated over
  ## 200 replicates it labelled 50% of an all-female cohort as sex mismatches,
  ## in 100% of runs. Nothing may be called here.
  set.seed(2)
  for (n in c(10, 20, 40)) {
    rf <- sexcall(femaleX(n), femaleY(n), paste0("f", seq_len(n)))
    expect_false(attr(rf, "called"))
    expect_true(all(is.na(rf$inferred_sex)))

    rm_ <- sexcall(maleX(n), maleY(n), paste0("m", seq_len(n)))
    expect_false(attr(rm_, "called"))
    expect_true(all(is.na(rm_$inferred_sex)))
  }
})

test_that("a single-sex cohort produces zero sex mismatches end to end", {
  ## The failure that mattered: correctly labelled samples reported as swaps.
  set.seed(3)
  n <- 20
  qc <- data.frame(
    sample_id = paste0("s", seq_len(n)),
    low_callrate = FALSE, low_intensity = FALSE,
    mean_intensity_raw = rlnorm(n, log(4000), 0.1),
    Reported_Sex = rep("F", n),
    sex_chrX_intensity = femaleX(n),
    sex_chrY_intensity = femaleY(n),
    stringsAsFactors = FALSE)

  fl <- flagsamples(qc, sexcol = "Reported_Sex")
  expect_equal(sum(fl$sex_mismatch, na.rm = TRUE), 0L)
  expect_true(all(is.na(fl$inferred_sex)))
  expect_false(any(grepl("sex_mismatch", fl$flag_reason)))
})

test_that("a heavily imbalanced cohort does not invent mismatches", {
  ## 38F + 2M scored 13.2% spurious mismatches under the old rule.
  set.seed(4)
  x <- c(femaleX(38), maleX(2)); y <- c(femaleY(38), maleY(2))
  truth <- c(rep("F", 38), rep("M", 2))
  qc <- data.frame(
    sample_id = paste0("s", seq_len(40)),
    low_callrate = FALSE, low_intensity = FALSE,
    mean_intensity_raw = rlnorm(40, log(4000), 0.1),
    Reported_Sex = truth,
    sex_chrX_intensity = x, sex_chrY_intensity = y,
    stringsAsFactors = FALSE)
  fl <- flagsamples(qc, sexcol = "Reported_Sex")
  expect_equal(sum(fl$sex_mismatch, na.rm = TRUE), 0L)
})

test_that("a genuine swap is still detected", {
  set.seed(5)
  n <- 15
  x <- c(femaleX(n), maleX(n)); y <- c(femaleY(n), maleY(n))
  reported <- c(rep("F", n), rep("M", n))
  reported[1] <- "M"                      # one mislabelled female
  qc <- data.frame(
    sample_id = paste0("s", seq_len(2 * n)),
    low_callrate = FALSE, low_intensity = FALSE,
    mean_intensity_raw = rlnorm(2 * n, log(4000), 0.1),
    Reported_Sex = reported,
    sex_chrX_intensity = x, sex_chrY_intensity = y,
    stringsAsFactors = FALSE)
  fl <- flagsamples(qc, sexcol = "Reported_Sex")
  expect_equal(sum(fl$sex_mismatch, na.rm = TRUE), 1L)
  expect_true(fl$sex_mismatch[1])
})

test_that("cohorts below sexmin and all-NA intensities refuse cleanly", {
  r <- sexcall(c(6000, 3200, 6100), c(800, 3000, 810), c("a", "b", "c"))
  expect_false(attr(r, "called"))
  expect_equal(nrow(r), 3L)

  r2 <- sexcall(rep(NA_real_, 20), rep(NA_real_, 20), paste0("s", 1:20))
  expect_false(attr(r2, "called"))
  expect_true(all(is.na(r2$inferred_sex)))
})

test_that("sexintensity strips EPICv2 replicate suffixes before matching", {
  ## The curated lists are bare cg identifiers, and this runs in the worker
  ## before any collapse, so an unstripped match returns nothing on EPICv2.
  ids <- methylQC:::.sesame_chrY_clean[1:20]
  xids <- methylQC:::.sesame_chrX_xlinked[1:20]
  mk <- function(pid) data.frame(
    Probe_ID = pid, MG = 100, MR = 100, UG = 100, UR = 100,
    stringsAsFactors = FALSE)

  plain <- methylQC:::sexintensity(mk(c(ids, xids)))
  expect_equal(plain$chrY, 400)
  expect_equal(plain$chrX, 400)

  v2 <- methylQC:::sexintensity(mk(paste0(c(ids, xids), "_TC21")))
  expect_equal(v2$chrY, 400)
  expect_equal(v2$chrX, 400)
})
