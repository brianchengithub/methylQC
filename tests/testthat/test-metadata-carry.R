test_that(".carry_meta attaches sheet columns so flagsamples can see them", {
  ## Regression test: qcmetrics() returns a frame built from the Stage 1
  ## diagnostics, which carries no sample sheet columns. Calling flagsamples()
  ## on it directly left reported_sex as NA for every sample, so sex_mismatch
  ## was silently always FALSE.
  f <- methylQC:::.carry_meta
  qcm <- data.frame(sample_id = c("a", "b"), call_rate = c(0.99, 0.99),
                    mean_intensity_raw = c(4000, 4100),
                    low_callrate = c(FALSE, FALSE),
                    low_intensity = c(FALSE, FALSE),
                    stringsAsFactors = FALSE)
  ss <- data.frame(sample_id = c("b", "a"),
                   Reported_Sex = c("F", "M"),
                   Age = c(30, 40), stringsAsFactors = FALSE)
  cols <- list(id = "sample_id", sex = "Reported_Sex", age = "Age",
               batch = NA_character_, cell = NA_character_, donor = NA_character_)

  out <- f(qcm, ss, cols)
  expect_true("Reported_Sex" %in% names(out))
  expect_equal(out$Reported_Sex, c("M", "F"))   # matched by id, not position

  ## An explicit sex call, so the comparison has something to compare against.
  si <- data.frame(sample_id = c("a", "b"), inferred_sex = c("M", "M"),
                   sex_unclear = c(FALSE, FALSE), sex_confidence = c(3, 3),
                   stringsAsFactors = FALSE)
  attr(si, "called") <- TRUE
  fl <- flagsamples(out, sexcol = "Reported_Sex", agecol = "Age", sexinfo = si)
  expect_equal(fl$reported_sex, c("M", "F"))
  expect_equal(fl$sex_mismatch, c(FALSE, TRUE))
  expect_true(fl$flagged[fl$sample_id == "b"])
})

test_that(".carry_meta does not clobber existing columns or duplicate them", {
  f <- methylQC:::.carry_meta
  qcm <- data.frame(sample_id = "a", Age = 99, stringsAsFactors = FALSE)
  ss <- data.frame(sample_id = "a", Age = 1, stringsAsFactors = FALSE)
  cols <- list(sex = NA_character_, age = "Age", batch = NA_character_,
               cell = NA_character_, donor = NA_character_)
  out <- f(qcm, ss, cols)
  expect_equal(out$Age, 99)                    # existing value preserved
  expect_equal(sum(names(out) == "Age"), 1L)
})

test_that("qcmetrics propagates the intensity thresholds it is given", {
  ## qcplots() passes a modified cfg; before the fix intflag() re-read the
  ## global options and ignored it.
  set.seed(9)
  st <- data.frame(sample_id = paste0("s", 1:60),
                   mean_intensity_raw = c(rlnorm(59, log(4000), 0.15), 400),
                   stringsAsFactors = FALSE)
  h <- matrix(0L, 1000L, 60L, dimnames = list(NULL, st$sample_id))
  h[1, ] <- 1000L                              # every probe detected

  strict <- qcmetrics(st, phist = h, pthresh = 0.05, intmad = 3, intfloor = NA)
  loose  <- qcmetrics(st, phist = h, pthresh = 0.05, intmad = 20, intfloor = NA)
  expect_true(sum(strict$low_intensity) >= 1L)
  expect_equal(sum(loose$low_intensity), 0L)   # a huge multiplier flags nothing
})

test_that("an unmeasurable sample is flagged even when the cohort is tiny", {
  ## intflag() disables the relative rule below n = 5, but a sample with no
  ## measurable intensity is still a failure and must not pass silently.
  r <- intflag(c(4000, NA, 0, 4200), nmad = 3, floor = NA)
  expect_equal(r$flag, c(FALSE, TRUE, TRUE, FALSE))
})
