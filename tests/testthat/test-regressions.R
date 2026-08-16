## One test per defect that made a headline feature non-functional in 3.0.0.

test_that("the QC cache carries phist and pfail even when Stage 1 is off disk", {
  ## 3.0.0: .load_stage1() read phist/pfail from a cache that does not exist on
  ## the first prep() run, so build_cache() stored NULLs and every later
  ## qcplots() died with "supply either detp or phist".
  ids <- paste0("cg", 1:50); sid <- paste0("s", 1:4)
  detp <- matrix(runif(200, 0, 0.02), 50, 4, dimnames = list(ids, sid))
  s1 <- list(probe_ids = ids, detp = detp, stats = NULL,
             phist = NULL, pfail = NULL, pgrid = NULL)
  cc <- methylQC:::build_cache(s1, "EPIC")
  expect_false(is.null(cc$phist))
  expect_false(is.null(cc$pfail))
  expect_equal(cc$n_samples, 4L)
  expect_equal(colnames(cc$phist), sid)
  expect_false(is.null(callrate_from_hist(cc$phist, 0.05)))
})

test_that("cached and direct call rates agree, including at the threshold", {
  ## Bin edges were half-open the wrong way, so a p-value exactly equal to the
  ## threshold passed on one path and failed on the other.
  p <- c(rep(0.05, 100), rep(0.2, 100), rep(0.001, 800))
  d <- matrix(p, ncol = 1, dimnames = list(NULL, "s1"))
  h <- matrix(methylQC:::phist_one(p), ncol = 1, dimnames = list(NULL, "s1"))
  for (thr in c(0.001, 0.01, 0.05, 0.2)) {
    expect_equal(unname(callrate_from_hist(h, thr)),
                 unname(callrate(d, thr)), tolerance = 1e-9,
                 info = paste("threshold", thr))
  }
})

test_that("a threshold finer than one histogram bin refuses instead of returning zero", {
  ## 3.0.0 returned a call rate of exactly 0 for every sample and wrote it to
  ## the sample sheet.
  h <- matrix(0L, 1000L, 2L, dimnames = list(NULL, c("a", "b")))
  h[1, ] <- 1000L
  expect_equal(unname(callrate_from_hist(h, 0.01)), c(1, 1))
  expect_null(callrate_from_hist(h, 0.0005))
  expect_error(qcmetrics(data.frame(sample_id = c("a", "b"),
                                    mean_intensity_raw = c(4000, 4100)),
                         phist = h, pthresh = 0.0005),
               "cannot compute call rate")
})

test_that("phist_one treats a missing p-value as a failure", {
  h <- methylQC:::phist_one(c(0.001, NA, 0.001))
  expect_equal(sum(h), 3L)
  expect_equal(h[length(h)], 1L)          # the NA landed in the top bin
})

test_that("the age panel selects the predicted age, not the outlier flag", {
  ## flagsamples() creates `age_outlier` before predictages() adds
  ## `age_horvath353`, and a bare grep("^age_") picked the logical first.
  qc <- data.frame(sample_id = c("a", "b"), low_callrate = FALSE,
                   low_intensity = FALSE, mean_intensity_raw = c(4000, 4100),
                   stringsAsFactors = FALSE)
  fl <- flagsamples(qc)
  fl$age_horvath353 <- c(45, 50)
  ac <- grep("^age_", names(fl), value = TRUE)
  ac <- setdiff(ac, grep("nprobes|nmodel|outlier", ac, value = TRUE))
  ac <- ac[vapply(ac, function(k) is.numeric(fl[[k]]), logical(1))]
  expect_equal(ac, "age_horvath353")
})

test_that("the bundled Horvath clock loads and predicts", {
  ## 3.0.0 fetched the model from two sesameData keys that do not exist, so
  ## epigenetic age was never computed on any human array.
  m <- clockmodel("Horvath353")
  expect_false(is.null(m))
  expect_equal(length(m$probes), 353L)
  expect_equal(length(m$coef), 353L)

  set.seed(6)
  b <- matrix(runif(353 * 3), 353, 3,
              dimnames = list(m$probes, c("a", "b", "c")))
  ages <- predictages(b)
  expect_false(is.null(ages))
  expect_equal(nrow(ages), 3L)
  expect_true(all(is.finite(ages$age_horvath353)))
  expect_equal(ages$age_horvath353_nprobes[1], 353L)
})

test_that("makemat ignores batch files left behind by a previous run", {
  ## 3.0.0 globbed every betas_*.rds, so a shorter re-run into the same output
  ## directory produced a matrix with phantom columns.
  d <- file.path(tempdir(), paste0("mm", as.integer(runif(1) * 1e6)))
  on.exit(unlink(d, recursive = TRUE))
  bd <- mqcpath(d, "batches", create = TRUE)
  dir.create(dirname(mqcpath(d, "sheet")), recursive = TRUE, showWarnings = FALSE)
  ids <- c("cg1", "cg2")

  utils::write.csv(data.frame(sample_id = c("new_s1", "new_s2")),
                   mqcpath(d, "sheet"), row.names = FALSE)
  saveRDS(matrix(1:4, 2, 2, dimnames = list(ids, c("new_s1", "new_s2"))),
          file.path(bd, "betas_0001.rds"))
  saveRDS(matrix(5:8, 2, 2, dimnames = list(ids, c("old_s3", "old_s4"))),
          file.path(bd, "betas_0002.rds"))

  expect_warning(m <- makemat(d, "betas"), "stale")
  expect_equal(colnames(m), c("new_s1", "new_s2"))
})

test_that("mqcset(x = NULL) restores the default instead of deleting the option", {
  on.exit(mqcreset())
  old <- mqcset(workers = 4L)
  expect_equal(mqcopts("workers")$workers, 4L)
  expect_null(old$workers)
  mqcset(workers = old$workers)            # the documented restore idiom
  expect_null(mqcopts("workers")$workers)  # not an "unknown option" error
})

test_that("intfloor accepts NA, as intflag documents", {
  on.exit(mqcreset())
  expect_silent(mqcset(intfloor = NA))
  expect_true(is.na(mqcopts("intfloor")$intfloor))
  expect_error(mqcset(intmad = NA), "must not be missing")
})

test_that("the EpiDISH reference names are all real datasets", {
  skip_if_not_installed("EpiDISH")
  have <- utils::data(package = "EpiDISH")$results[, "Item"]
  expect_true(all(methylQC:::.MQC_DISH_REFS %in% have))
  expect_error(mqcset(dishref = "breast"), "must be one of")
})

test_that("the gc() peak probe reads megabytes, not cell counts", {
  ## 3.0.0 hard-coded gc() column 6, which is "max used" in CELLS when R
  ## inserts a "limit (Mb)" column, and multiplied it by 2^20. Per-sample cost
  ## came out ~6 orders of magnitude high, pinning workers to 1.
  b <- methylQC:::.gc_peak_bytes()
  expect_true(is.na(b) || (b > 2^20 && b < 2^40))
})

test_that("pcainput tolerates a small cohort and one dead sample", {
  ## A fixed 1% tolerance means "zero failures" below n = 100, and one wholly
  ## failed array pushed every probe over it, emptying the PCA input.
  set.seed(7)
  np <- 400L; ns <- 10L
  ids <- sprintf("cg%04d", seq_len(np))
  b <- matrix(runif(np * ns), np, ns, dimnames = list(ids, paste0("s", 1:ns)))
  d <- matrix(0.001, np, ns, dimnames = dimnames(b))
  d[, ns] <- 0.9                              # one sample fails everywhere
  bt <- pcainput(b, d, stats::setNames(rep(FALSE, np), ids))
  expect_gt(nrow(bt), 300L)
  expect_false(anyNA(bt))
})

test_that("snpcheck bounds its output and still surfaces an expected-pair break", {
  set.seed(8)
  np <- 60L; ns <- 12L
  g <- matrix(runif(np * ns), np, ns,
              dimnames = list(paste0("rs", seq_len(np)), paste0("s", seq_len(ns))))
  g[, 2] <- g[, 1]                            # s1 and s2 are the same person
  donors <- stats::setNames(c("d1", "d1", paste0("d", 3:ns)), colnames(g))
  r <- snpcheck(g, donors = donors)
  expect_false(is.null(r))
  expect_lt(nrow(r), ns * (ns - 1) / 2)       # not the full pair list
  hit <- r[r$sample_a == "s1" & r$sample_b == "s2", ]
  expect_equal(nrow(hit), 1L)
  expect_true(hit$same_individual)
})

test_that("betadensity does not disturb the caller's random stream", {
  ## Build the input OUTSIDE the seeded region: constructing it with runif()
  ## inside would consume the stream itself and prove nothing.
  b <- matrix(runif(2000), 500, 4,
              dimnames = list(paste0("cg", 1:500), paste0("s", 1:4)))
  set.seed(99)
  want <- runif(3)
  set.seed(99)
  invisible(methylQC:::betadensity(b))
  expect_equal(runif(3), want)
})
