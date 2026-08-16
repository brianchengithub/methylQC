test_that("cleanmat treats NA detection p-values as failures, not passes", {
  ## v2 bug: `mat[detp > pthresh] <- NA` carries NA in the index, and R leaves
  ## NA-indexed elements untouched, so an unknown p-value silently passed.
  b <- matrix(1:6, 3, 2, dimnames = list(c("cg1", "cg2", "cg3"), c("s1", "s2")))
  storage.mode(b) <- "double"
  d <- matrix(c(0.01, NA, 0.9, 0.01, 0.01, 0.01), 3, 2, dimnames = dimnames(b))

  out <- cleanmat(b, detp = d, maskuse = "detection", pthresh = 0.05)
  expect_true(is.na(out["cg2", "s1"]))   # NA p-value -> masked
  expect_true(is.na(out["cg3", "s1"]))   # p > threshold -> masked
  expect_false(is.na(out["cg1", "s1"]))  # p <= threshold -> kept
})

test_that("an all-NA detection column masks the whole sample", {
  b <- matrix(1.0, 3, 2, dimnames = list(c("a", "b", "c"), c("s1", "s2")))
  d <- matrix(c(rep(NA_real_, 3), rep(0.01, 3)), 3, 2, dimnames = dimnames(b))
  out <- cleanmat(b, detp = d, maskuse = "detection", pthresh = 0.05)
  expect_equal(sum(is.na(out[, "s1"])), 3L)
  expect_equal(sum(is.na(out[, "s2"])), 0L)
})

test_that("maskuse selects the right mask", {
  b <- matrix(1.0, 2, 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  d <- matrix(c(0.9, 0.01, 0.01, 0.01), 2, 2, dimnames = dimnames(b))
  dm <- c(a = FALSE, b = TRUE)

  expect_equal(sum(is.na(cleanmat(b, d, dm, "none"))), 0L)
  expect_equal(sum(is.na(cleanmat(b, d, dm, "detection"))), 1L)
  expect_equal(sum(is.na(cleanmat(b, d, dm, "design"))), 2L)
  expect_equal(sum(is.na(cleanmat(b, d, dm, "both"))), 3L)
})

test_that("call rate counts missing probes against the sample", {
  ## v2 used na.rm = TRUE, so a sample with half its p-values missing scored
  ## a perfect 1.00.
  d <- matrix(c(0.01, 0.01, NA, NA), 2, 2, dimnames = list(NULL, c("good", "half")))
  cr <- callrate(d, 0.05)
  expect_equal(unname(cr["good"]), 1)
  expect_equal(unname(cr["half"]), 0)     # not 1
})

test_that("probefail counts NA as failure", {
  d <- matrix(c(0.01, NA, 0.01, 0.01), 2, 2, dimnames = list(c("p1", "p2"), NULL))
  pf <- probefail(d, 0.05)
  expect_equal(unname(pf["p1"]), 0)
  expect_equal(unname(pf["p2"]), 0.5)
})

test_that("mvals is the logit of beta and stays finite at the boundaries", {
  m <- mvals(c(0, 0.5, 1))
  expect_true(all(is.finite(m)))
  expect_equal(m[2], 0, tolerance = 1e-9)
  expect_lt(m[1], -19)
  expect_gt(m[3], 19)
})
