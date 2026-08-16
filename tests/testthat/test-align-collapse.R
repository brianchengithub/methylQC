test_that("align_to reorders by name and errors on a different probe set", {
  f <- methylQC:::align_to
  v <- c(b = 2, a = 1, c = 3)
  expect_equal(unname(f(v, c("a", "b", "c"), "v")), c(1, 2, 3))
  expect_equal(names(f(v, c("a", "b", "c"), "v")), c("a", "b", "c"))
  expect_error(f(v, c("a", "b", "zzz"), "v"), "missing")
  expect_equal(unname(f(v, c("a", "zzz"), "v", fill = TRUE)), c(1, TRUE))
})

test_that("conform_to catches transposed or reordered matrices by label", {
  f <- methylQC:::conform_to
  ref <- matrix(1:4, 2, 2, dimnames = list(c("p1", "p2"), c("s1", "s2")))
  scrambled <- ref[c(2, 1), c(2, 1), drop = FALSE]
  expect_equal(f(scrambled, ref, "m"), ref)

  ## same dimensions, different labels -- v2's identical(dim()) check passed
  bad <- matrix(1:4, 2, 2, dimnames = list(c("pX", "pY"), c("s1", "s2")))
  expect_error(f(bad, ref, "m"), "row \\(probe\\) sets differ")
})

test_that("stripv2 removes only genuine replicate suffixes", {
  expect_equal(stripv2("cg00004963_TC21"), "cg00004963")
  expect_equal(stripv2("cg00000029"), "cg00000029")
  expect_equal(stripv2("rs9363764"), "rs9363764")
  expect_equal(stripv2("ch.2.1234F"), "ch.2.1234F")
})

test_that("collapsev2 keeps betas, p-values and the design mask in step", {
  ids <- c("cg1_TC21", "cg1_TC22", "cg2_BC11")
  b <- matrix(c(0.1, 0.9, 0.5, 0.2, 0.8, 0.6), 3, 2,
              dimnames = list(ids, c("s1", "s2")))
  d <- matrix(c(0.001, 0.4, 0.001, 0.001, 0.4, 0.001), 3, 2,
              dimnames = dimnames(b))
  dm <- c(cg1_TC21 = FALSE, cg1_TC22 = TRUE, cg2_BC11 = FALSE)

  r <- collapsev2(b, d, dm, method = "minpval")
  expect_equal(rownames(r$betas), c("cg1", "cg2"))
  expect_equal(rownames(r$detp), rownames(r$betas))
  expect_equal(names(r$design), rownames(r$betas))
  ## cg1_TC21 has the lower median p-value, so it wins
  expect_equal(unname(r$betas["cg1", "s1"]), 0.1)
  expect_false(unname(r$design["cg1"]))
})

test_that("collapsev2 mean averages replicates", {
  ids <- c("cg1_TC21", "cg1_TC22")
  b <- matrix(c(0.2, 0.4, 0.6, 0.8), 2, 2, dimnames = list(ids, c("s1", "s2")))
  r <- collapsev2(b, method = "mean")
  expect_equal(unname(r$betas["cg1", "s1"]), 0.3)
  expect_equal(unname(r$betas["cg1", "s2"]), 0.7)
})

test_that("intflag finds planted outliers and stays silent on a clean cohort", {
  set.seed(1)
  x <- rlnorm(200, log(4200), 0.22)
  clean <- intflag(x, nmad = 3, floor = NA)
  expect_lte(sum(clean$flag), 2)          # a tight cohort flags essentially nothing

  x[1:3] <- c(700, 950, 1150)
  dirty <- intflag(x, nmad = 3, floor = NA)
  expect_true(all(dirty$flag[1:3]))
})

test_that("intflag reports a uniformly degraded cohort rather than flagging everyone", {
  set.seed(2)
  x <- rlnorm(200, log(1000), 0.22)
  r <- intflag(x, nmad = 3, floor = 1300)
  expect_true(r$cohort_low)               # the cohort is reported
  expect_lte(sum(r$flag), 5)              # but individuals are not all flagged
})

test_that("callrate_from_hist matches a direct computation", {
  set.seed(3)
  np <- 20000L
  p <- c(runif(np * 0.9, 0, 0.01), runif(np * 0.1, 0.05, 1))
  h <- matrix(methylQC:::phist_one(p), ncol = 1, dimnames = list(NULL, "s1"))
  for (thr in c(0.01, 0.05, 0.2)) {
    exact <- mean(p <= thr)
    cached <- unname(callrate_from_hist(h, thr))
    expect_equal(cached, exact, tolerance = 1 / methylQC:::.MQC_PBINS)
  }
})

test_that("probefail_from_grid returns NULL off the grid", {
  pf <- matrix(c(0L, 10L), 2, 1, dimnames = list(c("cg1", "cg2"), "p0.05"))
  expect_equal(unname(probefail_from_grid(pf, 10, 0.05, 0.05)), c(0, 1))
  expect_null(probefail_from_grid(pf, 10, 0.033, 0.05))
})
