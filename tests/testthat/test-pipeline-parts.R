test_that("mqcpath produces the documented layout", {
  d <- "/tmp/out"
  expect_equal(mqcpath(d, "betas"), "/tmp/out/data/matrices/betas_all.rds")
  expect_equal(mqcpath(d, "detp"), "/tmp/out/data/matrices/detP_all.rds")
  expect_equal(mqcpath(d, "designmask"), "/tmp/out/data/matrices/design_mask.rds")
  expect_equal(mqcpath(d, "sheet"), "/tmp/out/data/metadata/sample_sheet.csv")
  expect_equal(mqcpath(d, "cache"), "/tmp/out/qc/qccache.rds")
  expect_equal(mqcpath(d, "methods"), "/tmp/out/METHODS.txt")
  expect_error(mqcpath(d, "nope"), "unknown path key")
})

test_that("a v2 flat layout is rejected with a specific message", {
  d <- file.path(tempdir(), paste0("v2dir", as.integer(runif(1) * 1e6)))
  dir.create(d, recursive = TRUE)
  on.exit(unlink(d, recursive = TRUE))
  saveRDS(1, file.path(d, "betas_all.rds"))
  expect_error(methylQC:::mqccheckdir(d), "flat methylQC v2 layout")
})

test_that("options validate their inputs", {
  on.exit(mqcreset())
  expect_error(mqcset(nosuchoption = 1), "unknown methylQC option")
  expect_error(mqcset(detp = 5), "between 0 and 1")
  expect_error(mqcset(maskuse = "sometimes"), "must be one of")
  expect_error(mqcset(memfrac = 8), "between")
  old <- mqcset(detp = 0.01)
  expect_equal(mqcopts("detp")$detp, 0.01)
  mqcset(detp = old$detp)
  expect_equal(mqcopts("detp")$detp, 0.05)
})

test_that("atomic writes leave the previous file intact when the writer fails", {
  p <- file.path(tempdir(), "atomic.csv")
  on.exit(unlink(p))
  write_csv <- methylQC:::write_csv_atomic
  write_csv(data.frame(a = 1), p)
  expect_true(file.exists(p))
  before <- readLines(p)

  expect_error(methylQC:::write_atomic(NULL, p, function(o, f) stop("boom")))
  expect_equal(readLines(p), before)
  expect_length(list.files(dirname(p), pattern = "atomic\\.csv\\.tmp"), 0L)
})

test_that(".assign_cols is idempotent, unlike a repeated left_join", {
  f <- methylQC:::.assign_cols
  dst <- data.frame(sample_id = c("a", "b"), x = 1:2, stringsAsFactors = FALSE)
  src <- data.frame(sample_id = c("b", "a"), cell_CD4 = c(0.3, 0.7),
                    stringsAsFactors = FALSE)
  once <- f(dst, src, "sample_id")
  twice <- f(once, src, "sample_id")
  expect_equal(names(once), names(twice))          # no .x / .y duplication
  expect_equal(once$cell_CD4, c(0.7, 0.3))         # matched, not positional
  expect_equal(once, twice)
})

test_that("flagsamples never emits a phantom NA row and fails closed", {
  qc <- data.frame(
    sample_id = c("good", "bad", "unmeasurable"),
    call_rate = c(0.99, 0.50, NA_real_),
    mean_intensity_raw = c(4000, 4100, NA_real_),
    low_callrate = c(FALSE, TRUE, NA),
    low_intensity = c(FALSE, FALSE, NA),
    stringsAsFactors = FALSE)
  qc$low_callrate <- is.na(qc$call_rate) | qc$call_rate < 0.95
  qc$low_intensity <- is.na(qc$mean_intensity_raw)

  fl <- flagsamples(qc)
  expect_false(anyNA(fl$flagged))
  expect_true(fl$flagged[fl$sample_id == "bad"])
  expect_true(fl$flagged[fl$sample_id == "unmeasurable"])
  sub <- fl[fl$flagged, , drop = FALSE]
  expect_false(anyNA(sub$sample_id))               # v2 produced an <NA> row here
  expect_equal(nrow(sub), 2L)
})

test_that("normsex and parseage handle the usual encodings", {
  expect_equal(normsex(c("male", "F", "Female", "1", "2", "?")),
               c("M", "F", "F", "M", "F", NA))
  expect_equal(parseage(c("42", "7.5 years", "n/a")), c(42, 7.5, NA))
})

test_that("pcainput excludes sex and design probes and tolerates sparse NAs", {
  set.seed(4)
  np <- 500L; ns <- 10L
  ids <- c(sprintf("cg%04d", 1:490), sprintf("rs%04d", 1:5), sprintf("cgX%03d", 1:5))
  b <- matrix(runif(np * ns), np, ns, dimnames = list(ids, paste0("s", 1:ns)))
  d <- matrix(0.001, np, ns, dimnames = dimnames(b))
  d[1, 1] <- 0.9                              # one failure in 10% of samples
  dm <- stats::setNames(rep(FALSE, np), ids)
  dm[2] <- TRUE
  sexp <- ids[496:500]

  bt <- pcainput(b, d, dm, sexp)
  expect_false(any(grepl("^rs", rownames(bt))))
  expect_false(ids[2] %in% rownames(bt))      # design-masked dropped
  expect_false(any(sexp %in% rownames(bt)))   # sex chromosomes dropped
  expect_false(ids[1] %in% rownames(bt))      # 10% failure exceeds the 1% tolerance
  expect_false(anyNA(bt))
})

test_that("makemat reassembles batch files in order and aligns rows", {
  d <- file.path(tempdir(), paste0("mm", as.integer(runif(1) * 1e6)))
  bd <- mqcpath(d, "batches", create = TRUE)
  on.exit(unlink(d, recursive = TRUE))
  ids <- c("cg1", "cg2", "cg3")
  m1 <- matrix(1:6, 3, 2, dimnames = list(ids, c("s1", "s2")))
  m2 <- matrix(7:9, 3, 1, dimnames = list(rev(ids), "s3"))   # scrambled rows
  saveRDS(m1, file.path(bd, "betas_0001.rds"))
  saveRDS(m2, file.path(bd, "betas_0002.rds"))

  out <- makemat(d, "betas")
  expect_equal(rownames(out), ids)
  expect_equal(colnames(out), c("s1", "s2", "s3"))
  expect_equal(unname(out["cg1", "s3"]), 9L)  # realigned, not taken positionally
})

test_that("methodstext renders without a full run", {
  txt <- methodstext(list(platform = "EPIC", n_samples = 10L,
                          n_probes = 866553L, n_design_masked = 100000L))
  expect_true(any(grepl("^METHODS$", txt)))
  expect_true(any(grepl("ELBAR", txt)))
  expect_true(any(grepl("EPIC", txt)))
  expect_false(any(is.na(txt)))
})
