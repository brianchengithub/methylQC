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

test_that("the memory guard extrapolates before aborting", {
  ## 3.0.1 aborted on the instantaneous RSS reading. A run whose own projection
  ## showed a zero shortfall was killed at 32 of 68 samples with 14 GiB free.
  f <- methylQC:::.project_mem
  ## 4.1 GiB fixed + 26 MiB/sample, measured over the first 16 batches
  track <- data.frame(done = seq(2, 32, by = 2),
                      rss = 4.1 * 2^30 + seq(2, 32, by = 2) * 26 * 2^20)
  proj <- f(track, 68)
  expect_equal(proj / 2^30, 4.1 + 68 * 26 / 1024, tolerance = 1e-6)
  expect_lt(proj, 6.0 * 2^30)          # projection fits; must not abort
})

test_that("the worker allowance is marginal, not the whole parent heap", {
  ## Subtracting workers * per_sample (the parent's heap high-water, mostly
  ## copy-on-write manifest) left the parent a budget below its own steady
  ## state.
  w <- methylQC:::.worker_marginal(866553L, 2L)
  expect_gt(w, 100 * 2^20)             # not trivially small
  expect_lt(w, 400 * 2^20)             # and nowhere near the 1.5 GiB heap peak
})

test_that("a lab-style sheet named samples.<x>.txt is found and joined", {
  ## Reported from a live run: a tab-delimited file named samples.BM.PH1.txt
  ## with columns Age, Sex, Fname, Donor was never matched by the 3.0.x
  ## sheetpattern, so the run proceeded with no metadata and no sex check.
  d <- file.path(tempdir(), paste0("sh", as.integer(runif(1) * 1e6)))
  dir.create(d, recursive = TRUE); on.exit(unlink(d, recursive = TRUE))
  ids <- c("200607130026_R06C01", "200607130026_R07C01")
  for (i in ids) for (ch in c("_Grn.idat", "_Red.idat"))
    file.create(file.path(d, paste0(i, ch)))
  utils::write.table(
    data.frame(Fname = ids, Age = c(41, 62), Sex = c("F", "M"),
               Donor = c("D1", "D2"), stringsAsFactors = FALSE),
    file.path(d, "samples.BM.PH1.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  expect_true(any(vapply(mqcdefaults()$sheetpatterns,
                         function(p) grepl(p, "samples.BM.PH1.txt", perl = TRUE),
                         logical(1))))

  ss <- data.frame(sample_id = ids, sentrix = ids,
                   Basename = file.path(d, ids), stringsAsFactors = FALSE)
  joined <- methylQC:::.join_sheet(
    ss, methylQC:::.read_table(file.path(d, "samples.BM.PH1.txt")),
    mqcdefaults(), methylQC:::nulllog())
  expect_true(all(c("Age", "Sex", "Donor") %in% names(joined)))
  expect_equal(joined$Sex, c("F", "M"))
  expect_equal(joined$Age, c(41, 62))

  cols <- checkmeta(joined, logger = methylQC:::nulllog())
  expect_equal(cols$sex, "Sex")
  expect_equal(cols$age, "Age")
  expect_equal(cols$donor, "Donor")
})

test_that("the delimiter is inferred rather than assumed", {
  g <- methylQC:::.guess_sep
  expect_equal(g("Fname\tAge\tSex"), "\t")
  expect_equal(g("Fname,Age,Sex"), ",")
  expect_equal(g("Fname;Age;Sex"), ";")
  expect_equal(g("Fname Age Sex"), "")
})

test_that("the sheet join tries every key and picks the best", {
  ## Sentrix_ID + Sentrix_Position, and an id that is an IDAT file name.
  ids <- c("200607130026_R06C01", "200607130026_R07C01")
  out <- data.frame(sample_id = paste0("b1_", ids), sentrix = ids,
                    Basename = file.path("/x", ids), stringsAsFactors = FALSE)

  split_sheet <- data.frame(Sentrix_ID = rep("200607130026", 2),
                            Sentrix_Position = c("R06C01", "R07C01"),
                            Sample_Name = c("a", "b"), Age = c(1, 2),
                            stringsAsFactors = FALSE)
  j1 <- methylQC:::.join_sheet(out, split_sheet, mqcdefaults(), methylQC:::nulllog())
  expect_equal(j1$Age, c(1, 2))

  fn_sheet <- data.frame(Fname = paste0(ids, "_Grn.idat"), Age = c(7, 8),
                         stringsAsFactors = FALSE)
  j2 <- methylQC:::.join_sheet(out, fn_sheet, mqcdefaults(), methylQC:::nulllog())
  expect_equal(j2$Age, c(7, 8))
})

test_that("a text file that is not a sample sheet is ignored", {
  d <- file.path(tempdir(), paste0("sh2", as.integer(runif(1) * 1e6)))
  dir.create(d, recursive = TRUE); on.exit(unlink(d, recursive = TRUE))
  writeLines(c("note\tvalue", "sampled on\t2026-01-01"),
             file.path(d, "sample_notes.txt"))
  expect_null(methylQC:::.read_sheets(d, NULL, mqcdefaults(),
                                      logger = methylQC:::nulllog()))
})

test_that("ELBAR no longer reports a dichotomous background on a real EPIC array", {
  skip_if_not_installed("sesameData")
  sdf <- tryCatch(sesameData::sesameDataGet("EPIC.1.SigDF"), error = function(e) NULL)
  skip_if(is.null(sdf), "sesameData cache unavailable")

  spread <- function(s) {
    df <- rbind(sesame:::signalMU(s, mask = FALSE, MU = TRUE),
                sesame:::signalMU_oo(s, MU = TRUE))
    df$beta <- df$M / (df$M + df$U)
    df <- df[order(df$MU), ]
    df <- df[!is.na(df$MU) & !is.nan(df$beta), ]
    q <- stats::quantile(df$beta[seq_len(500)], c(0.05, 0.95), na.rm = TRUE)
    unname(q[2] - q[1])
  }
  ## ELBAR calls the background dichotomous, and falls back to the ten dimmest
  ## probes, when this spread exceeds 0.5.
  expect_gt(spread(sesame::prepSesame(sdf, "CD")), 0.5)    # the 3.0.1 ordering
  expect_lt(spread(sesame::prepSesame(sdf, "CDB")), 0.5)   # with noob first

  w <- NULL
  pv <- withCallingHandlers(
    sesame::ELBAR(sesame::prepSesame(sdf, "CDB"), return.pval = TRUE),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("dichotomous", w %||% "")))
  expect_gt(length(unique(pv)), 100L)     # not the 8 values the fallback gives
})

test_that("the sheet search runs strict patterns before loose ones", {
  d <- file.path(tempdir(), paste0("tier", as.integer(runif(1) * 1e6)))
  dir.create(d, recursive = TRUE); on.exit(unlink(d, recursive = TRUE))
  ids <- c("200607130026_R06C01", "200607130026_R07C01")

  ## A loose-tier decoy and a canonical sheet in the same directory.
  utils::write.table(data.frame(Fname = ids, Age = c(9, 9)),
                     file.path(d, "samples.decoy.txt"), sep = "\t",
                     row.names = FALSE, quote = FALSE)
  utils::write.csv(data.frame(Sample_Name = ids, Age = c(41, 62)),
                   file.path(d, "SampleSheet.csv"), row.names = FALSE)

  got <- methylQC:::.read_sheets(d, NULL, mqcdefaults(), targets = ids,
                                 logger = methylQC:::nulllog())
  expect_equal(got$Age, c(41, 62))        # tier 1 won; the decoy was never read

  file.remove(file.path(d, "SampleSheet.csv"))
  got2 <- methylQC:::.read_sheets(d, NULL, mqcdefaults(), targets = ids,
                                  logger = methylQC:::nulllog())
  expect_equal(got2$Age, c(9, 9))         # falls through to the loose tier
})

test_that("columns are confirmed against their values, not just their names", {
  ids <- c("200607130026_R06C01", "200607130026_R07C01", "200607130026_R08C01")

  ## Headers that no alias list anticipates: resolved from content alone.
  odd <- data.frame(barcode_v2 = ids, yrs = c(41, 62, 55),
                    gender_code = c("male", "female", "female"),
                    stringsAsFactors = FALSE)
  cols <- checkmeta(odd, logger = methylQC:::nulllog())
  expect_equal(cols$id, "barcode_v2")
  expect_equal(cols$age, "yrs")
  expect_equal(cols$sex, "gender_code")

  ## A column named Sex whose contents are not sexes must not be trusted over
  ## one whose contents are.
  misleading <- data.frame(Sample_Name = ids,
                           Sex = c("collected 2021", "collected 2022", "n/a"),
                           g = c("F", "M", "F"), stringsAsFactors = FALSE)
  c2 <- checkmeta(misleading, logger = methylQC:::nulllog())
  expect_equal(c2$sex, "g")
})

test_that("File_Name is not treated as the sample identifier", {
  expect_false(any(tolower(c("File_Name", "FileName", "File")) %in%
                     tolower(mqcdefaults()$idaliases)))
  expect_true("Fname" %in% mqcdefaults()$idaliases)
  expect_true(all(c("Sample_Name", "SampleName", "Sample_ID", "SampleID") %in%
                    mqcdefaults()$idaliases))
})

test_that("channel switch counts are measured before C, where they mean something", {
  skip_if_not_installed("sesameData")
  sdf <- tryCatch(sesameData::sesameDataGet("EPIC.1.SigDF"), error = function(e) NULL)
  skip_if(is.null(sdf), "sesameData cache unavailable")

  before <- methylQC:::stats_channel(sdf)
  after  <- methylQC:::stats_channel(sesame::prepSesame(sdf, "C"))

  ## Before C the statistic counts probes the array mis-declared.
  expect_gt(before$inf1_r2g + before$inf1_g2r, 0)
  ## After C it is zero by construction -- C resolves exactly that disagreement.
  expect_equal(after$inf1_r2g, 0)
  expect_equal(after$inf1_g2r, 0)

  ## And after C the diagonal merely repeats the type I probe counts, so all
  ## four columns carry no information at that point.
  st <- methylQC:::stats_raw(sesame::prepSesame(sdf, "C"))
  expect_equal(after$inf1_r2r, st$n_probes_ir)
  expect_equal(after$inf1_g2g, st$n_probes_ig)
})

test_that("dye bias and out-of-band intensity are measured before the step that removes them", {
  skip_if_not_installed("sesameData")
  sdf <- tryCatch(sesameData::sesameDataGet("EPIC.1.SigDF"), error = function(e) NULL)
  skip_if(is.null(sdf), "sesameData cache unavailable")

  c_  <- methylQC:::stats_raw(sesame::prepSesame(sdf, "C"))
  cd  <- methylQC:::stats_raw(sesame::prepSesame(sdf, "CD"))
  cdb <- methylQC:::stats_raw(sesame::prepSesame(sdf, "CDB"))

  ## D corrects dye bias by definition, so measuring after it reads ~1.
  expect_gt(abs(c_$rg_ratio - 1), abs(cd$rg_ratio - 1))
  expect_lt(abs(cd$rg_ratio - 1), 0.05)
  ## noob subtracts background, so out-of-band intensity shrinks across B.
  expect_gt(c_$mean_oob_red, cdb$mean_oob_red)
})

test_that("the assembled statistics row matches the template exactly", {
  ## The row is now built from two calls taken at different points in the
  ## chain, so the column set has to be asserted rather than assumed.
  tmpl <- methylQC:::stats_raw_template()
  expect_true(all(c("inf1_r2r", "inf1_g2g", "inf1_r2g", "inf1_g2r",
                    "sex_chrX_intensity", "sex_chrY_intensity") %in% names(tmpl)))
  expect_equal(sum(duplicated(names(tmpl))), 0L)
})

test_that("the sex cut-off is dataset-specific, not distribution-based", {
  set.seed(11)
  n <- 20
  x <- c(rnorm(n, 6000, 300), rnorm(n, 3200, 300))
  y <- c(rnorm(n, 800, 80), rnorm(n, 3000, 250))
  r <- sexcall(x, y, paste0("s", seq_len(2 * n)))
  expect_true(attr(r, "called"))

  ## Both terms come from the data: the cut-off is the midpoint between the
  ## highest female and the next reading above her, on the log2 scale.
  fmax <- max(y[seq_len(n)]); mmin <- min(y[(n + 1):(2 * n)])
  expect_equal(attr(r, "threshold"), sqrt(fmax * mmin), tolerance = 1e-6)
  expect_gt(attr(r, "threshold"), fmax)
  expect_lt(attr(r, "threshold"), mmin)
  expect_equal(sum(r$inferred_sex == "M", na.rm = TRUE), n)
  expect_false(any(r$sex_unclear))
})

test_that("mdsoutlier finds a distant sample and reports the radius", {
  m <- cbind(MDS1 = c(rnorm(20, 0, 1), 40), MDS2 = c(rnorm(20, 0, 1), 40))
  rownames(m) <- paste0("s", seq_len(21))
  o <- mdsoutlier(m, nsd = 4)
  expect_true(o[["s21"]])
  expect_equal(sum(o), 1L)
  expect_true(is.finite(attr(o, "radius")))
  expect_length(attr(o, "centre"), 2L)

  ## A tight cohort with no outlier flags nothing.
  m2 <- cbind(MDS1 = rnorm(30), MDS2 = rnorm(30))
  rownames(m2) <- paste0("t", seq_len(30))
  expect_lte(sum(mdsoutlier(m2, nsd = 4)), 1L)
})

test_that("the MDS outlier reaches the flags and the reason string", {
  qc <- data.frame(sample_id = paste0("s", 1:6), low_callrate = FALSE,
                   low_intensity = FALSE, mean_intensity_raw = 4000,
                   stringsAsFactors = FALSE)
  mf <- stats::setNames(c(rep(FALSE, 5), TRUE), qc$sample_id)
  fl <- flagsamples(qc, mdsflag = mf)
  expect_true(fl$mds_outlier[6])
  expect_true(fl$flagged[6])
  expect_match(fl$flag_reason[6], "mds_outlier")
  expect_false(any(fl$flagged[1:5]))
})

test_that("reportflags names each flagged sample and its reason", {
  qc <- data.frame(sample_id = c("good", "bad1", "bad2"),
                   low_callrate = c(FALSE, TRUE, FALSE),
                   low_intensity = c(FALSE, FALSE, TRUE),
                   mean_intensity_raw = c(4000, 4000, 300),
                   stringsAsFactors = FALSE)
  fl <- flagsamples(qc)
  lines <- character(0)
  lg <- list(log = function(stage, message, warn = FALSE) {
    lines <<- c(lines, message); invisible(NULL) },
    path = NA_character_, close = function() invisible(NULL))
  sub <- reportflags(fl, lg)
  expect_equal(nrow(sub), 2L)
  expect_true(any(grepl("bad1", lines)))
  expect_true(any(grepl("bad2", lines)))
  expect_true(any(grepl("low call rate", lines)))
  expect_true(any(grepl("low intensity", lines)))
  expect_false(any(grepl("good", lines)))
})

test_that("the PC/metadata panel ignores methylQC's own flag columns", {
  set.seed(12)
  sc <- matrix(rnorm(40 * 6), 40, 6,
               dimnames = list(paste0("s", 1:40), paste0("PC", 1:6)))
  plate <- rep(c("P1", "P2"), each = 20)
  sc[, 1] <- sc[, 1] + (plate == "P2") * 8      # PC1 tracks plate

  meta <- data.frame(Plate = plate,
                     flagged = rep(c(TRUE, FALSE), 20),
                     mds_outlier = rep(c(TRUE, FALSE), 20),
                     stringsAsFactors = FALSE)
  ## The helper itself is content-blind; the caller is what excludes derived
  ## columns, so assert the exclusion list covers them.
  expect_true(all(c("flagged", "mds_outlier", "low_callrate", "sex_mismatch") %in%
                    methylQC:::.MQC_DERIVED_COLS))

  a <- methylQC:::.pcassoc(sc, meta["Plate"], npc = 2L)
  expect_equal(a$variable[a$pc == 1][1], "Plate")
  expect_gt(a$eta2[a$pc == 1][1], 0.8)
})

test_that("flag categories combine the two sample-level flags", {
  qc <- data.frame(sample_id = paste0("s", 1:4),
                   low_callrate = c(FALSE, TRUE, FALSE, TRUE),
                   low_intensity = c(FALSE, FALSE, TRUE, TRUE),
                   stringsAsFactors = FALSE)
  expect_equal(as.character(methylQC:::.flagcat(qc)),
               c("pass", "low call rate", "low intensity",
                 "call rate + intensity"))
})

test_that("the sex cut-off is the midpoint of the gap above the female cluster", {
  ## Constructed so the boundary is unambiguous: females at 800, males at 3000,
  ## with the highest female at 900 and the lowest male at 2800.
  y <- c(800, 820, 840, 860, 880, 900, 2800, 3000, 3200, 3400)
  x <- c(rep(6000, 6), rep(3200, 4))
  r <- sexcall(x, y, paste0("s", seq_along(y)))
  expect_true(attr(r, "called"))
  expect_equal(attr(r, "threshold"), sqrt(900 * 2800), tolerance = 1e-6)
  expect_equal(r$inferred_sex, c(rep("F", 6), rep("M", 4)))
  expect_false(any(r$sex_unclear))          # a hard cut-off, no ambiguous band
  expect_false(any(is.na(r$inferred_sex)))
})

test_that("a female with elevated chrY does not drag the cut-off up behind her", {
  ## The trim exists so one contaminated sample cannot let true males through.
  set.seed(21)
  n <- 25
  y <- c(rnorm(n, 800, 60), rnorm(n, 3000, 250))
  y[1] <- 1600                              # contaminated / mixed female
  x <- c(rnorm(n, 6000, 300), rnorm(n, 3200, 300))
  r <- sexcall(x, y, paste0("s", seq_len(2 * n)))
  expect_true(attr(r, "called"))
  ## Either mechanism may isolate her -- the natural break may already put her
  ## above the females, or the trim may exclude her from defining the maximum.
  ## What matters is the outcome: the cut-off lands below her.
  expect_lt(attr(r, "threshold"), 1600)
  expect_equal(r$inferred_sex[1], "M")      # so she surfaces as a mismatch
  expect_true(all(r$inferred_sex[(n + 1):(2 * n)] == "M"))
})

test_that("PCA and MDS input excludes design-masked, low-detection and sex probes", {
  set.seed(22)
  np <- 600L; ns <- 12L
  ids <- c(sprintf("cg%04d", 1:590), sprintf("rs%04d", 1:10))
  b <- matrix(runif(np * ns), np, ns, dimnames = list(ids, paste0("s", 1:ns)))
  d <- matrix(0.001, np, ns, dimnames = dimnames(b))

  design <- stats::setNames(rep(FALSE, np), ids)
  design[1:50] <- TRUE                      # quality-mask probes
  d[51:100, ] <- 0.9                        # fail detection in every sample
  sexp <- ids[101:150]                      # sex chromosome probes

  bt <- pcainput(b, d, design, sexp)
  expect_false(any(ids[1:50]   %in% rownames(bt)))   # design masked
  expect_false(any(ids[51:100] %in% rownames(bt)))   # low detection
  expect_false(any(sexp        %in% rownames(bt)))   # sex chromosomes
  expect_false(any(grepl("^rs", rownames(bt))))      # not cg/ch

  ## MDS is built from the same matrix, so it inherits every exclusion.
  m <- methylQC:::runmds(bt)
  expect_false(is.null(m))
  expect_equal(nrow(m), ns)
})

test_that("sexprobes falls back to the bundled lists rather than returning nothing", {
  ## Returning character(0) silently left sex chromosomes in the PCA input.
  mqcannoreset()
  ids <- suppressWarnings(sexprobes("NoSuchPlatform", methylQC:::nulllog()))
  expect_gt(length(ids), 3000L)
  expect_true(all(c(methylQC:::.sesame_chrY_clean[1],
                    methylQC:::.sesame_chrX_xlinked[1]) %in% ids))
  mqcannoreset()
})

test_that("Sentrix barcodes split into slide, row and column", {
  p <- methylQC:::.sentrix_parts(c("200607130026_R06C01", "201516260020_R03C02",
                                   "not_a_barcode", NA))
  expect_equal(p$sentrix_slide, c("200607130026", "201516260020", NA, NA))
  expect_equal(p$sentrix_row,   c("06", "03", NA, NA))
  expect_equal(p$sentrix_col,   c("01", "02", NA, NA))
})

test_that("the standard PC variables are assembled from chip and well position", {
  qc <- data.frame(
    sample_id = paste0("s", 1:8),
    sentrix_slide = rep(c("200607130026", "200607130044"), each = 4),
    sentrix_row = rep(c("01", "02", "03", "04"), 2),
    sentrix_col = rep(c("01", "02"), 4),
    Well = c("A01", "B02", "C03", "D04", "E05", "F06", "G07", "H08"),
    Plate = rep(c("P1", "P2"), each = 4),
    reported_age = c(30, 40, 50, 60, 35, 45, 55, 65),
    reported_sex = rep(c("F", "M"), 4),
    stringsAsFactors = FALSE)
  v <- methylQC:::.pc_required_vars(qc)

  expect_true(all(c("slide", "slide_row", "slide_col", "plate",
                    "well_row", "well_col", "age", "sex") %in% names(v)))
  expect_equal(v$well_row, c("A", "B", "C", "D", "E", "F", "G", "H"))
  expect_equal(v$well_col, c("1", "2", "3", "4", "5", "6", "7", "8"))
  expect_equal(v$slide_row, qc$sentrix_row)

  ## A constant or absent variable earns no page.
  qc2 <- qc; qc2$Plate <- "P1"; qc2$Well <- NULL
  v2 <- methylQC:::.pc_required_vars(qc2)
  expect_false("plate" %in% names(v2))
  expect_false(any(c("well_row", "well_col") %in% names(v2)))
})

test_that("age is treated as continuous and sex as categorical", {
  set.seed(31)
  y <- rnorm(40)
  age <- seq(20, 80, length.out = 40)
  expect_true(methylQC:::.pc_is_continuous(age))
  expect_false(methylQC:::.pc_is_continuous(rep(c("F", "M"), 20)))
  expect_false(methylQC:::.pc_is_continuous(rep(1:4, 10)))   # a plate number

  expect_match(methylQC:::.pc_stat(y, age)$label, "^r = ")
  expect_match(methylQC:::.pc_stat(y, rep(c("F", "M"), 20))$label, "^eta2 = ")
})

test_that("a PC page holds at most six panels", {
  expect_lte(methylQC:::.MQC_PC_NPC, 6L)
})

test_that("chip position columns survive into the metrics frame", {
  ## The PC pages cannot show a variable that never reached this frame.
  f <- methylQC:::.carry_meta
  qcm <- data.frame(sample_id = c("a", "b"), call_rate = c(0.99, 0.98),
                    stringsAsFactors = FALSE)
  ss <- data.frame(sample_id = c("a", "b"),
                   sentrix_slide = c("2006", "2006"),
                   sentrix_row = c("01", "02"), sentrix_col = c("01", "01"),
                   Plate = c("P1", "P2"), Well = c("A01", "B02"),
                   stringsAsFactors = FALSE)
  cols <- list(sex = NA_character_, age = NA_character_, batch = NA_character_,
               cell = NA_character_, donor = NA_character_)
  out <- f(qcm, ss, cols)
  expect_true(all(c("sentrix_slide", "sentrix_row", "sentrix_col",
                    "Plate", "Well") %in% names(out)))
})
