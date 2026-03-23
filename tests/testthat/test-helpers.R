test_that("parse_peak_ranges handles underscore format", {
  peaks <- c("chr1_100_200", "chr2_300_400", "chrX_500_600")
  result <- parse_peak_ranges(peaks)

  expect_equal(nrow(result), 3)
  expect_equal(result$seqnames, c("chr1", "chr2", "chrX"))
  expect_equal(result$start, c(100L, 300L, 500L))
  expect_equal(result$end, c(200L, 400L, 600L))
})

test_that("parse_peak_ranges handles colon-dash format", {
  peaks <- c("chr1:100-200", "chr2:300-400")
  result <- parse_peak_ranges(peaks)

  expect_equal(nrow(result), 2)
  expect_equal(result$seqnames, c("chr1", "chr2"))
  expect_equal(result$start, c(100L, 300L))
  expect_equal(result$end, c(200L, 400L))
})

test_that("parse_peak_ranges handles dash format", {
  peaks <- c("chr1-100-200", "chr2-300-400")
  result <- parse_peak_ranges(peaks)

  expect_equal(nrow(result), 2)
  expect_equal(result$seqnames, c("chr1", "chr2"))
  expect_equal(result$start, c(100L, 300L))
  expect_equal(result$end, c(200L, 400L))
})

test_that("parse_peak_ranges adds chr prefix when missing", {
  peaks <- c("1_100_200", "2_300_400")
  result <- parse_peak_ranges(peaks)

  expect_equal(result$seqnames, c("chr1", "chr2"))
})

test_that("parse_peak_ranges errors on invalid format", {
  expect_error(parse_peak_ranges(c("invalid_format")), "must be one of")
  expect_error(parse_peak_ranges(NULL), "non-empty")
  expect_error(parse_peak_ranges(character(0)), "non-empty")
})

test_that("compute_topic_pvalues produces correct output structure", {
  set.seed(42)
  z_mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  colnames(z_mat) <- paste0("k", 1:5)

  result <- compute_topic_pvalues(z_mat)

  expect_true(is.list(result))
  expect_named(result, c("pvals", "qvals", "p_jk", "n_sig_per_topic"))
  expect_equal(dim(result$pvals), c(20, 5))
  expect_equal(dim(result$qvals), c(20, 5))
  expect_equal(dim(result$p_jk), c(20, 5))
  expect_equal(length(result$n_sig_per_topic), 5)
})

test_that("compute_topic_pvalues returns binary p_jk", {
  set.seed(42)
  z_mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  result <- compute_topic_pvalues(z_mat)

  expect_true(all(result$p_jk %in% c(0, 1)))
})

test_that("compute_topic_pvalues: strong signals produce significant peaks", {
  # Create z-scores with very strong signals in column 1
  z_mat <- matrix(0, nrow = 10, ncol = 2)
  z_mat[1:5, 1] <- 20  # extremely strong signal
  z_mat[6:10, 2] <- -2  # weak/negative

  result <- compute_topic_pvalues(z_mat, fdr_cutoff = 0.05)

  # Strong positive z-scores should produce significant peaks
  expect_true(sum(result$p_jk[, 1]) >= 1)
  # Negative z-scores should not
  expect_equal(sum(result$p_jk[, 2]), 0)
})

test_that("compute_topic_pvalues validates inputs", {
  expect_error(compute_topic_pvalues("not a matrix"), "must be a matrix")
  expect_error(compute_topic_pvalues(matrix(1, 2, 2), fdr_cutoff = 0), "between 0 and 1")
  expect_error(compute_topic_pvalues(matrix(1, 2, 2), fdr_cutoff = 1), "between 0 and 1")
})

test_that("compute_topic_pvalues handles NA values", {
  z_mat <- matrix(c(1, NA, 3, 4), nrow = 2, ncol = 2)
  result <- compute_topic_pvalues(z_mat)

  # NAs in z -> NA in pvals -> qval set to 1 -> p_jk = 0
  expect_equal(result$p_jk[1, 2], 0)
})

test_that("build_sldsc_command constructs munge command correctly", {
  cmd <- build_sldsc_command("munge",
    polyfun_path = "/path/to/polyfun",
    sumstats_gz = "/path/to/trait.sumstats.txt.gz",
    n = 60000,
    munged_out = "/output/trait_munged.parquet")

  expect_true(grepl("munge_polyfun_sumstats.py", cmd))
  expect_true(grepl("--sumstats", cmd))
  expect_true(grepl("--n 60000", cmd))
  expect_true(grepl("--out", cmd))
})

test_that("build_sldsc_command constructs annot command correctly", {
  cmd <- build_sldsc_command("annot",
    ldsc_path = "/path/to/ldsc",
    user_bed = "/path/to/peaks.bed",
    bim_file = "/path/to/chr1.bim",
    out_annot = "/output/trait.1.annot.gz")

  expect_true(grepl("make_annot.py", cmd))
  expect_true(grepl("--bed-file", cmd))
  expect_true(grepl("--bimfile", cmd))
})

test_that("build_sldsc_command constructs l2 command correctly", {
  cmd <- build_sldsc_command("l2",
    polyfun_path = "/path/to/polyfun",
    bfile_chr = "/path/to/1000G.1",
    hm3_snps = "/path/to/hm3.txt",
    out_annot = "/output/trait.1.annot.gz",
    out_prefix = "/output/trait.1")

  expect_true(grepl("ldsc.py", cmd))
  expect_true(grepl("--l2", cmd))
  expect_true(grepl("--thin-annot", cmd))
})

test_that("build_sldsc_command constructs h2 command for single chr", {
  cmd <- build_sldsc_command("h2",
    polyfun_path = "/path/polyfun",
    munged_out = "/out/munged.parquet",
    ref_ld = "/out/trait.1",
    baseline_dir = "/base/baselineLD.1",
    frqfile = "/frq/1000G.EUR.QC.1",
    w_ld = "/weights/weights.1",
    final_out = "/out/results/trait",
    use_chr = FALSE)

  expect_true(grepl("--h2", cmd))
  expect_true(grepl("--ref-ld", cmd))
  expect_false(grepl("--ref-ld-chr", cmd))
})

test_that("build_sldsc_command constructs h2 command for all chr", {
  cmd <- build_sldsc_command("h2",
    polyfun_path = "/path/polyfun",
    munged_out = "/out/munged.parquet",
    ref_ld = "/out/trait.",
    baseline_dir = "/base/baselineLD.",
    frqfile = "/frq/1000G.EUR.QC.",
    w_ld = "/weights/weights.",
    final_out = "/out/results/trait",
    use_chr = TRUE)

  expect_true(grepl("--ref-ld-chr", cmd))
  expect_true(grepl("--frqfile-chr", cmd))
  expect_true(grepl("--w-ld-chr", cmd))
})

test_that("build_sldsc_command quotes paths with spaces", {
  cmd <- build_sldsc_command("munge",
    polyfun_path = "/path with spaces/polyfun",
    sumstats_gz = "/my data/trait.txt.gz",
    n = 100,
    munged_out = "/my output/out.parquet")

  # shQuote wraps in single quotes on unix
  expect_true(grepl("'", cmd) || grepl('"', cmd))
})

test_that("build_sldsc_command errors on unknown step", {
  expect_error(build_sldsc_command("invalid"), "Unknown step")
})

test_that("scads_log prints when verbose=TRUE", {
  expect_output(scads_log("test message", verbose = TRUE), "test message")
})

test_that("scads_log suppresses when verbose=FALSE", {
  expect_silent(scads_log("test message", verbose = FALSE))
})
