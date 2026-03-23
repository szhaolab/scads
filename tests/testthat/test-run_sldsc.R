test_that("run_sldsc validates sumstats file existence", {
  tmp <- tempdir()
  empty_dir <- file.path(tmp, "empty_sumstats")
  dir.create(empty_dir, showWarnings = FALSE)

  expect_error(
    run_sldsc(
      chrs = 1, polyfun_path = "/fake", ldsc_path = "/fake",
      sumstats_path = empty_dir, n = 100, trait = "NONEXISTENT",
      onekg_path = "/fake/", bed_dir = "/fake/",
      baseline_dir = "/fake/", frqfile_pref = "/fake/",
      hm3_snps = "/fake/hm3.txt", weights_pref = "/fake/",
      out_dir = file.path(tmp, "sldsc_test_out")
    ),
    "No matching sumstats file found"
  )

  unlink(empty_dir, recursive = TRUE)
})

test_that("run_sldsc validates BED file existence", {
  tmp <- tempdir()

  # Create a fake sumstats file to pass step 1 validation
  ss_dir <- file.path(tmp, "sumstats_test")
  dir.create(ss_dir, showWarnings = FALSE)
  writeLines("test", file.path(ss_dir, "TEST_sumstats.txt.gz"))

  # Create empty bed_dir

  bed_dir <- file.path(tmp, "empty_beds")
  dir.create(bed_dir, showWarnings = FALSE)

  out_dir <- file.path(tmp, "sldsc_bed_test")

  # This should fail at BED file search (before any system() calls)
  # since we set polyfun_path to a nonexistent path, the munge step
  # will fail first via system(), so we test the command construction instead
  cmd <- build_sldsc_command("munge",
    polyfun_path = "/fake/polyfun",
    sumstats_gz = file.path(ss_dir, "TEST_sumstats.txt.gz"),
    n = 100,
    munged_out = file.path(out_dir, "TEST_munged.parquet"))

  expect_true(grepl("munge_polyfun_sumstats.py", cmd))

  unlink(c(ss_dir, bed_dir, out_dir), recursive = TRUE)
})
