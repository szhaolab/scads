test_that("get_cs validates inputs", {
  expect_error(get_cs(list(), "/fake/dir", "TRAIT", 5),
               "topic_res must contain 'Pmat'")

  expect_error(get_cs(list(Pmat = matrix(1)), "/fake/dir", "TRAIT", 5),
               "topic_res must contain 'Lmat'")

  expect_error(
    get_cs(list(Pmat = matrix(1), Lmat = matrix(1)), "/nonexistent/dir", "TRAIT", 5),
    "ldsc_res_dir does not exist"
  )

  tmp <- tempdir()
  expect_error(
    get_cs(list(Pmat = matrix(1), Lmat = matrix(1)), tmp, "", 5),
    "trait must be a non-empty string"
  )
})

test_that("get_cs computes cell scores from mock LDSC results", {
  skip_if_not_installed("ashr")
  skip_if_not_installed("GenomicRanges")

  # Create mock data
  nTopics <- 3
  nPeaks <- 100
  nCells <- 50

  # Peak names in chr_start_end format
  peak_names <- paste0("chr1:", seq(1000, by = 501, length.out = nPeaks),
                       "-", seq(1500, by = 501, length.out = nPeaks))

  # Create Pmat (binary peak-topic assignment)
  Pmat <- matrix(0, nrow = nPeaks, ncol = nTopics)
  Pmat[1:40, 1] <- 1
  Pmat[41:70, 2] <- 1
  Pmat[71:100, 3] <- 1
  rownames(Pmat) <- peak_names
  colnames(Pmat) <- paste0("k", 1:nTopics)

  # Create Lmat (cell-topic loading)
  Lmat <- matrix(runif(nCells * nTopics), nrow = nCells, ncol = nTopics)
  Lmat <- Lmat / rowSums(Lmat)
  colnames(Lmat) <- paste0("k", 1:nTopics)

  topic_res <- list(Pmat = Pmat, Lmat = Lmat)

  # Create mock LDSC results directory
  ldsc_dir <- file.path(tempdir(), "mock_ldsc")
  if (dir.exists(ldsc_dir)) unlink(ldsc_dir, recursive = TRUE)

  for (k in 1:nTopics) {
    res_dir <- file.path(ldsc_dir, paste0("k", k, "_output"), "results")
    dir.create(res_dir, recursive = TRUE)

    # Mock .results file
    mock_res <- data.frame(
      Category = paste0("L2_0"),
      `Prop._SNPs` = 0.01,
      `Prop._h2` = 0.05 * k,
      `Prop._h2_std_error` = 0.01,
      Enrichment = 1 + k * 0.5,
      Enrichment_std_error = 0.2,
      Coefficient = 1e-8,
      `Coefficient_std_error` = 1e-9,
      `Coefficient_z-score` = 2.0,
      check.names = FALSE
    )
    write.table(mock_res,
                file = file.path(res_dir, "TestTrait.results"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }

  # Run get_cs
  result <- get_cs(topic_res, ldsc_dir, "TestTrait", nTopics)

  expect_true(is.list(result))
  expect_named(result, c("cs", "z_cell", "p_cell", "ldsc_res_table", "ash_res", "cs_dat"))
  expect_equal(length(result$cs), nCells)
  expect_true(all(is.finite(result$cs)))
  expect_true(all(result$cs > 0))

  # cs_dat should contain M_i and N_i
  expect_true(!is.null(result$cs_dat$M_i))
  expect_true(!is.null(result$cs_dat$N_i))
  expect_equal(length(result$cs_dat$M_i), nCells)

  # Clean up
  unlink(ldsc_dir, recursive = TRUE)
})
