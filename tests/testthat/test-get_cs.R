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
  expect_named(result, c("cs", "z_cell", "p_cell", "var_cell", "Sigma_used",
                         "ldsc_res_table", "ash_res", "cs_dat"))
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

test_that("variance uses the full symmetric quadratic form", {
  # Two perfectly correlated topics with equal weight: the score equals a single
  # enrichment estimate, so Var(cs) = sigma^2 exactly. The one-sided (upper
  # triangle) sum returns 0.75 * sigma^2.
  u  <- c(0.5, 0.5)
  sg <- 1
  C  <- matrix(1, 2, 2)
  v  <- u * sg

  full  <- as.numeric(t(v) %*% C %*% v)
  upper <- as.numeric(t(v) %*% (C * upper.tri(C, diag = TRUE)) %*% v)

  expect_equal(full, sg^2)
  expect_equal(upper, 0.75 * sg^2)

  # Independent topics: the two forms agree.
  C0 <- diag(2)
  expect_equal(as.numeric(t(v) %*% C0 %*% v),
               as.numeric(t(v) %*% (C0 * upper.tri(C0, diag = TRUE)) %*% v))
})

test_that("p_cell keeps length I and honours `alternative`", {
  skip_if_not_installed("ashr")
  skip_if_not_installed("GenomicRanges")

  nTopics <- 3; nPeaks <- 100; nCells <- 50
  peak_names <- paste0("chr1:", seq(1000, by = 501, length.out = nPeaks),
                       "-", seq(1500, by = 501, length.out = nPeaks))
  Pmat <- matrix(0, nrow = nPeaks, ncol = nTopics)
  Pmat[1:40, 1] <- 1; Pmat[41:70, 2] <- 1; Pmat[71:100, 3] <- 1
  rownames(Pmat) <- peak_names; colnames(Pmat) <- paste0("k", 1:nTopics)

  set.seed(1)
  Lmat <- matrix(runif(nCells * nTopics), nrow = nCells, ncol = nTopics)
  Lmat <- Lmat / rowSums(Lmat); colnames(Lmat) <- paste0("k", 1:nTopics)
  topic_res <- list(Pmat = Pmat, Lmat = Lmat)

  ldsc_dir <- file.path(tempdir(), "mock_ldsc_sided")
  if (dir.exists(ldsc_dir)) unlink(ldsc_dir, recursive = TRUE)
  for (k in 1:nTopics) {
    res_dir <- file.path(ldsc_dir, paste0("k", k, "_output"), "results")
    dir.create(res_dir, recursive = TRUE)
    write.table(data.frame(
      Category = "L2_0", `Prop._SNPs` = 0.01, `Prop._h2` = 0.05 * k,
      `Prop._h2_std_error` = 0.01, Enrichment = 1 + k * 0.5,
      Enrichment_std_error = 0.2, Coefficient = 1e-8,
      `Coefficient_std_error` = 1e-9, `Coefficient_z-score` = 2.0,
      check.names = FALSE),
      file = file.path(res_dir, "TestTrait.results"),
      sep = "\t", row.names = FALSE, quote = FALSE)
  }

  two <- get_cs(topic_res, ldsc_dir, "TestTrait", nTopics)
  one <- get_cs(topic_res, ldsc_dir, "TestTrait", nTopics,
                alternative = "greater")

  # p_cell must align with cs / z_cell
  expect_equal(length(two$p_cell), nCells)
  expect_equal(length(two$p_cell), length(two$cs))
  expect_equal(length(two$p_cell), length(two$z_cell))

  # all topics enriched here, so z > 0 and one-sided p is half the two-sided p
  expect_true(all(two$z_cell > 0))
  expect_equal(one$p_cell, two$p_cell / 2, tolerance = 1e-10)

  # variance is returned and positive
  expect_equal(length(two$var_cell), nCells)
  expect_true(all(two$var_cell > 0))

  unlink(ldsc_dir, recursive = TRUE)
})

test_that("supplied Sigma overrides the annotation-correlation fallback", {
  skip_if_not_installed("ashr")
  skip_if_not_installed("GenomicRanges")

  nTopics <- 3; nPeaks <- 100; nCells <- 20
  peak_names <- paste0("chr1:", seq(1000, by = 501, length.out = nPeaks),
                       "-", seq(1500, by = 501, length.out = nPeaks))
  Pmat <- matrix(0, nrow = nPeaks, ncol = nTopics)
  Pmat[1:40, 1] <- 1; Pmat[41:70, 2] <- 1; Pmat[71:100, 3] <- 1
  rownames(Pmat) <- peak_names; colnames(Pmat) <- paste0("k", 1:nTopics)

  set.seed(2)
  Lmat <- matrix(runif(nCells * nTopics), nrow = nCells, ncol = nTopics)
  Lmat <- Lmat / rowSums(Lmat); colnames(Lmat) <- paste0("k", 1:nTopics)
  topic_res <- list(Pmat = Pmat, Lmat = Lmat)

  ldsc_dir <- file.path(tempdir(), "mock_ldsc_sigma")
  if (dir.exists(ldsc_dir)) unlink(ldsc_dir, recursive = TRUE)
  for (k in 1:nTopics) {
    res_dir <- file.path(ldsc_dir, paste0("k", k, "_output"), "results")
    dir.create(res_dir, recursive = TRUE)
    write.table(data.frame(
      Category = "L2_0", `Prop._SNPs` = 0.01, `Prop._h2` = 0.05 * k,
      `Prop._h2_std_error` = 0.01, Enrichment = 1 + k * 0.5,
      Enrichment_std_error = 0.2, Coefficient = 1e-8,
      `Coefficient_std_error` = 1e-9, `Coefficient_z-score` = 2.0,
      check.names = FALSE),
      file = file.path(res_dir, "TestTrait.results"),
      sep = "\t", row.names = FALSE, quote = FALSE)
  }

  # inflate the covariance 4x -> SE doubles -> |z| halves
  base  <- get_cs(topic_res, ldsc_dir, "TestTrait", nTopics)
  Sig   <- base$Sigma_used * 4
  wide  <- get_cs(topic_res, ldsc_dir, "TestTrait", nTopics, Sigma = Sig)

  expect_equal(wide$z_cell, base$z_cell / 2, tolerance = 1e-10)
  expect_error(get_cs(topic_res, ldsc_dir, "TestTrait", nTopics,
                      Sigma = matrix(1, 2, 2)), "must be a 3 x 3 matrix")

  unlink(ldsc_dir, recursive = TRUE)
})


test_that("ldsc_jackknife_cov refuses to run unless explicitly opted in", {
  expect_error(ldsc_jackknife_cov("/nonexistent", "TRAIT", 3),
               "EXPERIMENTAL")
})
