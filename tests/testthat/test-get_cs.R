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
  result <- get_cs(topic_res, ldsc_dir, "TestTrait", nTopics,
                   use_approximation = TRUE)

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

  two <- get_cs(topic_res, ldsc_dir, "TestTrait", nTopics,
                use_approximation = TRUE)
  one <- get_cs(topic_res, ldsc_dir, "TestTrait", nTopics,
                use_approximation = TRUE, alternative = "greater")

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
  base  <- get_cs(topic_res, ldsc_dir, "TestTrait", nTopics,
                  use_approximation = TRUE)
  Sig   <- base$Sigma_used * 4
  wide  <- get_cs(topic_res, ldsc_dir, "TestTrait", nTopics, Sigma = Sig)

  expect_equal(wide$z_cell, base$z_cell / 2, tolerance = 1e-10)
  expect_error(get_cs(topic_res, ldsc_dir, "TestTrait", nTopics,
                      Sigma = matrix(1, 2, 2)), "must be a 3 x 3 matrix")

  unlink(ldsc_dir, recursive = TRUE)
})



test_that("cell_type_heterogeneity: weights annihilate the constant vector", {
  set.seed(3); K <- 5; n <- 400
  Pmat <- matrix(0, 2000, K); for (k in 1:K) Pmat[((k-1)*400+1):(k*400), k] <- 1
  Lmat <- gtools::rdirichlet(n, rep(0.3, K))
  tr <- list(Pmat = Pmat, Lmat = Lmat)
  a_k <- colSums(Pmat); W <- sweep(Lmat, 2, a_k, "*"); W <- W / rowSums(W)
  expect_lt(max(abs(cov(W) %*% rep(1, K))), 1e-12)
})

test_that("cell_type_heterogeneity: the unknown common enrichment drops out", {
  set.seed(3); K <- 5; n <- 400
  Pmat <- matrix(0, 2000, K); for (k in 1:K) Pmat[((k-1)*400+1):(k*400), k] <- 1
  Lmat <- gtools::rdirichlet(n, rep(0.3, K))
  tr <- list(Pmat = Pmat, Lmat = Lmat)
  Sig <- diag(c(.5, .4, .6, .3, .45)^2)
  g <- rep("A", n)
  set.seed(9); eps <- as.vector(t(chol(Sig)) %*% rnorm(K))
  # E = c*1 + eps for two very different c must give the same p-value
  p1 <- cell_type_heterogeneity(tr, 1  + eps, Sig, g)$p
  p2 <- cell_type_heterogeneity(tr, 40 + eps, Sig, g)$p
  expect_equal(p1, p2, tolerance = 1e-8)
})

test_that("cell_type_heterogeneity: null p-values are uniform", {
  skip_if_not_installed("CompQuadForm")
  set.seed(3); K <- 5; n <- 400
  Pmat <- matrix(0, 2000, K); for (k in 1:K) Pmat[((k-1)*400+1):(k*400), k] <- 1
  Lmat <- gtools::rdirichlet(n, rep(0.3, K))
  tr <- list(Pmat = Pmat, Lmat = Lmat)
  Sig <- diag(c(.5, .4, .6, .3, .45)^2); R <- chol(Sig)
  g <- rep("A", n)
  set.seed(11)
  pv <- vapply(1:300, function(i)
    cell_type_heterogeneity(tr, 7 + as.vector(t(R) %*% rnorm(K)), Sig, g)$p, 0)
  expect_gt(suppressWarnings(ks.test(pv, "punif")$p.value), 0.01)
  expect_lt(abs(mean(pv < 0.05) - 0.05), 0.04)
})

test_that("cell_type_heterogeneity: input validation and small groups", {
  set.seed(3); K <- 3
  Pmat <- matrix(0, 300, K); for (k in 1:K) Pmat[((k-1)*100+1):(k*100), k] <- 1
  tr <- list(Pmat = Pmat, Lmat = gtools::rdirichlet(50, rep(0.3, K)))
  Sig <- diag(K)
  expect_error(cell_type_heterogeneity(tr, c(1, 2), Sig, rep("A", 50)),
               "enrichment must have length 3")
  expect_error(cell_type_heterogeneity(tr, rep(1, K), diag(2), rep("A", 50)),
               "Sigma must be a 3 x 3 matrix")
  expect_error(cell_type_heterogeneity(tr, rep(1, K), Sig, rep("A", 10)),
               "one entry per cell")
  expect_error(cell_type_heterogeneity(tr, rep(1, K), Sig, rep("A", 50),
                                       min_cells = 100), "no group has at least")
})

test_that("ldsc_jackknife_cov signature no longer requires an opt-in flag", {
  expect_false("allow_unvalidated" %in% names(formals(ldsc_jackknife_cov)))
  expect_true(all(c("baseline_prefix","frq_prefix") %in% names(formals(ldsc_jackknife_cov))))
})

test_that("get_cs defaults to the jackknife and refuses to guess", {
  expect_true("use_approximation" %in% names(formals(get_cs)))
  expect_false(eval(formals(get_cs)$use_approximation))   # jackknife is default
  # with neither Sigma, prefixes, nor the opt-out, it must stop rather than
  # silently fall back to the approximation
  expect_error(
    get_cs(list(Pmat = matrix(1), Lmat = matrix(1)), tempdir(), "TRAIT", 3),
    "needs the covariance")
})

test_that("get_cs(use_approximation = TRUE) still runs", {
  skip_if_not_installed("ashr")
  nTopics <- 3; nPeaks <- 100; nCells <- 40
  pn <- paste0("chr1:", seq(1000, by = 501, length.out = nPeaks),
               "-", seq(1500, by = 501, length.out = nPeaks))
  Pmat <- matrix(0, nPeaks, nTopics)
  Pmat[1:40, 1] <- 1; Pmat[41:70, 2] <- 1; Pmat[71:100, 3] <- 1
  rownames(Pmat) <- pn; colnames(Pmat) <- paste0("k", 1:nTopics)
  set.seed(1); Lmat <- matrix(runif(nCells * nTopics), nCells, nTopics)
  Lmat <- Lmat / rowSums(Lmat); colnames(Lmat) <- paste0("k", 1:nTopics)

  d <- file.path(tempdir(), "mock_approx"); unlink(d, recursive = TRUE)
  for (k in 1:nTopics) {
    rd <- file.path(d, paste0("k", k, "_output"), "results")
    dir.create(rd, recursive = TRUE)
    write.table(data.frame(Category = "L2_0", `Prop._SNPs` = 0.01,
      `Prop._h2` = 0.05 * k, `Prop._h2_std_error` = 0.01,
      Enrichment = 1 + k * 0.5, Enrichment_std_error = 0.2,
      Coefficient = 1e-8, `Coefficient_std_error` = 1e-9,
      `Coefficient_z-score` = 2.0, check.names = FALSE),
      file = file.path(rd, "TestTrait.results"), sep = "\t",
      row.names = FALSE, quote = FALSE)
  }
  r <- get_cs(list(Pmat = Pmat, Lmat = Lmat), d, "TestTrait", nTopics,
              use_approximation = TRUE)
  expect_equal(length(r$cs), nCells)
  expect_equal(length(r$p_cell), nCells)
  unlink(d, recursive = TRUE)
})
