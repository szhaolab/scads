# Mock S-LDSC directory: per-topic .results, plus (optionally) baseline and
# thin topic annotations for chromosomes 1-2.
mock_sldsc <- function(d, enrich, prop = rep(0.01, length(enrich)),
                       se = rep(0.2, length(enrich)), annot = NULL, cm = NULL) {
  unlink(d, recursive = TRUE)
  K <- length(enrich)
  for (k in seq_len(K)) {
    rd <- file.path(d, paste0("k", k, "_output"), "results")
    dir.create(rd, recursive = TRUE)
    write.table(data.frame(Category = "L2_0", `Prop._SNPs` = prop[k],
      `Prop._h2` = 0.05, `Prop._h2_std_error` = 0.01, Enrichment = enrich[k],
      Enrichment_std_error = se[k], Coefficient = 1e-8,
      `Coefficient_std_error` = 1e-9, `Coefficient_z-score` = 2.0,
      check.names = FALSE), file = file.path(rd, "TestTrait.results"),
      sep = "\t", row.names = FALSE, quote = FALSE)
  }
  if (!is.null(annot)) {
    for (chr in names(cm)) {
      n <- length(cm[[chr]])
      data.table::fwrite(data.frame(CHR = as.integer(chr), BP = seq_len(n),
                                    SNP = paste0("rs", chr, "_", seq_len(n)),
                                    CM = cm[[chr]], base = 1),
                         file.path(d, paste0("baseline.", chr, ".annot.gz")), sep = "\t")
      for (k in seq_len(K)) {
        ad <- file.path(d, paste0("k", k, "_output"), "annotations", "TestTrait")
        dir.create(ad, recursive = TRUE, showWarnings = FALSE)
        data.table::fwrite(data.frame(ANNOT = annot[[k]][[chr]]),
                           file.path(ad, paste0("TestTrait.", chr, ".annot.gz")), sep = "\t")
      }
    }
  }
  d
}

mock_topics <- function(nTopics = 3, nCells = 40, seed = 1) {
  nPeaks <- 99
  pn <- paste0("chr1:", seq(1000, by = 501, length.out = nPeaks),
               "-", seq(1500, by = 501, length.out = nPeaks))
  Pmat <- matrix(0, nPeaks, nTopics)
  for (k in seq_len(nTopics)) Pmat[(k - 1) * (nPeaks / nTopics) + 1:(nPeaks / nTopics), k] <- 1
  rownames(Pmat) <- pn; colnames(Pmat) <- paste0("k", seq_len(nTopics))
  set.seed(seed); Lmat <- matrix(runif(nCells * nTopics), nCells, nTopics)
  Lmat <- Lmat / rowSums(Lmat); colnames(Lmat) <- paste0("k", seq_len(nTopics))
  list(Pmat = Pmat, Lmat = Lmat)
}

test_that("annot_block_coverage counts 0.01 cM blocks per topic", {
  skip_if_not_installed("data.table")
  # chr1: 10 SNPs in 5 blocks (2 per block); chr2: 6 SNPs in 3 blocks -> 8 blocks
  cm <- list(`1` = rep(c(0.001, 0.011, 0.021, 0.031, 0.041), each = 2),
             `2` = rep(c(1.005, 1.015, 1.025), each = 2))
  annot <- list(
    list(`1` = c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0), `2` = c(0, 0, 0, 0, 0, 0)),  # 2 SNPs, 1 block
    list(`1` = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0), `2` = c(1, 0, 1, 0, 1, 0)))  # 8 SNPs, 8 blocks
  d <- mock_sldsc(file.path(tempdir(), "mock_cov"), enrich = c(2, 3),
                  annot = annot, cm = cm)
  cv <- annot_block_coverage(d, "TestTrait", 2, file.path(d, "baseline."), chrs = 1:2)
  expect_equal(attr(cv, "n_blocks_genome"), 8)
  expect_equal(cv$n_snps, c(2, 8))
  expect_equal(cv$n_blocks, c(1, 8))
  expect_equal(cv$pct_blocks, c(12.5, 100))
  expect_equal(cv$prop_snps, c(2, 8) / 16)

  # rows must match the baseline, and a genetic map is required
  annot[[1]][["2"]] <- c(0, 0)
  d <- mock_sldsc(d, enrich = c(2, 3), annot = annot, cm = cm)
  expect_error(annot_block_coverage(d, "TestTrait", 2, file.path(d, "baseline."), chrs = 1:2),
               "thin annotations matching")
  annot[[1]][["2"]] <- rep(0, 6); cm[["2"]] <- rep(0, 6)
  d <- mock_sldsc(d, enrich = c(2, 3), annot = annot, cm = cm)
  expect_error(annot_block_coverage(d, "TestTrait", 2, file.path(d, "baseline."), chrs = 1:2),
               "all zero")
  unlink(d, recursive = TRUE)
})

test_that("get_cs floors topics below min_pct_blocks, regardless of Prop._SNPs", {
  skip_if_not_installed("ashr"); skip_if_not_installed("GenomicRanges")
  tr <- mock_topics()
  # k1: 1% of SNPs but clustered (1.0% of blocks) -> small
  # k3: 0.3% of SNPs but spread (5% of blocks) -> kept (the old rule floored it)
  d <- mock_sldsc(file.path(tempdir(), "mock_floor"), enrich = c(3, 2, 4),
                  prop = c(0.01, 0.02, 0.003))
  expect_warning(r <- get_cs(tr, d, "TestTrait", 3, use_approximation = TRUE,
                             block_coverage = c(1.0, 10, 5)),
                 "1 topic\\(s\\) below .*k1 \\(1%\\).*floored")
  expect_equal(r$topic_qc$small, c(TRUE, FALSE, FALSE))
  expect_equal(r$topic_qc$floored, c(TRUE, FALSE, FALSE))
  expect_equal(r$topic_qc$pct_blocks, c(1.0, 10, 5))
  # a floored topic adds nothing to the score numerator: cs depends only on k2, k3
  expect_equal(r$cs_dat$M_i - r$cs_dat$N_i,
               as.numeric(tr$Lmat %*% (colSums(tr$Pmat) * (c(1, r$ash_res$PosteriorMean[2:3]) - 1))))

  # small_annot = "warn" keeps the topic and still warns
  expect_warning(w <- get_cs(tr, d, "TestTrait", 3, use_approximation = TRUE,
                             block_coverage = c(1.0, 10, 5), small_annot = "warn"),
                 "kept because small_annot")
  expect_equal(w$topic_qc$floored, c(FALSE, FALSE, FALSE))
  expect_false(isTRUE(all.equal(w$cs, r$cs)))

  # a data frame from annot_block_coverage() is accepted, and its length checked
  expect_warning(get_cs(tr, d, "TestTrait", 3, use_approximation = TRUE,
                        block_coverage = data.frame(pct_blocks = c(1.0, 10, 5))), "k1")
  expect_error(suppressWarnings(get_cs(tr, d, "TestTrait", 3, use_approximation = TRUE,
                                       block_coverage = c(1, 2))), "one value per topic")
  unlink(d, recursive = TRUE)
})

test_that("get_cs falls back to Prop._SNPs, and min_prop_snps is deprecated", {
  skip_if_not_installed("ashr"); skip_if_not_installed("GenomicRanges")
  tr <- mock_topics()
  d <- mock_sldsc(file.path(tempdir(), "mock_fallback"), enrich = c(3, 2, 4),
                  prop = c(0.01, 0.02, 0.003))
  # no baseline_prefix and no block_coverage: old 0.5% rule, with a message
  expect_message(expect_warning(r <- get_cs(tr, d, "TestTrait", 3, use_approximation = TRUE),
                                "k3"), "Block coverage unavailable")
  expect_equal(r$topic_qc$floored, c(FALSE, FALSE, TRUE))
  expect_true(all(is.na(r$topic_qc$pct_blocks)))

  # coverage that cannot be computed (e.g. missing annotation files) also falls back
  expect_message(expect_warning(get_cs(tr, d, "TestTrait", 3, Sigma = diag(0.04, 3),
                                       baseline_prefix = file.path(d, "nope.")), "k3"),
                 "Could not compute 0.01 cM block coverage")

  expect_warning(expect_warning(o <- get_cs(tr, d, "TestTrait", 3, use_approximation = TRUE,
                                            min_prop_snps = 0.015), "deprecated"), "k1, k3")
  expect_equal(o$topic_qc$floored, c(TRUE, FALSE, TRUE))
  unlink(d, recursive = TRUE)
})

test_that("floored topics contribute no variance when Sigma is supplied", {
  skip_if_not_installed("ashr"); skip_if_not_installed("GenomicRanges")
  tr <- mock_topics()
  d <- mock_sldsc(file.path(tempdir(), "mock_sigma"), enrich = c(3, 2, 4))
  Sig <- matrix(0.02, 3, 3); diag(Sig) <- c(0.09, 0.04, 0.16)
  expect_warning(r <- get_cs(tr, d, "TestTrait", 3, Sigma = Sig, block_coverage = c(1, 10, 5)), "k1")
  W <- sweep(tr$Lmat, 2, colSums(tr$Pmat), "*") / sum(colSums(tr$Pmat))
  S0 <- Sig; S0[1, ] <- 0; S0[, 1] <- 0
  expect_equal(r$var_cell, rowSums((W %*% S0) * W))
  # with every topic passing, the supplied Sigma is used unchanged
  ok <- get_cs(tr, d, "TestTrait", 3, Sigma = Sig, block_coverage = c(5, 10, 5))
  expect_equal(ok$var_cell, rowSums((W %*% Sig) * W))
  unlink(d, recursive = TRUE)
})

test_that("get_cs computes block coverage from baseline_prefix", {
  skip_if_not_installed("ashr"); skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("data.table")
  tr <- mock_topics()
  # 22 chromosomes x 10 blocks of one SNP each = 220 blocks
  chrs <- as.character(1:22)
  cm <- setNames(lapply(chrs, function(c) seq(0.001, by = 0.01, length.out = 10)), chrs)
  one <- function(v) setNames(lapply(chrs, function(c) v), chrs)
  annot <- list(one(c(1, rep(0, 9))),                                     # 22/220 = 10%
                c(list(`1` = c(1, rep(0, 9))), one(rep(0, 10))[-1]),       # 1/220  -> small
                one(rep(c(1, 0), 5)))                                     # 50%
  d <- mock_sldsc(file.path(tempdir(), "mock_prefix"), enrich = c(3, 2, 4),
                  annot = annot, cm = cm)
  Sig <- diag(0.04, 3)
  expect_warning(r <- get_cs(tr, d, "TestTrait", 3, Sigma = Sig,
                             baseline_prefix = file.path(d, "baseline.")), "k2 \\(0.45%\\)")
  expect_equal(r$topic_qc$pct_blocks, c(10, 100 / 220, 50))
  expect_equal(r$topic_qc$floored, c(FALSE, TRUE, FALSE))
  unlink(d, recursive = TRUE)
})
