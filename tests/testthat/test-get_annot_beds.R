test_that("get_annot_beds creates BED files for each topic", {
  # Create mock topic results
  n_peaks <- 50
  n_topics <- 3

  Fmat <- matrix(runif(n_peaks * n_topics), nrow = n_peaks, ncol = n_topics)
  colnames(Fmat) <- paste0("k", 1:n_topics)
  rownames(Fmat) <- paste0("chr1_", seq(1000, by = 500, length.out = n_peaks),
                           "_", seq(1500, by = 500, length.out = n_peaks))

  # Binary assignment matrix: each peak assigned to one topic
  Pmat <- matrix(0, nrow = n_peaks, ncol = n_topics)
  Pmat[1:20, 1] <- 1
  Pmat[21:35, 2] <- 1
  Pmat[36:50, 3] <- 1
  colnames(Pmat) <- colnames(Fmat)
  rownames(Pmat) <- rownames(Fmat)

  topics_res <- list(Fmat = Fmat, Pmat = Pmat)

  out_dir <- tempdir()
  test_dir <- file.path(out_dir, "test_beds")
  if (dir.exists(test_dir)) unlink(test_dir, recursive = TRUE)

  bed_dirs <- get_annot_beds(topics_res, output_dir = test_dir)

  # Check that output is a list with one entry per topic

  expect_equal(length(bed_dirs), n_topics)
  expect_named(bed_dirs, paste0("k", 1:n_topics))

  # Check that BED files were created
  for (k in paste0("k", 1:n_topics)) {
    bed_file <- file.path(bed_dirs[[k]], paste0(k, "_annotations.bed"))
    expect_true(file.exists(bed_file), info = paste("BED file missing for", k))
  }

  # Check BED file content for k1
  bed_k1 <- read.table(file.path(bed_dirs[["k1"]], "k1_annotations.bed"),
                        header = FALSE, sep = "\t")
  expect_equal(nrow(bed_k1), 20)
  expect_equal(ncol(bed_k1), 3)
  expect_true(all(grepl("^chr", bed_k1$V1)))

  # Clean up
  unlink(test_dir, recursive = TRUE)
})

test_that("get_annot_beds handles colon-dash peak format", {
  Fmat <- matrix(1, nrow = 5, ncol = 2)
  Pmat <- matrix(1, nrow = 5, ncol = 2)
  rownames(Fmat) <- c("chr1:100-200", "chr2:300-400", "chr3:500-600",
                       "chr4:700-800", "chr5:900-1000")
  rownames(Pmat) <- rownames(Fmat)
  colnames(Fmat) <- colnames(Pmat) <- paste0("k", 1:2)

  topics_res <- list(Fmat = Fmat, Pmat = Pmat)
  out_dir <- file.path(tempdir(), "test_colondash")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)

  bed_dirs <- get_annot_beds(topics_res, output_dir = out_dir)

  bed_k1 <- read.table(file.path(bed_dirs[["k1"]], "k1_annotations.bed"),
                        header = FALSE, sep = "\t")
  expect_equal(bed_k1$V1[1], "chr1")
  expect_equal(bed_k1$V2[1], 100)

  unlink(out_dir, recursive = TRUE)
})

test_that("get_annot_beds errors when rownames are missing", {
  Fmat <- matrix(1, nrow = 3, ncol = 2)
  Pmat <- matrix(1, nrow = 3, ncol = 2)
  # No rownames set
  topics_res <- list(Fmat = Fmat, Pmat = Pmat)

  expect_error(get_annot_beds(topics_res, output_dir = tempdir()),
               "Row names must contain peak ranges")
})

test_that("get_annot_beds warns on empty topics", {
  Fmat <- matrix(runif(10), nrow = 5, ncol = 2)
  Pmat <- matrix(0, nrow = 5, ncol = 2)  # no peaks pass cutoff
  rownames(Fmat) <- paste0("chr1_", 1:5 * 100, "_", 1:5 * 100 + 100)
  rownames(Pmat) <- rownames(Fmat)
  colnames(Fmat) <- colnames(Pmat) <- paste0("k", 1:2)

  topics_res <- list(Fmat = Fmat, Pmat = Pmat)
  out_dir <- file.path(tempdir(), "test_empty_topic")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)

  expect_warning(
    get_annot_beds(topics_res, output_dir = out_dir),
    "No peaks passed the cutoff"
  )

  unlink(out_dir, recursive = TRUE)
})
