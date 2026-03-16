#' Run fastTopics and Generate Topic Annotations
#'
#' This function runs the fastTopics model on the provided count matrix,
#' performs differential expression analysis, computes p-values, and returns
#' the factor matrices along with differential expression results.
#'
#' @param count_matrix A numeric matrix or sparse matrix of counts with cells as rows and peaks as columns.
#' @param nTopics An integer specifying the number of topics to fit. Default is 10.
#' @param n_s Number of samples for the DE analysis (default: 1000).
#' @param n_c Number of cores for parallel processing in DE analysis (default: 1).
#' @param ... Additional arguments passed to `fastTopics::fit_topic_model()`.
#' @return A list containing the factor matrices `Fmat` and `Lmat`, differential expression results `de_res`, and p-values matrix `p_jk`.
#' @details
#' The function ensures that the count matrix has cells as rows and peaks as columns.
#' It then runs the fastTopics model, performs differential expression analysis,
#' computes local false discovery rates (lfdr), and calculates p-values for each peak-topic pair.
#' If `locfdr::locfdr` fails, it binarizes `Fmat` based on a threshold and computes `p_jk`.
#' @examples
#' \dontrun{
#' # Assuming 'counts' is your count matrix
#' results <- run_fastTopics(counts, nTopics = 5)
#' }
#' @export
run_fastTopics <- function(count_matrix, nTopics = 10, n_s = 1000, n_c = 1,
                           baseline_method = "constant",
                           bl_celltype_cells = NULL, 
                           bl_celltype_peak_file = NULL,
                           fdr_cutoff = 0.05,
                           outdir, ...) {
  
  
  if (is.null(rownames(count_matrix))) {
    rownames(count_matrix) <- paste0("cell_", 1:nrow(count_matrix))
  }
  
  cat("\nRemoving columns with all zeroes\n")
  count_matrix_filtered <- count_matrix[, colSums(count_matrix) > 0]
  cat("\nNew dimensions\n", dim(count_matrix_filtered))
  rm(count_matrix)
  
  # Run fastTopics model
  cat("\nRunning fastTopics\n")
  cat("\nStart time: ")
  print(Sys.time())

  # if using a previously generate fastTopics object
  # (1) model fit only
  # fastTopics_fit <- readRDS(file.path(outdir, "fasTopics_fit.rds"))
  # (2) model fit and DE results
  # out1 <- readRDS(file.path(outdir, "run_fastTopics_res.rds"))
  # fastTopics_fit <- out1$fastTopics_fit
  # saveRDS(fastTopics_fit, file.path(outdir, "fasTopics_fit.rds"))
  
  # run fastTopics
  fastTopics_fit <- fastTopics::fit_topic_model(count_matrix_filtered,
                                               k = nTopics, ...)
  saveRDS(fastTopics_fit, file.path(outdir, "fasTopics_fit.rds"))

  cat("\nStop time: ")
  print(Sys.time())
  
  # Extract factor matrices
  Fmat <- fastTopics_fit$F  # Peaks x topics matrix (feature loadings)
  Lmat <- fastTopics_fit$L  # Cells x topics matrix (topic proportions)
  
  # Determine baseline via 3 methods 
  cat("\nDetermining baseline\n")
  if (baseline_method == "gc"){
    
    gc_baseline_res <- get_gc_baseline(count_matrix_filtered, Lmat, Fmat, outdir, plot=FALSE) # peak-by-topic; require the column names to be peak IDs with chr_start_end format
    baseline <- gc_baseline_res$lambda_jk # mu0 / N_k
    print(baseline[1:5,1:nTopics])
    cat("Average GC baseline by topics: ", paste(round(colMeans(baseline, na.rm = TRUE), 10)))
    
  } else if (baseline_method == "average"){
    
      total_peaks <- ncol(count_matrix_filtered)
      total_reads <- sum(count_matrix_filtered)
      n_cells <- nrow(count_matrix_filtered)
      reads_per_cell <- total_reads/n_cells
      baseline <- (total_reads)/(total_peaks)/n_cells/reads_per_cell
      cat("Baseline: ", baseline)

  } else if (baseline_method == "estimate") {
    
      baseline <<- get_average_bg(count_matrix_filtered,
                               cell_type = bl_celltype_cells,
                               cell_type_peaks = bl_celltype_peak_file)
      
      cat("Baseline: ", baseline)
    
  } else {
    
      baseline <<- 1e-7
      cat("Baseline: ", baseline)
  }
  
  
  # Run differential expression analysis
  s <- Matrix::rowSums(count_matrix_filtered)  # Vector of total counts per cell

  cat("\nRunning fastTopics-DE\n")
  cat("\nStart time: ")
  print(Sys.time())
  de_res <- de_analysis2(
    fit = fastTopics_fit,
    X = count_matrix_filtered,
    s = s,
    lfc.stat = "vsnull",
    shrink.method = "none",
    control = list(ns = n_s, nc = n_c, minval = 1e-50),
    f0 = baseline
  )
  # de_res <- out1$de_res
  print(summary(de_res$z))

  # Binarize Fmat per topic using z-score from DE and apply FDR
  total_bins <- ceiling(3e9 / 500) # number of genomic bins that could be tested
  # nrow(Fmat) * nTopics # actual number of tests
  ## Get one-side p-value
  pvals <- 1 - pnorm(de_res$z)
  # pvals <- 2 * pnorm(-abs(de_res$z)) # two-sided - previously used
  ## Get adjusted p‑values
  qvals <- apply(pvals, 2, function(col) p.adjust(col,
                                                  method = "fdr",
                                                  n = total_bins))
  qvals[is.na(qvals)] <- 1 # set missing to largest p-value
  p_jk <- ifelse(qvals < fdr_cutoff, 1, 0)
  p_jk[p_jk<0] <- 0
  print(colSums(p_jk))

  # Clean up 
  rm(count_matrix_filtered)
  
  # Return results as a list
  return(list(
    fastTopics_fit = fastTopics_fit,
    Fmat = Fmat,
    Lmat = Lmat,
    de_res = de_res,
    Pmat = p_jk,  # Matrix of p-values (peaks x topics)
    baseline = baseline
  ))
}
