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
run_fastTopics <- function(count_matrix, nTopics = 10, n_s = 1000, n_c = 1, ...) {
  
  # Ensure that count_matrix has cells as rows and peaks as columns
  # Transpose if necessary 
  # Assuming num_peaks > num_cells
  if (nrow(count_matrix) > ncol(count_matrix)) {
    count_matrix <- Matrix::t(count_matrix)
  }
  print(dim(count_matrix))
  
  # Run fastTopics model
  cat("\nRunning fastTopics\n")
  cat("\nStart time: ")
  print(Sys.time())
  fastTopics_fit <- fastTopics::fit_topic_model(count_matrix, 
                                                k = nTopics, ...)
  
  cat("\nStop time: ")
  print(Sys.time())
  
  # Extract factor matrices
  Fmat <- fastTopics_fit$F  # Peaks x topics matrix (feature loadings)
  Lmat <- fastTopics_fit$L  # Cells x topics matrix (topic proportions)
  
  # Determine per-topic cutoffs using KDE midpoint
  cat("\nDetermining per-topic cutoffs using KDE\n")
  cutoffs <- numeric(nTopics)  # Initialize vector to store cutoffs
  
  for (topic in 1:nTopics) {
    cat(sprintf("Processing Topic %d...\n", topic))
    F_topic <- Fmat[, topic]
    
    # Find cutoff using KDE
    cutoff_topic <- 10^find_kde_midpoint(log10(F_topic))
    cutoffs[topic] <- cutoff_topic
    cat(sprintf("Cutoff for Topic %d: %.5e\n", topic, cutoff_topic))
  }
  
  # # Determine global cutoff using KDE midpoint 
  # cat("\nDetermining global cutoff using KDE\n")
  # global_cutoff <- 10^find_kde_midpoint(log10(c(Fmat)))
  # cat("\nGlobal cutoff: ", global_cutoff)
  
  # # Run differential expression analysis
  # s <- Matrix::rowSums(count_matrix)  # Vector of total counts per cell
  # 
  # cat("\nRunning fastTopics-DE\n")
  # cat("\nStart time: ")
  # print(Sys.time())
  # de_res <- de_analysis2(
  #   fit = fastTopics_fit,
  #   X = count_matrix,
  #   s = s,
  #   lfc.stat = "vsnull",
  #   shrink.method = "none",
  #   control = list(ns = n_s, nc = n_c, minval = 1e-50),
  #   f0 = global_cutoff  # Global KDE cutoff
  # )
  # 
  # print(summary(de_res$z))
  # 
  # # 'de_res$z' is a matrix of z-scores with dimensions peaks x topics
  # # Remove peaks with any NA values for z-scores in topics
  # keep_peaks <- which(complete.cases(de_res$z))
  # keep_z <- de_res$z[keep_peaks, ]  # Select rows corresponding to keep_peaks
  # cat("\nNumber of peaks with complete z-scores: ", length(keep_peaks), "\n")
  # 
  # # Convert z-score matrix to a vector for locfdr input
  # z_scores <- as.vector(keep_z)
  # print(summary(z_scores))
  # 
  # # Initialize locfdr_vals_full with NAs
  # locfdr_vals_full <- rep(NA, length(z_scores))
  # 
  # # Try applying locfdr to all z-scores - Instead of locfdr, try GMM (see scads_p1c2)
  # locfdr_failed <- FALSE  # Flag to track if locfdr failed
  # cat("\nApplying locfdr to all z-scores\n")
  # cat("\nStart time: ")
  # print(Sys.time())
  # 
  # locfdr_res <- tryCatch(
  #   locfdr::locfdr(z_scores, nulltype = 2, plot = 0),
  #   error = function(e) {
  #     cat("\nlocfdr failed with error:\n")
  #     print(e)
  #     locfdr_failed <<- TRUE  # Set the flag to TRUE
  #     return(NULL)
  #   }
  # )
  # 
  # if (!locfdr_failed) {
  #   locfdr_vals <- locfdr_res$fdr  # Extract lfdr values
  #   locfdr_vals_full <- locfdr_vals
  #   
  #   # Reshape locfdr_vals_full back to a matrix with the same dimensions as keep_z
  #   locfdr_mat <- matrix(
  #     locfdr_vals_full,
  #     nrow = nrow(keep_z),
  #     ncol = ncol(keep_z),
  #     byrow = FALSE
  #   )
  #   
  #   # Compute p-values as 1 - lfdr for each peak-topic pair
  #   p_jk <- 1 - locfdr_mat  # Matrix of p-values (peaks x topics)
  #   
  # } else {
  #   cat("\nlocfdr failed. Creating p_jk using per-topic binarization based on KDE cutoffs.\n")
  #   # Adjust Fmat to match keep_peaks
  #   Fmat_kept <- Fmat[keep_peaks, ]
  #   # Initialize p_jk matrix
  #   p_jk <- matrix(0, nrow = nrow(Fmat_kept), ncol = ncol(Fmat_kept))
  #   colnames(p_jk) <- colnames(Fmat_kept)
  #   rownames(p_jk) <- rownames(Fmat_kept)
  #   
  #   # Binarize Fmat per topic using the corresponding cutoff
  #   for (topic in 1:nTopics) {
  #     threshold <- cutoffs[topic]
  #     p_jk[, topic] <- ifelse(Fmat_kept[, topic] > threshold, 1, 0)
  #     cat(sprintf("Binarized Topic %d with threshold %.5e\n", topic, threshold))
  #   }
  # }
  # 
  # # Adjust Fmat to match p_jk
  # Fmat <- Fmat[keep_peaks, ]
  
  # # Binarize Fmat with global_cutoff 
  # p_jk <- matrix(0, nrow = nrow(Fmat), ncol = ncol(Fmat))
  # p_jk <- ifelse(Fmat > global_cutoff, 1, 0)
  
    # Binarize Fmat per topic using the corresponding cutoff
    p_jk <- matrix(0, nrow = nrow(Fmat), ncol = ncol(Fmat))
    for (topic in 1:nTopics) {
      threshold <- cutoffs[topic]
      p_jk[, topic] <- ifelse(Fmat[, topic] > threshold, 1, 0)
      cat(sprintf("Binarized Topic %d with threshold %.5e\n", topic, threshold))
    }
  
  # Return results as a list
  return(list(
    fastTopics_fit = fastTopics_fit,
    Fmat = Fmat,
    Lmat = Lmat,
    # de_res = de_res,
    Pmat = p_jk,  # Matrix of p-values (peaks x topics)
    cutoffs = cutoffs
  ))
}
