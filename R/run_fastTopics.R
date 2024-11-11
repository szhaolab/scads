#' Run fastTopics and Generate Topic Annotations
#'
#' This function runs the fastTopics model on the provided count matrix,
#' performs differential expression analysis, computes p-values, and returns
#' the factor matrices along with differential expression results.
#'
#'
#' @param count_matrix A numeric matrix or sparse matrix of counts with cells as rows and peaks as columns.
#' @param nTopics An integer specifying the number of topics to fit. Default is 10.
#' @param ... Additional arguments passed to \code{fastTopics::fit_topic_model()}.
#' @return A list containing the factor matrices \code{Fmat} and \code{Lmat}, differential expression results \code{de_res}, and p-values matrix \code{p_jk}.
#' @details
#' The function ensures that the count matrix has cells as rows and peaks as columns.
#' It then runs the fastTopics model, performs differential expression analysis,
#' computes local false discovery rates (lfdr), and calculates p-values for each peak-topic pair.
#' @examples
#' \dontrun{
#' # Assuming 'counts' is your count matrix
#' results <- run_fastTopics(counts, nTopics = 5)
#' }
#' @export
#' 
run_fastTopics <- function(count_matrix, nTopics = 10, n_s=1000, n_c=1, ...) {
  
  # Ensure that count_matrix has cells as rows and peaks as columns
  # Transpose if necessary
  if (nrow(count_matrix) > ncol(count_matrix)) {
    count_matrix <- Matrix::t(count_matrix)
  }
  print(dim(count_matrix))
  
  # # Convert to sparse matrix if not already
  # if (!inherits(count_matrix, "dgCMatrix")) {
  #   count_matrix_sp <- as(count_matrix, "dgCMatrix")
  # } else {
  #   count_matrix_sp <- count_matrix
  # }
  
  # Run fastTopics model
  # count_matrix should have cells as rows and peaks as columns
  # '...' allows passing additional arguments to the model fitting function
  cat("\nRunning fastTopics\n")
  cat("\nStart time: ")
  print(Sys.time())
  fastTopics_fit <- fastTopics::fit_topic_model(count_matrix, 
                                                k = nTopics, ...)
  
  cat("\nStop time: ")
  print(Sys.time())
  
  # Extract factor matrices
  Fmat <- fastTopics_fit$F  # Topics x peaks matrix (feature loadings)
  Lmat <- fastTopics_fit$L  # Cells x topics matrix (topic proportions)
  
  # # Return results as a list
  # return(list(
  #   Fmat = Fmat,
  #   Lmat = Lmat, 
  #   fit = fastTopics_fit
  # ))
  
  # Run differential expression analysis
  # 's' represents the total counts per cell
  s <- Matrix::rowSums(count_matrix)  # Vector of total counts per cell

  cat("\nRunning fastTopics-DE\n")
  cat("\nStart time: ")
  print(Sys.time())
  de_res <- de_analysis2(
    fit = fastTopics_fit,
    X = count_matrix,
    s = s,
    lfc.stat = "vsnull",
    shrink.method = "none",
    control = list(ns = n_s,
                   nc = n_c)
  )

  # 'de_res$z' is a matrix of z-scores with dimensions peaks x topics
  # Remove peaks with any NA values for z-scores in topics
  keep_peaks <- which(complete.cases(de_res$z))
  keep_z <- de_res$z[complete.cases(de_res$z), ]
  # Convert z-score matrix to a vector for locfdr input
  z_scores <- as.vector(keep_z)
  print(summary(z_scores))
  # hist(z_scores, breaks = 50, main = "Histogram of z_scores", xlab = "z_scores")
  
  # # Return results as a list
  # return(list(
  #   Fmat = Fmat,
  #   Lmat = Lmat,
  #   de_res = de_res
  # ))
  
  # Apply locfdr to compute local false discovery rates
  cat("\nApply locfdr\n")
  cat("\nStart time: ")
  print(Sys.time())
  locfdr_res <- locfdr::locfdr(z_scores, nulltype = 2, plot=0, df=30)  # 'plot = 0' suppresses plotting
  locfdr_vals <- locfdr_res$fdr  # Extract lfdr values

  # Reshape locfdr_vals back to a matrix with the same dimensions as de_res$z
  locfdr_mat <- matrix(
    locfdr_vals,
    nrow = length(keep_peaks),
    ncol = ncol(de_res$z),
    byrow = FALSE
  )

  # Compute p-values as 1 - lfdr for each peak-topic pair
  p_jk <- 1 - locfdr_mat  # Matrix of p-values (peaks x topics)
  
  # Adjust Fmat to match p_jk
  Fmat <- Fmat[keep_peaks, ]

  # Return results as a list
  return(list(
    Fmat = Fmat,
    Lmat = Lmat,
    de_res = de_res,
    Pmat = p_jk
  ))
}
