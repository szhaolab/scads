
#'  Calculate GC Baseline Correction Parameters for Topic-Specific Count Data 
#'
#' @description
#' This function is the main function that computes GC-content baseline correction parameters (lambda values)
#' for each genomic region and topic combination. It processes count matrices with
#' topic modeling results, estimates GC effects for each topic separately, and
#' returns normalized baseline parameters that can be used for GC bias correction.
#'
#' @param count_matrix A SummarizedExperiment or similar object containing a count
#'   matrix accessible via \code{assay()}. Rows represent genomic peaks and columns
#'   represent cells or samples.
#' @param Lmat A matrix of topic loadings where rows correspond to cells (matching
#'   column names in \code{count_matrix}) and columns represent topics.
#' @param genome A character string specifying the genome build. Must be a valid
#'   BSgenome package name (e.g., "hg19", "hg38", "mm10"). Default is "hg19".
#' @param peakwidth An integer specifying the width of genomic peaks in base pairs.
#'   Default is 501.
#' @param emtrace Logical; if \code{TRUE}, prints iteration details during EM
#'   algorithm convergence for each topic. Default is \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, prints progress messages during GC
#'   content calculation and data preparation. Default is \code{FALSE}.
#' @param plot Logical; if \code{TRUE}, generates diagnostic plots showing GC
#'   effects for each topic. Default is \code{FALSE}.
#' @param gcrange A numeric vector of length 2 specifying the GC content range
#'   to include in model fitting. Default is \code{c(0.3, 0.8)}.
#' @param mu0,mu1 Numeric initial values for background and foreground mean read
#'   counts. Defaults are \code{mu0 = 1} and \code{mu1 = 50}.
#' @param theta0,theta1 Numeric initial shape parameters for negative binomial
#'   distribution. Default to \code{mu0} and \code{mu1} respectively.
#' @param p Numeric initial mixture proportion for foreground regions. Default is 0.02.
#' @param converge Numeric convergence threshold for EM algorithm. Default is 1e-3.
#' @param max_pts Integer specifying maximum number of points to plot. Default is 100000.
#' @param max_line Integer specifying maximum number of points for fitted curve lines.
#'   Default is 5000.
#'
#' @return A list containing:
#' \describe{
#'   \item{lambda_jk}{A data frame where rows represent genomic regions (peaks) and
#'     columns represent topics. Each entry is the normalized baseline parameter
#'     (lambda) for that region-topic combination.}
#'   \item{N_k}{A named numeric vector containing the total read count for each topic.}
#'   \item{gc_k}{A named list where each element contains the GC content vector for
#'     all regions, organized by topic.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Multiplies the count matrix by the topic loading matrix to obtain
#'     topic-specific read counts.
#'   \item Calculates GC content for each genomic peak region.
#'   \item Fits adaptive GC effect models for each topic using negative binomial
#'     regression with an EM algorithm.
#'   \item Predicts baseline mean values (mu0) for each region and topic.
#'   \item Normalizes by total topic counts to obtain lambda parameters.
#' }
#'
#' The lambda values can be used to correct for GC bias in downstream analyses
#' such as differential accessibility or peak calling.
#'
#' @seealso \code{\link{prep_gcEffects_input}}, \code{\link{adp_gcEffects}},
#'   \code{\link{pred_baseline_mu0}}
#'
#' @importFrom SummarizedExperiment assay
#'
#' @examples
#' \dontrun{
#'   gc_baseline_res <- get_gc_baseline(count_matrix = count_mat_raw, Lmat = L_mat_topic)
#' }
#'
#' @export

get_gc_baseline <- function(count_matrix, Lmat, Fmat, outdir,
                             genome="hg19",
                             peakwidth = 501,
                             emtrace = FALSE,
                             verbose = FALSE,
                             plot = FALSE,
                             # gctype=c("ladder","tricube"),
                             gcrange=c(0.3,0.8),
                             # model=c('nbinom','poisson'),
                             mu0=1,mu1=50,theta0=mu0,theta1=mu1,
                             p=0.02,converge=1e-3,
                             max_pts = 100000,
                             max_line = 5000) {
  
  # get count_mat_tp matrix --> matrix of read count: peak x topic
  # count_mat <- assay(count_matrix)
  # peak_names <- count_mat@Dimnames[[1]]
  if (is(count_matrix, "SummarizedExperiment")) {
    count_mat <- assay(count_matrix)
  } else {
    count_mat <- count_matrix
  }
  
  # if (nrow(count_mat) > ncol(count_mat)) {
  #   count_mat <- Matrix::t(count_mat)
  # }
  message("Count matrix dimensions: ", nrow(count_mat), " x ", ncol(count_mat))

  # Get peak names from columns (peaks are in columns)
  peak_names <- colnames(count_mat)
  # Only convert peak names if necessary (check format first)
  if (any(grepl(":", peak_names))) {
    # If using colon format like "chr1:10234-10734"
    peak_names <- gsub(":", "_", peak_names)
    peak_names <- gsub("-", "_", peak_names)
  }
  # Otherwise peak_names are already in correct format: e.g. "chr1_10234_10734"
  
  # Check for NAs in peak names
  if (any(is.na(peak_names))) {
    stop("Found NA values in peak names. Check colnames of count_matrix.")
  }
  
  # Select cells in count matrix that only appear in matrix L
  # unique_cells_in_L <- rownames(Lmat)
  # Make sure dimensions align
  # count_mat should be: cells (rows) × peaks (columns)
  # Lmat should be: cells (rows) × topics (columns)
  # if (!all(unique_cells_in_L %in% rownames(count_mat))) {
  #   stop("Not all cells in Lmat are present in count_matrix rownames")
  # }
  # count_mat_in_L <- count_mat[unique_cells_in_L, , drop=FALSE]
  # Matrix multiplication: (peaks × cells)^T %*% (cells × topics) = (peaks × topics)
  # count_mat_tp <- Matrix::t(count_mat_in_L) %*% Lmat
  
  # ## Obtain pseudobulk (memory usage is massive)
  # rate_matrix <- tcrossprod(Lmat, Fmat)
  # # c_k <- colSums(count_mat * tcrossprod(Lmat[, k, drop=F], F_scaled[,k , drop = F])/rate_matrix)
  # weights <- count_mat / rate_matrix
  # rm(rate_matrix)
  # count_mat_tp <- Fmat * (Matrix::t(weights) %*% Lmat)
  # saveRDS(count_mat_tp, file.path(outdir, "count_mat_tp.rds"))
  
  ## Obtain pseudobulk (save memory)
  # 1. Extract the non-zero coordinates from count_mat
  # summary() on a sparse matrix gives a dataframe with:
  # i (row index), j (col index), x (value)
  cat("\nGet Triplets")
  triplets <- Matrix::summary(as(count_mat, "dgCMatrix"))
  # 2. Compute the 'rate' ONLY at these non-zero coordinates
  # Instead of a full matrix multiplication, we do a "lookup and dot product"
  # We extract the specific rows of L and F corresponding to each non-zero count
  # L_subset: (N_nonzero x K)
  # F_subset: (N_nonzero x K)
  cat("\nGet L_subset")
  L_subset <- Lmat[triplets$i, , drop = FALSE] 
  cat("\nGet F_subset")
  F_subset <- Fmat[triplets$j, , drop = FALSE] 
  # Compute dot product for each pair (equivalent to the specific cell in tcrossprod)
  cat("\nCalc rate values")
  rate_values <- rowSums(L_subset * F_subset)
  # 3. Perform the division just for these values
  # w_ij = count_ij / rate_ij
  cat("\nCalc weight values")
  weight_values <- triplets$x / rate_values
  # 4. Reconstruct the sparse 'weights' matrix
  # We create a new sparse matrix using the original indices but the new values
  weights <- Matrix::sparseMatrix(
    i = triplets$i,
    j = triplets$j,
    x = weight_values,
    dims = dim(count_mat)
  )
  # 6. Final Calculation
  cat("\nGet pseudobulk count matrix")
  count_mat_tp <- Fmat * (Matrix::t(weights) %*% Lmat)
  saveRDS(count_mat_tp, file.path(outdir, "count_mat_tp.rds"))
  rm(weights, weight_values, rate_values, F_subset, L_subset, count_mat)
  
  # prepare input
  gcEffects_input_res_topic_list <- list()
  for (k in 1:ncol(count_mat_tp)) {
    gcEffects_input_res_topic <- prep_gcEffects_input(rc = count_mat_tp[,k], 
                                                      peak_names = peak_names,
                                                      gctype = "ladder",
                                                      genome = genome,
                                                      peakwidth = peakwidth,
                                                      verbose = verbose)
    gcEffects_input_res_topic_list[[k]] <- gcEffects_input_res_topic
  }
  cat("...... GC content correction inputs are prepared. \n")
  
  # get adaptive gc effects
  adp_gcEffects_res_topic_list <- list()
  for (k in 1:length(gcEffects_input_res_topic_list)) {
    cat("...... Estimating GC effects for topic ", colnames(count_mat_tp)[k],"\n")
    cur_gc_res <- gcEffects_input_res_topic_list[[k]]
    # print(cur_gc_res)
    gc <- cur_gc_res$gc
    # print(gc)
    region <- cur_gc_res$region
    # print(region)
    rc <- cur_gc_res$rc
    # print(rc)
    
    # Add these diagnostic prints:
  cat("Topic ", k, "\n")
  cat("......... GC summary: ", summary(cur_gc_res$gc), "\n")
  cat("......... NA count in GC: ", sum(is.na(cur_gc_res$gc)), "\n")
  cat("......... NA count in RC: ", sum(is.na(cur_gc_res$rc)), "\n")
  cat("......... Inf count in GC: ", sum(is.infinite(cur_gc_res$gc)), "\n")
  cat("......... Inf count in RC: ", sum(is.infinite(cur_gc_res$rc)), "\n")
    
    adp_gcEffects_res_topic <- adp_gcEffects(gc=gc,
                                             region=region,
                                             rc=rc,
                                             model=('poisson'), # model=('nbinom'),
                                             emtrace = emtrace,
                                             plot = plot,
                                             gcrange = gcrange,
                                             peakwidth = peakwidth,
                                             mu0=mu0,mu1=mu1,theta0=theta0,theta1=theta1,
                                             p=p,converge=converge,
                                             max_pts = max_pts,
                                             max_line = max_line)
    cat("...... Plotted Topic ", k, "\n")
    adp_gcEffects_res_topic_list[[k]] <- adp_gcEffects_res_topic
  }
  cat("...... Baseline models for are fitted. \n")
  
  # get predicted baseline mu0 for all regions
  pred_mu0_list <- list()
  pred_mu1_list <- list()
  for (k in 1:length(gcEffects_input_res_topic_list)) {
    pred_mu0_list[[k]] <- pred_baseline_mu0(adp_gcEffects_res_topic_list[[k]], 
                                            gcEffects_input_res_topic_list[[k]]$gc)
    pred_mu1_list[[k]] <- predict(adp_gcEffects_res_topic_list[[k]]$lmns1, 
                                  data.frame(gc = gcEffects_input_res_topic_list[[k]]$gc), 
                                  type = "response")
  }
  names(pred_mu0_list) <- colnames(count_mat_tp)
  names(pred_mu1_list) <- colnames(count_mat_tp)
  
  # get lambda
  lambda_j_list <- list()
  N_k_list <- c()
  gc_list <- c()
  # rc_list <- c()
  for (k in 1:length(pred_mu0_list)) {
    N_k <- sum(count_mat_tp[,k])
    lambda_j_list[[k]] <- pred_mu0_list[[k]] / N_k
    ## N_k_list <- c(N_k, N_k_list) ## indexing is reversed - k5, k4, k3, k2, k1
    N_k_list <- c(N_k_list, N_k)
    gc_list[[k]] <- gcEffects_input_res_topic_list[[k]]$gc
    # rc_list[[k]] <- gcEffects_input_res_topic_list[[k]]$rc
  }
  names(lambda_j_list) <- colnames(count_mat_tp)
  names(N_k_list) <- colnames(count_mat_tp)
  names(gc_list) <- colnames(count_mat_tp)
  lambda_jk <- as.data.frame(do.call(cbind, lambda_j_list))
  rownames(lambda_jk) <- peak_names
  
  rm(count_mat_tp)
  return(list(
    lambda_jk = lambda_jk,
    mu0_jk = as.data.frame(do.call(cbind, pred_mu0_list)),
    mu1_jk = as.data.frame(do.call(cbind, pred_mu1_list)),
    N_k = N_k_list,
    gc_k = gc_list
    # rc_k = rc_list
  ))
}


#' Prepare GC content input for GC-effect correction
#'
#' This function computes the GC content across genomic peak regions, applying a weighting
#' scheme (either uniform "ladder" or smoothed "tricube") over each region of fixed width.
#' It prepares a list of GC content, genomic regions, and the provided read count information
#' for downstream GC-effect correction or modeling.
#'
#' @param rc A numeric vector, matrix, or data frame representing read counts (e.g., from ATAC-seq or ChIP-seq).
#'   This can be sparse or dense, and should align with the provided \code{peak_names}.
#' @param peak_names A character vector of peak names, typically in the format
#'   \code{"chr_start_end"}, e.g., \code{"chr1_12345_12845"}.
#' @param genome A character string specifying the genome build to use. Must be a valid
#'   BSgenome package name (e.g., \code{"hg19"}, \code{"hg38"}, \code{"mm10"}). Default is \code{"hg19"}.
#' @param peakwidth An integer specifying the width of each peak (default: \code{501}).
#'   Used to determine the smoothing window for GC content calculation.
#' @param gctype A character string specifying the GC weighting scheme.
#'   Options are \code{"ladder"} (uniform weights) or \code{"tricube"} (smooth kernel weights).
#'   Default is \code{"ladder"}.
#'
#' @details
#' The function extracts sequences for the given genomic regions using \code{getSeq()}
#' from a BSgenome object and computes the GC content in each region by summing weights
#' at positions corresponding to G or C nucleotides. The result is scaled to sum to 1.
#'
#' @return A list containing:
#' \describe{
#'   \item{gc}{A numeric vector of GC content values for each region.}
#'   \item{region}{A \code{GRanges} object representing the genomic regions.}
#'   \item{rc}{The input read count data, returned for convenience.}
#' }
#'
#' @importFrom BSgenome getBSgenome
#' @importFrom Biostrings getSeq vmatchPattern startIndex
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' peaks <- c("chr1_100000_100500", "chr2_200000_200500")
#' rc <- c(50, 75)
#' res <- prep_gcEffects_input(rc = rc, peak_names = peaks, genome = "hg19", gctype = "tricube")
#' head(res$gc)
#' }
#'
#' @export

prep_gcEffects_input <- function(rc, 
                                 # read_count_mat_sparse,
                                 peak_names,
                                 genome="hg19",
                                 peakwidth = 501,
                                 gctype=c("ladder","tricube"),
                                 verbose = TRUE){
  ### input sanity check
  genome <- getBSgenome(genome)
  gctype <- match.arg(gctype)
  halfbdw <- floor(peakwidth/2)
  if(gctype=="ladder"){
    weight <- c(rep(1,peakwidth))
    weight <- weight/sum(weight)
  }else if(gctype=="tricube"){
    w <- halfbdw
    weight <- (1-abs(seq(-w,w)/w)^3)^3
    weight <- weight/sum(weight)
  }
  
  ### prepare peak regions
  parts  <- do.call(rbind, strsplit(peak_names, "_", fixed = TRUE))
  chrs   <- parts[, 1]
  starts <- as.integer(parts[, 2])
  ends   <- as.integer(parts[, 3])

  region <- GRanges(chrs,IRanges(start=starts,end=ends))
  
  ### effective gc content1
  if (verbose) {
    cat("...... Calculating GC content", "\n")
  }
  nr <- region
  seqs <- getSeq(genome,nr) # slow
  gcpos <- startIndex(vmatchPattern("S", seqs, fixed="subject")) # a list of the positions within each window that are G or C
  
  # gc <- round(sapply(gcpos,function(x) sum(weight[x])),3) # commented out 3/9/2026
  # Compute GC per-peak using actual sequence length (handles variable-width peaks)
  gc <- round(sapply(seq_along(seqs), function(i) {
    seq_len <- width(seqs)[i]
    w <- rep(1, seq_len) / seq_len   # uniform weight over actual peak length
    pos <- gcpos[[i]]
    pos_valid <- pos[pos <= seq_len]  # safety guard (should always be true)
    sum(w[pos_valid])
  }), 3)
  
  ### round all values in rc to integers
  if (verbose) {
    cat("...... Rounding read counts to nearest integers", "\n")
  }
  rc <- round(rc)
  
  return(list(gc=gc,
              region=region,
              rc=rc))
}



#' @title Adaptive GC Effects Estimation for ChIP-seq Read Counts
#'
#' @description
#' Estimates GC-content-dependent effects on ChIP-seq read counts using an
#' expectation-maximization (EM) algorithm with generalized linear models.
#' Regions are modeled as a mixture of background and foreground components
#' (e.g. non-peak and peak-like regions), each following either a Poisson or
#' negative binomial distribution. GC effects are fitted separately for each
#' component using natural spline regression. Optional visualization shows
#' fitted curves and mixture probabilities.
#'
#' @param gc A numeric vector of GC content values for genomic windows or regions.
#'
#' @param region A genomic region object or identifier (not directly used in
#' modeling but retained for consistency or downstream reference).
#'
#' @param rc A numeric vector of read counts corresponding to the same
#' regions as \code{gc}.
#'
#' @param peakwidth A non-negative integer specifying the expected ChIP-seq
#' binding width (default 501). Used for labeling and downstream interpretation.
#'
#' @param plot Logical; if \code{TRUE} (default), generates a scatter plot
#' showing GC content vs read counts, colored by posterior mixture
#' probabilities, and overlays fitted foreground (red) and background (blue)
#' curves.
#'
#' @param gcrange A numeric vector of length 2 specifying the GC-content range
#' to include for model fitting. Regions outside this range are ignored.
#' Default is \code{c(0.3, 0.8)}.
#'
#' @param emtrace Logical; if \code{TRUE} (default), prints log-likelihood and
#' convergence progress at each EM iteration.
#'
#' @param model Character string specifying the distribution model for read
#' counts. Supported options are \code{"nbinom"} (default) and \code{"poisson"}.
#'
#' @param mu0,mu1 Numeric values giving the initial mean read counts for
#' background and foreground components, respectively. Defaults are
#' \code{mu0 = 1}, \code{mu1 = 50}.
#'
#' @param theta0,theta1 Numeric values specifying initial shape parameters
#' for negative binomial models of background and foreground. Used only when
#' \code{model = "nbinom"}.
#'
#' @param p Numeric value specifying the initial mixture proportion of
#' foreground regions. Default is 0.02.
#'
#' @param converge Numeric value specifying the EM convergence threshold.
#' Iteration stops when the relative log-likelihood change is below this
#' threshold. Default is 1e-3.
#'
#' @param max_pts Integer specifying the maximum number of points plotted
#' (default 100000).
#'
#' @param max_line Integer specifying the maximum number of points used to
#' draw fitted curves (default 5000).
#'
#' @return A list containing:
#' \item{gc}{GC-content values at which GC effects were estimated.}
#' \item{lmns0}{Fitted GLM object for the background component.}
#' \item{lmns1}{Fitted GLM object for the foreground component.}
#' \item{z}{Posterior probabilities of each region belonging to the foreground component.}
#' \item{mu0, mu1}{Predicted read counts (fitted means) at each GC content for
#' background and foreground components, respectively.}
#' \item{mu0med0, mu1med1}{Medians of fitted background and foreground signals
#' in their respective component regions.}
#' \item{mu0med1, mu1med0}{Cross-median values: background prediction in
#' foreground-like regions, and vice versa.}
#' \item{g}{A \code{ggplot} object (if \code{plot=TRUE}) showing fitted GC effects.}
#' \item{dat}{Data frame containing the filtered counts and GC content used for fitting.}
#'
#' @details
#' The algorithm iteratively updates posterior probabilities (E step) and fits
#' GC-dependent regression models for both mixture components (M step). The GC
#' dependence is modeled using a natural spline with 2 degrees of freedom.
#' 
#' The fitted curves can reveal whether sequencing depth or read coverage is
#' biased by GC content, separately for peak-enriched and background regions.
#'
#' @import MASS
#' @importFrom splines ns
#' @importFrom stats glm dpois dnbinom predict
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' gc <- runif(10000, 0.2, 0.9)
#' rc <- rpois(10000, lambda = exp(5 * (gc - 0.5)^2))
#' res <- adp_gcEffects(gc = gc, rc = rc, model = "poisson", plot = TRUE)
#' 
#' # visualize the estimated GC effects
#' print(res$g)
#' }
#' 
#' 
adp_gcEffects <- function(gc=gc,
                          region=region,
                          rc=rc,
                          # coverage,bdwidth,flank=NULL,
                          peakwidth = 501, 
                          plot=TRUE,
                      # sampling=c(0.05,1),#supervise=GRanges(),
                      gcrange=c(0.3,0.8),emtrace=TRUE,
                      model=c('nbinom','poisson'),
                      mu0=1,mu1=50,theta0=mu0,theta1=mu1,
                      p=0.02,converge=1e-3,
                      max_pts = 100000,
                      max_line = 5000
                      # genome="hg19",gctype=c("ladder","tricube")
                      ){
  
  ### em algorithms
  # gc <- rep(gc,2)
  # rc <- c(rcfwd,rcrev)
  
  # browser()
  idx <- gc>=gcrange[1] & gc<=gcrange[2] & !is.na(rc) & !is.na(gc) #QX
  # idx <- !is.na(gc) & !is.na(rc) & gc >= gcrange[1] & gc <= gcrange[2] & is.finite(gc) & is.finite(rc) #LY
  dat <- data.frame(y=rc[idx], gc=gc[idx]) 
  
  # Add one safety check:
  if(nrow(dat) < 100) {
    stop("Insufficient data points after filtering (n=", nrow(dat), 
         "). Check your data quality and gcrange parameter.")
  }
  
  # dat0 <- data.frame(y=rc,gc=gc)
  if(model=='poisson'){
    logp1 <- dpois(dat$y, lambda = mu1, log = TRUE)
    logp0 <- dpois(dat$y, lambda = mu0, log = TRUE)
  }else{
    # browser()
    logp1 <- dnbinom(dat$y, size=theta1, mu=mu1, log = TRUE)
    logp0 <- dnbinom(dat$y, size=theta0, mu=mu0, log = TRUE)
  }
  z <- 1/(1+exp(logp0-logp1)*(1-p)/p)
  llf <- sum(z*(logp1+log(p)) + (1-z)*(logp0+log(1-p)))
  llgap <- llf
  i <- 0
  # rm(rcfwd,rcrev,rc,gc,idx)
  while(abs(llgap) > (abs(llf) * converge) && i < 100){
    p <- (2+sum(z))/(2*2+length(z))
    dat1 <- dat[z>=0.5,]
    dat0 <- dat[z<0.5,]
    if(model=='poisson'){
      lmns0 <- glm(y ~ ns(gc, df = 2), data=dat0, family="poisson")
      lmns1 <- glm(y ~ ns(gc, df = 2), data=dat1, family="poisson")
      predY0 <- predict(lmns0, data.frame(gc = dat$gc),type="response")
      predY1 <- predict(lmns1, data.frame(gc = dat$gc),type="response")
      logp1 <- dpois(dat$y, lambda = predY1, log = TRUE)
      logp0 <- dpois(dat$y, lambda = predY0, log = TRUE)
    }else{
      # slow
      lmns0 <- glm.nb(y ~ ns(gc, df = 2), data=dat0, init.theta=theta0)
      lmns1 <- glm.nb(y ~ ns(gc, df = 2), data=dat1, init.theta=theta1)
      predY0 <- predict(lmns0, data.frame(gc = dat$gc),type="response")
      predY1 <- predict(lmns1, data.frame(gc = dat$gc),type="response")
      theta1 <- lmns1$theta
      theta0 <- lmns0$theta
      logp1 <- dnbinom(dat$y, size=theta1, mu=predY1, log = TRUE)
      logp0 <- dnbinom(dat$y, size=theta0, mu=predY0, log = TRUE)
    }
    z <- 1/(1+exp(logp0-logp1)*(1-p)/p)
    if(sum(z>=0.5) < length(gc)*0.0005 | sum(z<0.5) < length(gc)*0.0005)
      break;
    lli <- sum(z*(logp1+log(p)) + (1-z)*(logp0+log(1-p)))
    llgap <- lli - llf
    llf <- lli
    i <- i + 1
    if(emtrace)
      cat("......... Iteration",i,'\tll',llf,'\tincrement',llgap,'\n')
  }
  
  samptype <- ""
  if(plot){
    tmp <- nrow(dat)
    idx0 <- sample.int(tmp,min(100000,tmp))
    rbPal <- colorRampPalette(c('skyblue','pink'))
    color <- rbPal(20)[as.numeric(cut(z[idx0],breaks = 20))]
    pdf(file=paste0(outdir, sprintf("/GC_plot_%s.pdf", idx0)))
    plot(dat$gc[idx0],dat$y[idx0]+0.5,col=color,xlim=gcrange,pch=20,
         main=paste(model,samptype),
         xlab='Effective GC content',ylab="Read counts",log='y',yaxt='n')
    idx00 <- sample.int(tmp,min(5000,tmp))
    idx00 <- idx00[order(dat$gc[idx00])]
    lines(dat$gc[idx00],predY1[idx00]+0.5,col='red',lwd=3)
    lines(dat$gc[idx00],predY0[idx00]+0.5,col='blue',lwd=3)
    axis(side=2, at=c(0,2^(0:10))+0.5, labels=c(0,2^(0:10)))
    dev.off()
  }
  
  ### gc effects
  gcbase <- round(seq(0,1,0.001),3)
  gcbias <- list(gc=gcbase,
                 lmns0 = lmns0, # add the predicted model to the output
                 lmns1 = lmns1,
                 z = z,
                 mu0=predict(lmns0,data.frame(gc = gcbase),type="response"),
                 mu1=predict(lmns1,data.frame(gc = gcbase),type="response"),
                 mu0med0=median(predY0[z<0.5]),
                 mu1med1=median(predY1[z>=0.5]),
                 mu0med1=median(predY0[z>=0.5]),
                 mu1med0=median(predY1[z<0.5]),
                 # g = g,
                 dat = dat)
  
  gcbias
}

#' Predict Baseline Mean Values (mu0) from GC Content
#'
#' This function predicts the baseline (background component) mean read counts
#' for specified GC content values using a fitted model from \code{adp_gcEffects()}.
#' It extracts the background GLM model (lmns0) and generates predictions for
#' the provided GC content vector.
#'
#' @param adp_gcEffects_res A list object returned by \code{adp_gcEffects()},
#'   which contains fitted GLM models including \code{lmns0} (the background model).
#' @param gc A numeric vector of GC content values (between 0 and 1) for which
#'   baseline predictions are desired. (The processed gc content information from \code{prep_gcEffects_input()}.)
#'
#' @return A numeric vector of predicted baseline mean values (mu0) corresponding
#'   to each element in \code{gc}.
#'
#' @details
#' This is a helper function typically used after fitting GC effects with
#' \code{adp_gcEffects()}. It uses the background component model (lmns0) to
#' predict expected read counts at given GC content levels, which can be used
#' for normalization or bias correction.
#'
#' @seealso \code{\link{adp_gcEffects}}, \code{\link{get_gc_baseline}}
#'
#' @export
#' 
pred_baseline_mu0 <- function(adp_gcEffects_res, gc) {
  # Predict baseline mu0 values for given GC content values
  pred_mu0 <- predict(adp_gcEffects_res$lmns0, 
                           data.frame(gc = gc), 
                           type = "response")
  return(pred_mu0)
}




