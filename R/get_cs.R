#' Compute Single-Cell Level Scores (c_i)
#'
#' This function computes the single-cell level scores (\eqn{c_i}) based on the fastTopics model outputs
#' and other required parameters.
#'
#' @param topic_res Results from run_fastTopics 
#'  - p_jk Matrix of probabilities of peaks in topics (peaks x topics), obtained from the F matrix of fastTopics.
#'  - l_ik Matrix of topic proportions for each cell (cells x topics), obtained from the L matrix of fastTopics.
#' @param ldsc_res Enrichment results from run_ldsc 
#'  - tau_Ck Vector of per-topic heritability contributions (\eqn{\tau_{C_k}}) of length K (number of topics).
#' @return A vector of single-cell level scores (\eqn{c_i}) of length I (number of cells).
#' @details
#' The function computes the scores based on the following corrected approximation:
#' \deqn{
#' c_i = \frac{ \left( \frac{E(M_i)}{E(N_i)} \left[ 1 - \frac{\text{Cov}(M_i, N_i)}{E(M_i) E(N_i)} + \frac{\text{Var}(N_i)}{[E(N_i)]^2} \right] \right) }{ \left( \frac{ \sum_{j=1}^{J} h_j }{ J } \right) }
#' }
#' where:
#' \itemize{
#'   \item \eqn{E(M_i)} and \eqn{E(N_i)} are the expected values.
#'   \item \eqn{\text{Cov}(M_i, N_i)} is the covariance between \eqn{M_i} and \eqn{N_i}.
#'   \item \eqn{\text{Var}(N_i)} is the variance of \eqn{N_i}.
#' }
#' @export
get_cs <- function(topic_res, ldsc_res, trait) {
  
  # Ensure that the input matrices and vectors have compatible dimensions
  # p_jk: J x K (peaks x topics)
  # l_ik: I x K (cells x topics)
  # tau_Ck: vector of length K
  
  p_jk <- topic_res$Pmat
  l_ik <- topic_res$Lmat
  
  # tau_res <- read.table(file.path(ldsc_res, sprintf('results/%s_enrichment.results', trait), header = T)
  # tau_res <- read.table("/dartfs-hpc/rc/home/b/f005d5b/Zhao-labshare/liyang/organoid_stimulus_proj/ATAC/Combined_samples_202405/filtered/bam2vcf/UT/ASoC_samp18/polyfun/ASoC_11samp_WGS_outputs/results/CRC_EUR_enrichment.results", header = T)
  tau_Ck <- c(1, 152, 1, 1, 1)
  
  # Number of peaks (J), topics (K), and cells (I)
  J <- nrow(p_jk)
  K <- ncol(p_jk)
  I <- nrow(l_ik)
  
  # Check dimensions
  if (ncol(l_ik) != K) {
    stop("Number of topics in l_ik and p_jk must be the same.")
  }
  if (length(tau_Ck) != K) {
    stop("Length of tau_Ck must be equal to the number of topics (K).")
  }
  
  # Since v_j is the number of variants and each peak corresponds to one variant,
  # we set v_j to a vector of ones.
  v_j <- rep(1, J)
  
  # 1. Compute h_j: heritability contribution from peak j
  # h_j = sum over k of (tau_Ck[k] * p_jk[j, k])
  h_j <- rowSums(sweep(p_jk, 2, tau_Ck, '*'))  # h_j is length J
  
  # 2. Compute E(I_{ij}) matrix: E_Iij = l_ik %*% t(p_jk)
  E_Iij <- l_ik %*% t(p_jk)  # Resulting matrix is I x J
  
  # 3. Compute E(M_i) and E(N_i)
  # E(M_i) = rowSums( E_Iij * h_j )
  E_Mi <- rowSums(E_Iij * matrix(h_j, nrow = I, ncol = J, byrow = TRUE))
  
  # E(N_i) = rowSums( E_Iij )
  E_Ni <- rowSums(E_Iij)
  
  # 4. Compute Var(I_{ij}) = E(I_{ij}) * (1 - E(I_{ij}))
  Var_Iij <- E_Iij * (1 - E_Iij)  # I x J matrix
  
  # 5. Compute Cov(M_i, N_i)
  # Cov(M_i, N_i) = rowSums( Var_Iij * h_j )
  Cov_Mi_Ni <- rowSums(Var_Iij * matrix(h_j, nrow = I, ncol = J, byrow = TRUE))
  
  # 6. Compute Var(N_i)
  # Var(N_i) = rowSums( Var_Iij )
  Var_Ni <- rowSums(Var_Iij)
  
  # 7. Compute Numerator
  Numerator <- (E_Mi / E_Ni) * (1 - (Cov_Mi_Ni / (E_Mi * E_Ni)) + (Var_Ni / (E_Ni^2)))
  
  # 8. Compute Denominator
  Denominator <- sum(h_j) / J  # Since sum(v_j) = J
  
  # 9. Compute c_i
  c_i <- Numerator / Denominator  # Vector of length I (cells)
  
  # Return the single-cell scores
  return(c_i)
}
