#' Compute Single-Cell Level Scores (c_i)
#'
#' @param topic_res Results from run_fastTopics, with Pmat (J×K) and Lmat (I×K)
#' @param ldsc_res_dir Directory containing k*_output/results/Trait.results
#' @param trait Trait name used in LDSC results files.
#' @param nTopics Number of topics (K).
#' @return A list with
#'   - cs: vector of length I (cells)
#'   - var_cs: variance of each cs[i]
#'   - ldsc_res_table: the raw tau_display table from LDSC
#' @export
get_cs <- function(topic_res, ldsc_res_dir, trait, nTopics) {
  
  # 1) unpack
  p_jk <- topic_res$Pmat   # J x K
  l_ik <- topic_res$Lmat   # I x K
  
  # 2) read in LDSC results per topic
  tau        <- vector("list", nTopics)
  tau_display <- NULL
  for(k in seq_len(nTopics)){
    res_f <- file.path(ldsc_res_dir,
                       paste0("k", k, "_output"),
                       "results",
                       paste0(trait, ".results"))
    if(!file.exists(res_f)){
      warning("Missing ", res_f)
      tau[[k]] <- NA
      next
    }
    res <- read.table(res_f, header=TRUE, sep="\t", check.names=FALSE)
    if(nrow(res)<1){
      warning("Empty ", res_f)
      tau[[k]] <- NA
      next
    }
    res$Category <- paste0("k", k)
    tau[[k]] <- res$Prop._h2[1]
    tau_display <- rbind(tau_display, res[1,])
  }
  tau_Ck <- unlist(tau)
  e_Ck   <- tau_display$Enrichment
  se_e   <- tau_display$Enrichment_std_error
  
  # 3) per-topic peak counts
  a_k <- colSums(p_jk)  # length K
  
  # 4) compute M_i and N_i for each cell
  #    M_i = Σ_k  l_ik[i,k] * (a_k * e_Ck)[k]
  #    N_i = Σ_k  l_ik[i,k] * a_k[k]
  M_i <- as.numeric(l_ik %*% (a_k * e_Ck))
  N_i <- as.numeric(l_ik %*% a_k)
  
  # 5) single-cell score
  cs <- M_i / N_i
  
  # 6) propagate SE(e_Ck) via delta‐method:
  var_e   <- se_e^2
  #    Var(M_i) = Σ_k [l_ik[i,k]*a_k[k]]^2 * var_e[k]
  V_Mi <- rowSums( (l_ik * a_k)^2 *
                     matrix(var_e,
                            nrow=nrow(l_ik),
                            ncol=length(var_e),
                            byrow=TRUE) )
  var_cs <- V_Mi / (N_i^2)
  
  return(list(
    cs            = cs,
    var_cs        = var_cs,
    ldsc_res_table = tau_display,
    cs_dat = list(M_i=M_i, N_i=N_i, Var_Mi = V_Mi))
    
  )
}


#' #' Compute Single-Cell Level Scores (c_i)
#' #'
#' #' This function computes the single-cell level scores (\eqn{c_i}) based on the fastTopics model outputs
#' #' and other required parameters.
#' #'
#' #' @param topic_res Results from run_fastTopics
#' #'  - Pmat: Matrix of probabilities of peaks in topics (peaks x topics), obtained from the F matrix of fastTopics.
#' #'  - Lmat: Matrix of topic proportions for each cell (cells x topics), obtained from the L matrix of fastTopics.
#' #' @param ldsc_res_dir Directory containing folders for each topic with enrichment results from run_ldsc
#' #'  - tau_Ck: Vector of per-topic heritability contributions (\eqn{\tau_{C_k}}) of length K (number of topics).
#' #' @param trait Trait name used in LDSC results files.
#' #' @param nTopics Number of topics in the model.
#' #' @return A list containing:
#' #'   - cs: A vector of single-cell level scores (\eqn{c_i}) of length I (number of cells).
#' #'   - ldsc_res_table: A data frame of LDSC results for each topic.
#' #' @details
#' #' The function computes the scores based on the corrected approximation.
#' #' @export
#' get_cs <- function(topic_res, ldsc_res_dir, trait, nTopics) {
#'   
#'   # Extract matrices from topic_res
#'   p_jk <- topic_res$Pmat  # J x K matrix
#'   l_ik <- topic_res$Lmat  # I x K matrix
#'   
#'   # Initialize an empty list to store per-snp heritability (tau) values
#'   tau <- list()
#'   tau_display <- list()
#'   
#'   # Loop over each topic
#'   for (k in 1:nTopics) {
#'     # Construct the folder name for the current topic
#'     topic_folder <- paste0("k", k, "_output")
#'     
#'     # Construct the full path to the .results file
#'     results_file <- file.path(ldsc_res_dir, topic_folder, "results", sprintf("%s.results", trait))
#'     
#'     # Check if the file exists
#'     if (file.exists(results_file)) {
#'       # Read the .results file
#'       res <- read.table(
#'         results_file,
#'         header = TRUE,
#'         sep = "\t",
#'         stringsAsFactors = FALSE,
#'         check.names = FALSE
#'       )
#'       res$Category <- paste0("k", k)
#'       
#'       # Ensure that the data has at least one row
#'       if (nrow(res) >= 1) {
#'         
#'         # Extract the tau value (per-snp h2 for annot) for the topic annotation
#'         # Store the enrichment value in the tau list
#'         tau[[k]] <- res$Prop._h2[1]
#'         
#'         # Store the LDSC results for display
#'         if (is.null(tau_display)) {
#'           tau_display <- res[1, ]
#'         } else {
#'           tau_display <- rbind(tau_display, res[1, ])
#'         }
#'         
#'       } else {
#'         warning(paste("File", results_file, "does not have at least one row. Skipping."))
#'         tau[[k]] <- NA  # Assign NA if there are not enough rows
#'       }
#'     } else {
#'       warning(paste("File", results_file, "does not exist. Skipping."))
#'       tau[[k]] <- NA  # Assign NA if the file does not exist
#'     }
#'   }
#'   
#'   # Convert tau to a numeric vector
#'   tau_Ck <- unlist(tau)
#'   
#'   # Print the tau vector
#'   print(tau_display)
#'   
#'   # Number of peaks (J), topics (K), and cells (I)
#'   J <- nrow(p_jk)
#'   K <- ncol(p_jk)
#'   I <- nrow(l_ik)
#'   
#'   # Check dimensions
#'   if (ncol(l_ik) != K) {
#'     stop("Number of topics in l_ik and p_jk must be the same.")
#'   }
#'   if (length(tau_Ck) != K) {
#'     stop("Length of tau_Ck must be equal to the number of topics (K).")
#'   }
#'   
#'   # Since v_j is the number of variants and each peak corresponds to one variant,
#'   # we set v_j to a vector of ones.
#'   v_j <- rep(1, J)
#'   
#'   # get enrichment per topic
#'   e_Ck <- c(tau_display$Enrichment)
#'   
#'   # Step 1: Compute variants in open peaks
#'   v_j_p_jk <- p_jk * v_j  # Matrix of m by k
#'   
#'   # Step 2: Weight enrichment by variants in open peaks
#'   v_j_p_jk_e_Ck <- e_Ck * t(v_j_p_jk) # Matrix of k by m
#'   
#'   # Step 3: Compute E(M_i)
#'   E_Mi <- rowSums(l_ik %*% v_j_p_jk_e_Ck)  # Vector of length I
#'   
#'   # Step 4: Compute E(N_i)
#'   E_Ni <- rowSums(l_ik %*% t(v_j_p_jk))  # Vector of length I
#'   
#'   # Step 5: Compute Cov(M_i, N_i)
#'   S_k_Cov1 <- colSums(p_jk * tau_Ck * v_j)
#'   First_Term_Cov_Mi_Ni <- l_ik %*% S_k_Cov1
#'   S_kk_Cov2 <- t(p_jk * tau_Ck * v_j) %*% p_jk
#'   Second_Term_Cov_Mi_Ni <- rowSums((l_ik %*% S_kk_Cov2) * l_ik)
#'   Cov_Mi_Ni <- First_Term_Cov_Mi_Ni - Second_Term_Cov_Mi_Ni
#'   
#'   # p_ij <- l_ik %*% t(p_jk)
#'   # Cov_Mi_Ni <- h_j * v_j * p_ij * (1 - p_ij)
#'   
#'   # Step 6: Compute Var(N_i)
#'   S_k_N <- colSums(p_jk * (v_j^2))
#'   First_Term_Var_Ni <- l_ik %*% S_k_N
#'   S_kk_Var_Ni <- t(p_jk * (v_j^2)) %*% p_jk
#'   Second_Term_Var_Ni <- rowSums((l_ik %*% S_kk_Var_Ni) * l_ik)
#'   Var_Ni <- First_Term_Var_Ni - Second_Term_Var_Ni
#'   
#'   # Step 7: Compute Var(M_i)
#'   S_k_M1 <- colSums(p_jk * tau_Ck * e_Ck * v_j)
#'   First_Term_Var_Mi <- l_ik %*% S_k_M1
#'   S_kk_M2 <- t(p_jk * tau_Ck * e_Ck * v_j) %*% p_jk
#'   Second_Term_Var_Mi <- rowSums((l_ik %*% S_kk_M2) * l_ik)
#'   Var_Mi <- First_Term_Var_Mi - Second_Term_Var_Mi
#'   
#'   # Step 8: Compute cs = E(M_i / N_i)
#'   cs <- (E_Mi / E_Ni) * (1 - (Cov_Mi_Ni / (E_Mi * E_Ni)) + (Var_Ni / (E_Ni^2)))
#'   
#'   # Step 9: Compute Var(M_i / N_i) via delta-method approx 
#'   var_cs <- (Var_Mi / E_Ni^2) - 
#'     2 * ( E_Mi / (E_Ni^3) ) * Cov_Mi_Ni + 
#'     ( (E_Mi^2) / (E_Ni^4) ) * Var_Ni 
#'   
#'   # precompute a_k
#'   a_k <- colSums(p_jk * v_j)   # length K
#'   
#'   # variances of enrichment
#'   var_e <- (tau_display$Enrichment_std_error)^2
#'   
#'   # for each cell i compute V_{E M,i}
#'   # l_ik is I×K, a_k is length K, so (l_ik * a_k) is I×K
#'   # square it, multiply by var_e, then rowSum
#'   V_EM <- rowSums( (l_ik * matrix(a_k, nrow=I, ncol=K, byrow=TRUE))^2
#'                    * matrix(var_e, nrow=I, ncol=K, byrow=TRUE) )
#'   
#'   # finally add into your var_cs
#'   var_cs_total <- var_cs + V_EM / (E_Ni^2)
#'   
#'   # Return the single-cell scores and LDSC results table
#'   return(list(cs = cs,
#'               var_cs = var_cs, 
#'               ldsc_res_table = tau_display,
#'               cs_dat = list(J = J, K =K, I=I, v_j=v_j,
#'                             E_Mi=E_Mi, E_Ni=E_Ni, Cov_Mi_Ni=Cov_Mi_Ni, 
#'                             Var_Ni=Var_Ni, Var_Mi=Var_Mi)))
#' }
#' 
#' 
