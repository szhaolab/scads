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
  nTopics <- ncol(p_jk)
  
  # 2) read in LDSC results per topic
  tau_display <- NULL
  for(k in seq_len(nTopics)){
    res_f <- file.path(ldsc_res_dir,
                       paste0("k", k, "_output"),
                       "results",
                       paste0(trait, ".results"))
    if(!file.exists(res_f)){
      warning("Missing ", res_f)
      next
    }
    res <- read.table(res_f, header=TRUE, sep="\t", check.names=FALSE)
    if(nrow(res)<1){
      warning("Empty ", res_f)
      next
    }
    res$Category <- paste0("k", k)
    if(is.null(tau_display)){
      tau_display <- res[1,]
    } else {
      tau_display <- rbind(tau_display, res[1,])
    }
  }
  
  # 3) per-topic peak counts (assuming one variant per peak)
  a_k <- colSums(p_jk)  # length K
  
  # 4) compute M_i and N_i for each cell
  #    M_i = Σ_k  l_ik[i,k] * (a_k * e_Ck)[k]
  #    N_i = Σ_k  l_ik[i,k] * a_k[k]
  
  # first filter out negative enrichment 
  e_Ck   <- ifelse(tau_display$Enrichment < 0, 0, tau_display$Enrichment)
  se_e   <- ifelse(tau_display$Enrichment < 0, 0, tau_display$Enrichment_std_error)
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


# library(data.table)
# library(GenomicRanges)

#' Compute Single‑Cell Topic Enrichment Scores (c_i)
#'
#' This function computes per‑cell enrichment scores \(c_i\) by combining
#' fastTopics outputs (peak‑by‑topic probabilities and cell‑by‑topic loadings)
#' with per‑topic LDSC enrichment (or heritability τ) estimates and optional
#' GWAS variant counts per peak. Both first‑order (only uncertainty in LDSC
#' enrichment) and full second‑order (ratio variance via delta‑method)
#' approximations are supported.
#'
#' @param topic_res A list returned by \code{run_fastTopics}, containing
#'   \code{Pmat} (J × K peak × topic probability matrix) and
#'   \code{Lmat} (I × K cell × topic loadings matrix).
#' @param ldsc_res_dir Path to the directory where, for each topic k,
#'   the file \code{k<k>_output/results/<trait>.results} holds
#'   LDSC output with columns \code{Enrichment} and
#'   \code{Enrichment_std_error}.
#' @param trait Character string; the trait name used to locate
#'   \code{<trait>.results} within each topic’s LDSC results folder.
#' @param sumstats_dir Optional character: path to a folder of
#'   \code{<trait>_sumstats.txt.gz} summary‑stats files. Required if
#'   \code{v_option="gwas"}.
#' @param v_option Character, one of \code{"one"} (assume exactly one
#'   variant per peak) or \code{"gwas"} (count variants per peak
#'   from GWAS summary‑stats via \code{sumstats_dir}).
#' @param weight_by Character, one of \code{"enrichment"} (use LDSC
#'   \code{Enrichment} estimates) or \code{"tau"} (use raw
#'   \code{Prop._h2} τ estimates) as peak weights.
#' @param cs_approx Character, one of \code{"first"} (delta‑method
#'   propagating only LDSC enrichment SE) or \code{"second"} (full
#'   ratio variance using E[M], E[N], Var(M), Var(N), Cov(M,N)).
#'
#' @return A list with
#'   \item{cs}{Numeric vector of length I, the per‑cell scores \(c_i\).}
#'   \item{var_cs}{Numeric vector of the same length, the variance of \(c_i\).}
#'   \item{ldsc_res_table}{Data.frame of per‑topic LDSC results (enrichment & SE).}
#'   \item{cs_dat}{A list with
#'     \item{M_i}{Vector of the numerator \(M_i\) used in the score.}
#'     \item{N_i}{Vector of the denominator \(N_i\).}
#'     \item{Var_Mi}{The variance of \(M_i\) under the chosen approximation.}
#'   }
#'
#' @export
get_cs2 <- function(topic_res,
                    ldsc_res_dir,
                    trait,
                    sumstats_dir     = NULL,
                    v_option         = "one", # c("one","gwas"),
                    weight_by        = "enrichment", # c("enrichment","tau"),
                    cs_approx        = "first" # c("first","second")
){
  ## ——————————————————————————————————————
  ## 1. Unpack inputs and set up
  v_option  <- match.arg(v_option)
  weight_by <- match.arg(weight_by)
  cs_approx <- match.arg(cs_approx)
  
  p_jk <- topic_res$Pmat    # J × K
  l_ik <- topic_res$Lmat    # I × K
  I    <- nrow(l_ik)
  J    <- nrow(p_jk)
  K    <- ncol(p_jk)
  
  ## ——————————————————————————————————————
  ## 2. Read LDSC results to get tau_Ck, e_Ck, se_e_Ck
  tau_display <- NULL
  # tau_Ck   <- numeric(K)
  # e_Ck     <- numeric(K)
  # se_e_Ck  <- numeric(K)
  for(k in seq_len(K)){
    res_f <- file.path(ldsc_res_dir,
                       paste0("k", k, "_output"),
                       "results",
                       paste0(trait, ".results"))
    if(!file.exists(res_f)){
      warning("Missing ", res_f); tau_Ck[k] <- NA; next
    }
    d <- fread(res_f)
    if(nrow(d)<1){
      warning("Empty ", res_f); tau_Ck[k] <- NA; next
    }
    d$Category <- paste0("k", k)
    tau_display <- rbind(tau_display, d[1,])
    # tau_Ck[k]  <- d$Prop._h2[1]
    # e_Ck[k]    <- d$Enrichment[1]
    # se_e_Ck[k] <- d$Enrichment_std_error[1]
  }
  # filter out negative enrichment 
  e_Ck  <- ifelse(tau_display$Enrichment < 0, 0, tau_display$Enrichment)
  se_e  <- ifelse(tau_display$Enrichment < 0, 0, tau_display$Enrichment_std_error)
  tau_Ck  <- ifelse(tau_display$Enrichment < 0, 0, tau_display$Prop._h2)
  w_k <- if(weight_by=="enrichment") e_Ck else tau_Ck
  
  ## ——————————————————————————————————————
  ## 3. Compute v_j (variants per peak)
  if(v_option=="one"){
    v_j <- rep(1L, J)
  } else {
    # parse peaks from rownames of Pmat
    peaks    <- rownames(p_jk)
    peaks_dt <- as.data.table(tstrsplit(peaks, "_",
                                        keep = 1:3,
                                        col.names = c("seqnames","start","end")))
    peaks_dt[, `:=`(
      start = as.integer(start),
      end   = as.integer(end),
      peak  = peaks
    )]
    peaks_gr <- GRanges(peaks_dt$seqnames,
                        IRanges(peaks_dt$start, peaks_dt$end),
                        peak = peaks_dt$peak)
    
    # read GWAS sumstats positions
    sumstats_gz <- file.path(sumstats_dir,
                             sprintf("%s_sumstats.txt.gz", trait))
    ss <- fread(cmd = paste("zcat", sumstats_gz),
                select = c("chr","pos"))
    ss[, seqnames := paste0("chr", chr)]
    snps_gr <- GRanges(ss$seqnames,
                       IRanges(ss$pos, ss$pos))
    
    # overlap counts
    v_j <- countOverlaps(peaks_gr, snps_gr)
  }
  
  
  ## ——————————————————————————————————————
  ## 6. Return requested result
  if(cs_approx == "first"){
    ## 4. First‐order approximation (uncertainty only from se_e_Ck)
    # 4a) M_i & N_i
    a_k <- colSums(p_jk * v_j)            # sum_j p_jk * v_j  → length K
    M_i <- as.numeric(l_ik %*% (a_k * w_k))
    N_i   <- as.numeric(l_ik %*% a_k)
    cs_1  <- M_i / N_i
    cat("First-order CS approx: ", summary(cs_1))
    
    # 4b) var(M_i) from var(e_Ck) and propagate to var(cs)
    var_e    <- se_e_Ck^2
    V_Mi   <- rowSums( (l_ik * a_k)^2 *
                           matrix(var_e, nrow=I, ncol=K, byrow=TRUE) )
    var_cs_1 <- V_Mi / (N_i^2)
    cat("First-order CS variance approx: ", summary(var_cs_1))
    
    return(list(
      cs            = cs_1,
      var_cs        = var_cs_1,
      ldsc_res_table = tau_display,
      method       = "first-order",
      cs_dat = list(M_i=M_i, N_i=N_i, Var_Mi = V_Mi)))
    
  } else {
    
    # ## 5. Second‐order (full delta‐method) approximation
    # 5a) First-moments (M_i & N_i)
    a_k <- colSums(p_jk * v_j)    #  ∑_j v_j p_jk
    E_Mi <- as.numeric(l_ik %*% (a_k * w_k)) #  ∑_j v_j p_jk w_k
    E_Ni   <- as.numeric(l_ik %*% a_k)
    
    # 5b) Var(N_i)
    Var_Ni <- numeric(I)  # result vector of length I
    for (i in 1:I) {
      var_sum <- 0
      for (j in 1:J) {
        # sum_k L[i,k] * P[j,k]
        pij <- sum(l_ik[i, ] * p_jk[j, ])
        var_sum <- var_sum + v_j^2 * pij * (1 - pij)
      }
      Var_Ni[i] <- var_sum
    }

    
    # Now n_ij as a sparse Matrix:
    n_ij <- l_ik %*% t(P)             # I × J, computed via sparse×dense
    m_ij <- l_ik %*% t(P %*% Diagonal(x = w_k))


    # 5c) Var(N_i) and Cov(M_i,N_i)
    Var_Ni    <- rowSums( (v_j^2) * (n_ij * (1 - n_ij)) )
    Cov_Mi_Ni <- rowSums( (v_j^2) * (m_ij  * (1 - n_ij)) )

    # Calculate cell score
    cs_2 <- (E_Mi / E_Ni) *
      (1
       - Cov_Mi_Ni / (E_Mi * E_Ni)
       + Var_Ni      / (E_Ni^2)
      )
    cat("Second-order CS approx: ", summary(cs_2))

    # 5d) Var(M_i)
    #  — full Var(M_i) including enrichment‐uncertainty  —
    mu2_tau  <- e_Ck^2       # length K
    var2_tau <- se_e_Ck^2    # length K
    # Build the two J×K matrices:
    # S1[j,k] = mu2_tau[k] * p_jk[j,k] * (1 - p_jk[j,k])
    # S2[j,k] = var2_tau[k]  * p_jk[j,k]
    S1 <- p_jk * (mu2_tau * (1 - p_jk))   # elementwise
    S2 <- p_jk * var2_tau
    # Sum them: each peak j contributes
    #    sum_k l_ik[i,k] * (S1[j,k] + S2[j,k])
    # so we form the I×J matrix of those terms:
    Mvar_terms <- S1 + S2                # still J×K
    Mvar_ij    <- l_ik %*% t(Mvar_terms)  # I×J
    # Finally weight by v_j^2 and sum over peaks j:
    Var_Mi_full <- rowSums( (v_j^2) * Mvar_ij )

    # 5e) delta‐method for Var(M_i/N_i)
    var_cs_2 <- (Var_Mi_full / E_Ni^2) -
      2*(E_Mi/(E_Ni^3))*Cov_Mi_Ni +
      (E_Mi^2/(E_Ni^4))*Var_Ni
    cat("Second-order CS var approx: ", summary(var_cs_2))
    
    # Another way?
    ## — fast O(JK + IK^2) computation of Var(N_i) and Cov(M_i,N_i) —
    # A) First-moments
    a_k <- colSums(p_jk * v_j)    #  ∑_j v_j p_jk
    # b_k <- a_k * w_k     #  ∑_j v_j p_jk w_k
    b_k <- colSums(sweep(p_jk, 2, e_Ck, "*") * v_j)

    # B) Second-moments for Var(N_i) and Cov(M_i,N_i)
    c1    <- colSums(v_j^2 * p_jk)               # ∑_j v_j² p_jk
    c1_CN <- colSums(v_j^2 * p_jk * w_k)         # ∑_j v_j² p_jk w_k=

    D2    <- crossprod(p_jk * v_j^2, p_jk)        # K×K: ∑_j v_j² p_jk p_jk′
    D2_CN <- crossprod((p_jk * w_k) * v_j^2, p_jk)

    # C) Full Var(M_i) including τ‐uncertainty
    mu2_tau  <- e_Ck^2
    var2_tau <- se_e_Ck^2
    # c_k_M    <- colSums(v_j^2 * (mu2_tau * p_jk * (1 - p_jk) +
    #                                var2_tau * p_jk))
   

    # D) collapse to per-cell (I×1) in O(IK²)
    E_Ni       <- as.numeric(l_ik %*% a_k)
    E_Mi       <- as.numeric(l_ik %*% b_k)

    Term1_N    <- as.numeric(l_ik %*% c1)
    Term2_N    <- rowSums((l_ik %*% D2)    * l_ik)
    Var_Ni     <- Term1_N - Term2_N       # ≥0

    Term1_CN   <- as.numeric(l_ik %*% c1_CN)
    Term2_CN   <- rowSums((l_ik %*% D2_CN) * l_ik)
    Cov_Mi_Ni  <- Term1_CN - Term2_CN

    # Var_Mi_full<- as.numeric(l_ik %*% c_k_M)
    S1 <- sweep(p_jk * (1 - p_jk), 2, mu2_tau, "*")
    S2 <- sweep(p_jk, 2, var2_tau, "*")
    Mvar_terms <- S1 + S2
    Mvar_ij <- l_ik %*% t(Mvar_terms)
    Var_Mi_full <- rowSums(v_j^2 * Mvar_ij)

    # E) ratio‐moment and its delta‐method variance
    cs_2    <- (E_Mi / E_Ni) * (
      1
      - Cov_Mi_Ni / (E_Mi * E_Ni)
      + Var_Ni      / (E_Ni^2)
    )
    var_cs_2<- (Var_Mi_full / E_Ni^2) -
      2*(E_Mi/(E_Ni^3))*Cov_Mi_Ni +
      (E_Mi^2/(E_Ni^4))*Var_Ni
    
    
    return(list(
      cs            = cs_2,
      var_cs        = var_cs_2,
      ldsc_res_table = tau_display,
      method       = "second-order",
      cs_dat = list(E_Mi         = E_Mi,
                    E_Ni         = E_Ni,
                    Var_Mi       = Var_Mi_full,
                    Var_Ni       = Var_Ni,
                    Cov_Mi_Ni    = Cov_Mi_Ni)))
  }
  
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
