#' Compute Single-Cell Level Scores (c_i)
#'
#' @param topic_res Results from run_fastTopics_res, with Pmat (J×K) and Lmat (I×K)
#' @param ldsc_res_dir Directory containing k*_output/results/Trait.results
#' @param trait Trait name used in LDSC results files.
#' @param nTopics Number of topics (K).
#' @return A list with
#'   - cs: cell scores (vector of length I)
#'   - z_cell: z-score accounting for variance 
#'   - p_cell: p-value for assessing score reliability
#'   - ldsc_res_table: the raw tau_display table from S-LDSC
#'   - ash_res: results from running ASH on S-LDSC enrichment estimates 
#'   - cs_dat: intermediate estimates for calculating cell score 
#'   
#' @export
get_cs <- function(topic_res, ldsc_res_dir, trait, nTopics) {

  # --- Input validation ---
  if (is.null(topic_res$Pmat)) stop("topic_res must contain 'Pmat' (peaks x topics binary matrix)")
  if (is.null(topic_res$Lmat)) stop("topic_res must contain 'Lmat' (cells x topics loading matrix)")
  if (!dir.exists(ldsc_res_dir)) stop("ldsc_res_dir does not exist: ", ldsc_res_dir)
  if (!is.character(trait) || nchar(trait) == 0) stop("trait must be a non-empty string")

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
  
  cat("\nRunning ASH on enrichment estimates")
  tau_display_ash <- ash(tau_display$Enrichment, tau_display$Enrichment_std_error)
  
  # 3) per-topic peak counts (assuming one variant per peak)
  a_k <- colSums(p_jk)  # length K
  
  # instead filter by annotation size 
  e_Ck <- ifelse((tau_display$Enrichment < 0 | tau_display$Prop._SNPs < 0.005), 
                 1, tau_display$Enrichment)
  se_e <- ifelse((tau_display$Enrichment < 0 | tau_display$Prop._SNPs < 0.005), 
                 0, tau_display$Enrichment_std_error)
  e_Ck_ash <- ifelse((tau_display_ash$result$PosteriorMean < 0 | tau_display$Prop._SNPs < 0.005), 
                     1, tau_display_ash$result$PosteriorMean)
  
  # 4) single-cell score
  M_i <- as.numeric(l_ik %*% (a_k * e_Ck_ash))
  N_i <- as.numeric(l_ik %*% a_k)
  cs <- M_i / N_i
  
  # 5) calc variance and get z-score and p-value for the cell
  
  # weight each topic 
  l_ik <- sweep(l_ik, 2, a_k, FUN = "*")/sum(a_k)
  v_ik <- sweep(l_ik, 2, se_e, FUN = "*")
  
  cat("\nAccounting for correlation")
  #--- accounting for correlation of topic annotations in LDSC.---
  # get number of total peaks based on peak width 
  gr_peaks <- GRanges(rownames(p_jk))
  med_peak_width <- median(GenomicRanges::width(gr_peaks))
  cat("\nMedian peak width: ", med_peak_width, "\n")
  
  num_total_peaks <- 3e9 / med_peak_width
  p_jk <- rbind(p_jk, matrix(0, nrow = num_total_peaks - nrow(p_jk), ncol = nTopics))
  corm <- cor(p_jk, p_jk)
  corm.upper<- corm * upper.tri(corm, diag = TRUE)
  # z_cell <- l_ik %*% (e_Ck-1)/ sqrt(diag(v_ik %*% corm.upper %*% t(v_ik)))
    # 1. Calculate the intermediate product (Cells x Topics)
    intermediate_prod <- v_ik %*% corm.upper 
    # 2. Calculate the diagonal without forming the full N x N matrix
    # (v_ik * intermediate_prod) performs element-wise multiplication
    v_ik_diag <- rowSums(intermediate_prod * v_ik)
    # 3. Calculate your Z-score
    z_cell <- (l_ik %*% (e_Ck - 1)) / sqrt(v_ik_diag)
    
  p_cell <- 2 * pnorm(-abs(z_cell))
  # p_cell <- 1 - pnorm(z_cell) 
  p_cell <- p_cell[is.finite(p_cell) & !is.na(p_cell)]
  
  return(list(
    cs            = cs,
    z_cell        = z_cell,
    p_cell        = p_cell,
    ldsc_res_table = tau_display,
    ash_res        = tau_display_ash$result,
    cs_dat = list(M_i=M_i, N_i=N_i))
  )
}
