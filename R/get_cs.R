#' Compute Single-Cell Level Scores (c_i)
#'
#' The cell score is a weighted average of the per-topic enrichment estimates,
#' `cs_i = sum_k w_ik ehat_k` with `w_ik = l_ik s_k / sum_k l_ik s_k`, so the
#' weights sum to one and the score equals 1 when no topic is enriched. The
#' z-score tests `cs_i = 1` (no enrichment in the cell's open regions relative to
#' genome-wide).
#'
#' Two filtering rules are applied to the S-LDSC estimates before they enter the
#' score. A topic whose annotation covers < `min_prop_snps` of the genome, or
#' whose estimated enrichment is negative, is floored to `e_k = 1` with
#' `w_k = 0`. Such topics contribute nothing to the numerator and nothing to the
#' variance, but they remain in the weight denominator, so the effect is to
#' dilute the score toward 1 rather than to renormalise over the retained topics.
#' This is deliberate and conservative.
#'
#' @param topic_res Results from run_fastTopics_res, with Pmat (J×K) and Lmat (I×K)
#' @param ldsc_res_dir Directory containing k*_output/results/Trait.results
#' @param trait Trait name used in LDSC results files.
#' @param nTopics Number of topics (K).
#' @param Sigma Optional K x K covariance matrix of the enrichment estimates, as
#'   returned by [ldsc_jackknife_cov()]. When supplied it is used directly for
#'   the cell score variance. When `NULL` (default) the variance falls back to
#'   the annotation-correlation approximation
#'   `Cov(e_k, e_k') ~ w_k w_k' Cor(A_k, A_k')`, which captures only peak overlap
#'   and ignores the correlation induced by the shared GWAS, LD reference and
#'   baseline-LD covariates.
#' @param alternative Either `"two.sided"` (default, preserving previous
#'   behaviour) or `"greater"` for a one-sided test of enrichment. Under
#'   `"two.sided"`, cells that are significantly *depleted* (z < 0) also pass
#'   FDR; if the scores are presented as disease relevance, `"greater"` is
#'   usually the intended test.
#' @param min_prop_snps Annotations covering less than this fraction of the
#'   genome are floored (see Details). Default 0.005.
#' @return A list with
#'   - cs: cell scores (vector of length I)
#'   - z_cell: z-score accounting for variance
#'   - p_cell: p-value for assessing score reliability, length I, `NA` where the
#'     variance is zero or non-finite
#'   - var_cell: cell score variance used in the denominator of z
#'   - Sigma_used: the K x K covariance actually used
#'   - ldsc_res_table: the raw tau_display table from S-LDSC
#'   - ash_res: results from running ASH on S-LDSC enrichment estimates 
#'   - cs_dat: intermediate estimates for calculating cell score 
#'   
#' @seealso [ldsc_jackknife_cov()]
#' @export
get_cs <- function(topic_res, ldsc_res_dir, trait, nTopics,
                   Sigma = NULL,
                   alternative = c("two.sided", "greater"),
                   min_prop_snps = 0.005) {

  alternative <- match.arg(alternative)

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
  # Flooring rule: negative enrichment, or an annotation smaller than
  # min_prop_snps of the genome, contributes e_k = 1 and w_k = 0. See Details.
  drop_k <- tau_display$Enrichment < 0 | tau_display$`Prop._SNPs` < min_prop_snps
  e_Ck <- ifelse(drop_k, 1, tau_display$Enrichment)
  se_e <- ifelse(drop_k, 0, tau_display$Enrichment_std_error)
  e_Ck_ash <- ifelse((tau_display_ash$result$PosteriorMean < 0 |
                        tau_display$`Prop._SNPs` < min_prop_snps),
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
  clean_names      <- normalize_peak_names(rownames(p_jk))
  valid            <- grepl("^[^:]+:[0-9]+-[0-9]+(:[+\\-\\*])?$", clean_names)
  p_jk             <- p_jk[valid, , drop = FALSE]
  rownames(p_jk)   <- clean_names[valid]
  gr_peaks <- GRanges(rownames(p_jk))
  med_peak_width <- median(GenomicRanges::width(gr_peaks))
  cat("\nMedian peak width: ", med_peak_width, "\n")
  
  num_total_peaks <- 3e9 / med_peak_width
  p_jk <- rbind(p_jk, matrix(0, nrow = num_total_peaks - nrow(p_jk), ncol = nTopics))
  corm <- cor(p_jk, p_jk)

  # Var(sum_k a_k e_k) = sum_k sum_k' a_k a_k' Cov(e_k, e_k') -- the FULL double
  # sum over ordered pairs. Writing it one-sided over k' <= k requires a factor
  # of 2 on the off-diagonal terms; omitting that factor understates the
  # variance whenever the annotations are positively correlated. Check: two
  # perfectly correlated topics with equal weights give Var = sigma^2, which the
  # full form recovers and the one-sided form returns as 0.75 * sigma^2.
  if (!is.null(Sigma)) {
    if (!is.matrix(Sigma) || any(dim(Sigma) != c(nTopics, nTopics))) {
      stop("Sigma must be a ", nTopics, " x ", nTopics, " matrix")
    }
    v_ik_diag <- rowSums((l_ik %*% Sigma) * l_ik)
  } else {
    v_ik_diag <- rowSums((v_ik %*% corm) * v_ik)
  }

  z_cell <- as.vector((l_ik %*% (e_Ck - 1)) / sqrt(v_ik_diag))

  # Keep length I: a cell whose topics were all floored has zero variance and
  # gets NA, rather than being dropped (which would shift every subsequent
  # p-value onto the wrong cell).
  p_cell <- rep(NA_real_, length(z_cell))
  ok <- is.finite(z_cell)
  p_cell[ok] <- if (alternative == "two.sided") {
    2 * pnorm(-abs(z_cell[ok]))
  } else {
    pnorm(z_cell[ok], lower.tail = FALSE)
  }
  
  return(list(
    cs            = cs,
    z_cell        = z_cell,
    p_cell        = p_cell,
    var_cell      = v_ik_diag,
    Sigma_used    = if (is.null(Sigma)) outer(se_e, se_e) * corm else Sigma,
    ldsc_res_table = tau_display,
    ash_res        = tau_display_ash$result,
    cs_dat = list(M_i=M_i, N_i=N_i))
  )
}
