#' Compute Single-Cell Level Scores (c_i)
#'
#' The cell score is a weighted average of the raw S-LDSC enrichment estimates,
#' `cs_i = sum_k w_ik ehat_k` with `w_ik = l_ik s_k / sum_k l_ik s_k`, so the
#' weights sum to one and the score equals 1 when no topic is enriched. The
#' enrichments are used without shrinkage. The
#' z-score tests `cs_i = 1` (no enrichment in the cell's open regions relative to
#' genome-wide).
#'
#' Two filtering rules are applied to the S-LDSC estimates before they enter the
#' score. A topic whose estimated enrichment is negative, or whose annotation is
#' too small for S-LDSC's p-values to control type I error (below), is floored
#' to `e_k = 1` with `w_k = 0`. Such topics contribute nothing to the numerator
#' and nothing to the variance, but they remain in the weight denominator, so
#' the effect is to dilute the score toward 1 rather than to renormalise over
#' the retained topics. This is deliberate and conservative.
#'
#' @section Small annotations:
#' For small annotations S-LDSC's block-jackknife p-values do not control type I
#' error: the tail of very small p-values is inflated (Tashman et al. 2021).
#' Following that paper, an annotation is judged by the percentage of the
#' genome's 0.01 cM blocks it overlaps, which reflects both its size and how
#' clustered it is (see [annot_block_coverage()]). Topics overlapping fewer than
#' `min_pct_blocks` percent (default 1.7, the paper's threshold for the one-sided
#' and enrichment tests) are reported in a warning and, with
#' `small_annot = "floor"` (default), floored. With `small_annot = "warn"` they
#' are kept, and the warning is the only safeguard.
#'
#' Block coverage needs the reference SNPs' genetic map positions, read from
#' `baseline_prefix`, or a precomputed `block_coverage`. Without either, the
#' previous rule is used instead, flooring topics with S-LDSC `Prop._SNPs` below
#' 0.5%, and a message says so.
#'
#' The reliability of the cell ranking is summarised by the resolvable-error rank
#' FDR ([rank_fdr()]), computed among the significant cells (one-sided BH <
#' `fdr_alpha`) from the enrichments and `Sigma`.
#'
#' @param topic_res Results from run_fastTopics_res, with Pmat (J×K) and Lmat (I×K)
#' @param ldsc_res_dir Directory containing k*_output/results/Trait.results
#' @param trait Trait name used in LDSC results files.
#' @param nTopics Number of topics (K).
#' @param baseline_prefix,frq_prefix Prefixes of the baseline annotation and
#'   PLINK frequency files, matching the `--ref-ld-chr` and `--frqfile-chr` of
#'   the S-LDSC runs. These are used to build the covariance of the enrichment
#'   estimates with [ldsc_jackknife_cov()], which is the default.
#' @param Sigma Optional K x K covariance matrix of the enrichment estimates,
#'   supplied directly instead of being computed.
#' @param use_approximation Use the annotation-correlation approximation
#'   `Cov(e_k, e_k') ~ w_k w_k' Cor(A_k, A_k')` instead of the jackknife
#'   covariance. Default `FALSE`; the jackknife is preferred, and one of
#'   `Sigma`, the reference prefixes, or this flag must be given.
#' @param alternative Either `"greater"` (default), a one-sided test of
#'   enrichment (cell score > 1), or `"two.sided"`, under which cells that are
#'   significantly *depleted* (z < 0) also pass FDR. Before 0.9.0 the default
#'   was `"two.sided"`.
#' @param min_pct_blocks Minimum percentage of the genome's 0.01 cM blocks a
#'   topic annotation must overlap (see Small annotations). Default 1.7.
#' @param block_coverage Optional precomputed coverage: the data frame returned
#'   by [annot_block_coverage()], or a numeric vector of `pct_blocks` of length
#'   K. Computed from `baseline_prefix` when `NULL`.
#' @param small_annot What to do with topics below `min_pct_blocks`:
#'   `"floor"` (default) floors them as described in Details; `"warn"` keeps
#'   them. Both warn.
#' @param min_prop_snps Deprecated. If given, the previous rule (floor topics
#'   with S-LDSC `Prop._SNPs` below this fraction) is used instead of block
#'   coverage.
#' @param fdr_alpha BH level defining the significant cells ranked by the rank
#'   FDR. Default 0.05.
#' @param fdr_draws Number of enrichment draws for the rank FDR. Default 1000.
#'   Set to 0 to skip it.
#' @param fdr_seed Random seed for the rank FDR draws. Default 1.
#' @param chrs Chromosomes the S-LDSC runs used, for the jackknife covariance and
#'   the block coverage. Default `1:22`.
#' @return A list with
#'   - cs: cell scores (vector of length I)
#'   - z_cell: z-score accounting for variance
#'   - p_cell: p-value for assessing score reliability, length I, `NA` where the
#'     variance is zero or non-finite
#'   - var_cell: cell score variance used in the denominator of z
#'   - Sigma_used: the K x K covariance actually used
#'   - ldsc_res_table: the raw tau_display table from S-LDSC
#'   - rank_fdr: resolvable-error rank FDR of the significant cells, as returned by
#'     [rank_fdr()] (`curve`, `mean_fdr`, `fdr_top10`, `n_sig`, `sig`); `NULL`
#'     if `fdr_draws = 0`
#'   - cs_dat: intermediate estimates for calculating cell score, including the
#'     normalised cell topic weights `W` and the floored-topic indicator `drop`
#'   - topic_qc: one row per topic with the S-LDSC enrichment, `Prop._SNPs`,
#'     `pct_blocks`, `small` (below the size threshold), `negative` (negative
#'     enrichment) and `floored`
#'
#' @references Tashman KC, Cui R, O'Connor LJ, Neale BM, Finucane HK (2021).
#'   Significance testing for small annotations in stratified LD-Score
#'   regression. medRxiv 2021.03.13.21249938.
#' @seealso [ldsc_jackknife_cov()], [annot_block_coverage()], [rank_fdr()]
#' @export
get_cs <- function(topic_res, ldsc_res_dir, trait, nTopics,
                   baseline_prefix = NULL, frq_prefix = NULL,
                   Sigma = NULL,
                   use_approximation = FALSE,
                   alternative = c("greater", "two.sided"),
                   min_pct_blocks = 1.7,
                   block_coverage = NULL,
                   small_annot = c("floor", "warn"),
                   min_prop_snps = NULL,
                   fdr_alpha = 0.05, fdr_draws = 1000, fdr_seed = 1,
                   chrs = 1:22) {

  alternative <- match.arg(alternative)
  small_annot <- match.arg(small_annot)


  # --- Input validation ---
  if (is.null(topic_res$Pmat)) stop("topic_res must contain 'Pmat' (peaks x topics binary matrix)")
  if (is.null(topic_res$Lmat)) stop("topic_res must contain 'Lmat' (cells x topics loading matrix)")
  if (!dir.exists(ldsc_res_dir)) stop("ldsc_res_dir does not exist: ", ldsc_res_dir)
  if (!is.character(trait) || nchar(trait) == 0) stop("trait must be a non-empty string")

  # Covariance of the enrichment estimates. The block jackknife is the default;
  # the annotation-correlation approximation must be asked for explicitly.
  if (is.null(Sigma) && !isTRUE(use_approximation)) {
    if (is.null(baseline_prefix) || is.null(frq_prefix))
      stop("get_cs() needs the covariance of the enrichment estimates.\n",
           "  Supply baseline_prefix and frq_prefix (matching this run's ",
           "--ref-ld-chr and --frqfile-chr) to build it with ",
           "ldsc_jackknife_cov(),\n  or pass a precomputed Sigma, ",
           "or set use_approximation = TRUE to use ",
           "se_k se_m Cor(A_k, A_m) instead.", call. = FALSE)
    jk <- ldsc_jackknife_cov(ldsc_res_dir = ldsc_res_dir, trait = trait,
                             nTopics = nTopics,
                             baseline_prefix = baseline_prefix,
                             frq_prefix = frq_prefix, chrs = chrs, verbose = FALSE)
    if (!jk$pass)
      warning("jackknife covariance failed its acceptance check; ",
              "verify baseline_prefix and frq_prefix match this run's S-LDSC call")
    Sigma <- jk$Sigma
  }

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
  
  # Annotation size: share of the genome's 0.01 cM blocks overlapped
  # (Tashman et al. 2021), or the deprecated Prop._SNPs rule. See Small annotations.
  pct_blocks <- rep(NA_real_, nTopics)
  if (!is.null(min_prop_snps)) {
    warning("min_prop_snps is deprecated; topic size is now judged by 0.01 cM ",
            "block coverage (min_pct_blocks). Using Prop._SNPs < ", min_prop_snps,
            " as requested.", call. = FALSE)
    small_k <- tau_display$`Prop._SNPs` < min_prop_snps
    size_rule <- paste0("Prop._SNPs < ", min_prop_snps)
  } else {
    if (is.null(block_coverage) && !is.null(baseline_prefix))
      block_coverage <- tryCatch(
        annot_block_coverage(ldsc_res_dir = ldsc_res_dir, trait = trait,
                             nTopics = nTopics, baseline_prefix = baseline_prefix,
                             chrs = chrs),
        error = function(e) {
          message("Could not compute 0.01 cM block coverage: ", conditionMessage(e))
          NULL })
    if (!is.null(block_coverage)) {
      pct_blocks <- if (is.data.frame(block_coverage)) block_coverage$pct_blocks else
        as.numeric(block_coverage)
      if (length(pct_blocks) != nTopics)
        stop("block_coverage must have one value per topic (", nTopics, ")")
      small_k <- pct_blocks < min_pct_blocks
      size_rule <- paste0("< ", min_pct_blocks, "% of 0.01 cM blocks")
    } else {
      message("Block coverage unavailable (no baseline_prefix or block_coverage); ",
              "flooring topics with Prop._SNPs < 0.005 instead.")
      small_k <- tau_display$`Prop._SNPs` < 0.005
      size_rule <- "Prop._SNPs < 0.005"
    }
  }
  if (any(small_k))
    warning(sum(small_k), " topic(s) below the annotation-size threshold (",
            size_rule, "): ",
            paste0("k", which(small_k),
                   if (all(is.na(pct_blocks))) "" else
                     paste0(" (", round(pct_blocks[small_k], 2), "%)"),
                   collapse = ", "),
            ". S-LDSC p-values for annotations this small may not control type I ",
            "error (Tashman et al. 2021); ",
            if (small_annot == "floor") "these topics are floored to e = 1." else
              "they are kept because small_annot = \"warn\".", call. = FALSE)

  # Flooring rule: negative enrichment, or a small annotation when
  # small_annot = "floor", contributes e_k = 1 and w_k = 0. See Details.
  floor_small <- small_k & small_annot == "floor"
  drop_k <- tau_display$Enrichment < 0 | floor_small
  e_Ck <- ifelse(drop_k, 1, tau_display$Enrichment)
  se_e <- ifelse(drop_k, 0, tau_display$Enrichment_std_error)
  topic_qc <- data.frame(topic = paste0("k", seq_len(nTopics)),
                         Enrichment = tau_display$Enrichment,
                         Prop._SNPs = tau_display$`Prop._SNPs`,
                         pct_blocks = pct_blocks, small = small_k,
                         negative = tau_display$Enrichment < 0,
                         floored = drop_k, check.names = FALSE)
  
  # 4) single-cell score, from the raw (floored) enrichments
  M_i <- as.numeric(l_ik %*% (a_k * e_Ck))
  N_i <- as.numeric(l_ik %*% a_k)
  cs <- M_i / N_i
  
  # 5) calc variance and get z-score and p-value for the cell
  
  # weight each topic 
  l_ik <- sweep(l_ik, 2, a_k, FUN = "*")/sum(a_k)
  W <- l_ik / rowSums(l_ik)   # normalised cell topic weights, for the rank FDR
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
    # floored topics contribute nothing to the variance (see Details)
    Sigma_f <- Sigma
    Sigma_f[drop_k, ] <- 0; Sigma_f[, drop_k] <- 0
    v_ik_diag <- rowSums((l_ik %*% Sigma_f) * l_ik)
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
  
  # 6) resolvable-error rank FDR among the significant cells
  Sigma_used <- if (is.null(Sigma)) outer(tau_display$Enrichment_std_error,
                                          tau_display$Enrichment_std_error) * corm else Sigma
  fdr_res <- NULL
  if (fdr_draws > 0) {
    cat("\nResolvable-error rank FDR (", fdr_draws, " draws)", sep = "")
    fdr_res <- rank_fdr(W, tau_display$Enrichment, Sigma_used, small = floor_small,
                        alpha = fdr_alpha, n_draws = fdr_draws, seed = fdr_seed)
    if (is.na(fdr_res$mean_fdr)) {
      cat("\nFewer than two significant cells: rank FDR not defined\n")
    } else {
      cat(sprintf("\n%d significant cells; estimated FDR top 10%%: %.2f, mean over top 1-50%%: %.2f\n",
                  fdr_res$n_sig, fdr_res$fdr_top10, fdr_res$mean_fdr))
    }
  }

  return(list(
    cs            = cs,
    z_cell        = z_cell,
    p_cell        = p_cell,
    var_cell      = v_ik_diag,
    Sigma_used    = Sigma_used,
    ldsc_res_table = tau_display,
    rank_fdr       = fdr_res,
    cs_dat = list(M_i=M_i, N_i=N_i, W=W, drop=drop_k),
    topic_qc       = topic_qc)
  )
}
