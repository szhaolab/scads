#' Covariance of S-LDSC Enrichment Estimates from Block-Jackknife Delete Values
#'
#' Reconstructs the `K x K` sampling covariance of the per-topic heritability
#' enrichment estimates from the block-jackknife delete values that S-LDSC writes
#' when run with `--print-delete-vals`.
#'
#' This replaces the approximation
#' `Cov(e_k, e_m) ~ se_k se_m Cor(A_k, A_m)`, which substitutes the correlation
#' between annotations for the correlation between their LD scores. That
#' approximation captures only peak overlap: estimates from separate S-LDSC runs
#' are also correlated because they share the GWAS summary statistics, the LD
#' reference panel and the baseline-LD covariates, so two topics with disjoint
#' peak sets still have correlated estimates.
#'
#' @section When it matters:
#' For the cell score the approximation is usually adequate, because the score's
#' weight vector has a large component along the constant direction where the two
#' agree closely. For [cell_type_heterogeneity()] it is not: that statistic is
#' constructed to be orthogonal to the constant direction, which is exactly where
#' the approximation error concentrates.
#'
#' @section Method:
#' Per LD block `b`, for topic `k`,
#' \deqn{\hat E_{k,(b)} = \frac{\sum_j \mathrm{overlap}[k,j]\,\tau_{j,(b)}}{\sum_j M_j\,\tau_{j,(b)}} \Big/ \mathrm{Prop.\_SNPs}[k]}
#' and Sigma is the ordinary block-jackknife covariance of those vectors. Two
#' details are required for this to reproduce S-LDSC's own output:
#' \itemize{
#'   \item **Overlap correction.** S-LDSC reports overlap-corrected per-category
#'     heritability (`_overlap_output`), weighting each category by its overlap
#'     with every other annotation. The uncorrected ratio
#'     `M_k tau_k / sum_j M_j tau_j` is a different quantity, and annotations
#'     overlap heavily.
#'   \item **MAF restriction.** When `--frqfile-chr` is supplied S-LDSC uses
#'     `M_5_50`, so the overlap must be computed over MAF 5-50% variants only.
#' }
#' With both applied the reconstruction matches the reported `Enrichment` and
#' `Enrichment_std_error` to four decimal places; `check` reports that agreement.
#'
#' @param ldsc_res_dir Directory holding `k*_output/` for each topic, as produced
#'   by [run_sldsc()].
#' @param trait Trait name used in the S-LDSC output filenames.
#' @param nTopics Number of topics (K).
#' @param baseline_prefix Prefix of the baseline annotation files, matching the
#'   `--ref-ld-chr` of the runs being reconstructed, e.g.
#'   `".../baselineLD_v2.2/baselineLD."`.
#' @param frq_prefix Prefix of the PLINK frequency files matching
#'   `--frqfile-chr`, e.g. `".../1000G.EUR.QC."`. Used for the MAF 5-50%
#'   restriction.
#' @param chrs Chromosomes to use. Defaults to 1:22.
#' @param tol_enrichment,tol_se Relative tolerances for the acceptance check.
#' @param verbose Print per-chromosome progress and the acceptance table.
#' @return A list with
#'   \item{Sigma}{the `K x K` covariance of the enrichment estimates}
#'   \item{check}{per-topic comparison against S-LDSC's reported values}
#'   \item{pass}{`TRUE` if every topic is within tolerance}
#'   \item{enrichment,se}{the reported point estimates and standard errors}
#'   \item{E_blocks}{`n_blocks x K` matrix of per-block enrichments}
#' @seealso [get_cs()], [cell_type_heterogeneity()]
#' @export
ldsc_jackknife_cov <- function(ldsc_res_dir, trait, nTopics,
                               baseline_prefix, frq_prefix,
                               chrs = 1:22,
                               tol_enrichment = 0.01, tol_se = 0.05,
                               verbose = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE))
    stop("data.table is required for ldsc_jackknife_cov()")
  fread <- data.table::fread

  kdir <- function(k) file.path(ldsc_res_dir, paste0("k", k, "_output"))
  res <- lapply(seq_len(nTopics), function(k) {
    f <- file.path(kdir(k), "results", paste0(trait, ".results"))
    if (!file.exists(f)) stop("missing S-LDSC results: ", f)
    utils::read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  })
  nA <- nrow(res[[1]])

  # overlap[k, ] and the focal annotation's own M, accumulated across chromosomes.
  # The baseline is read once per chromosome and reused across topics; never
  # materialise cbind(focal, baseline), since crossprod(f, [f|B]) = c(f.f, f'B).
  OV <- matrix(0, nTopics, nA); Mbase <- numeric(nA - 1); Mfoc <- numeric(nTopics)
  for (chr in chrs) {
    b <- fread(paste0(baseline_prefix, chr, ".annot.gz"), showProgress = FALSE)
    q <- fread(paste0(frq_prefix, chr, ".frq"), showProgress = FALSE)
    maf  <- q$MAF[match(b$SNP, q$SNP)]
    keep <- !is.na(maf) & maf >= 0.05 & maf <= 0.5
    bm <- as.matrix(b[, -(1:4)])[keep, , drop = FALSE]
    rm(b, q); gc(verbose = FALSE)
    if (ncol(bm) != nA - 1)
      stop("chr", chr, ": baseline has ", ncol(bm), " annotations but .results has ",
           nA, " categories; baseline_prefix does not match this run's --ref-ld-chr")
    Mbase <- Mbase + colSums(bm)
    for (k in seq_len(nTopics)) {
      f <- fread(file.path(kdir(k), "annotations", trait,
                           paste0(trait, ".", chr, ".annot.gz")), showProgress = FALSE)
      fv <- as.numeric(f[[1]])[keep]
      OV[k, ]  <- OV[k, ] + c(sum(fv * fv), as.vector(crossprod(fv, bm)))
      Mfoc[k]  <- Mfoc[k] + sum(fv)
      rm(f, fv)
    }
    if (verbose) { cat("chr", chr, "\r", sep = ""); utils::flush.console() }
    rm(bm); gc(verbose = FALSE)
  }
  if (verbose) cat("\n")
  Mtot <- Mbase[1]                                 # the all-ones 'base' column

  E_blocks <- NULL; chk <- NULL
  for (k in seq_len(nTopics)) {
    r  <- res[[k]]; p1 <- r$`Prop._SNPs`[1]
    ovp <- OV[k, ] / Mtot
    Mp  <- c(Mfoc[k], Mbase) / Mtot                # M differs per topic in its first element
    tau <- r$Coefficient
    E_full <- sum(ovp * tau) / sum(Mp * tau) / p1
    dv <- as.matrix(utils::read.table(
      file.path(kdir(k), "results", paste0(trait, ".part_delete"))))
    eb <- as.vector(dv %*% ovp) / as.vector(dv %*% Mp) / p1
    nb <- length(eb)
    se_jk <- sqrt((nb - 1) / nb * sum((eb - mean(eb))^2))
    chk <- rbind(chk, data.frame(
      topic = k,
      enrichment_reported = r$Enrichment[1], enrichment_jackknife = E_full,
      se_reported = r$Enrichment_std_error[1], se_jackknife = se_jk,
      pass = abs(E_full / r$Enrichment[1] - 1) < tol_enrichment &&
             abs(se_jk / r$Enrichment_std_error[1] - 1) < tol_se))
    E_blocks <- if (is.null(E_blocks)) matrix(eb, ncol = 1) else cbind(E_blocks, eb)
  }

  nb  <- nrow(E_blocks)
  Ec  <- sweep(E_blocks, 2, colMeans(E_blocks), "-")
  Sigma <- ((nb - 1) / nb) * crossprod(Ec)
  dimnames(Sigma) <- list(paste0("k", seq_len(nTopics)), paste0("k", seq_len(nTopics)))

  if (verbose) {
    cat("acceptance: ", sum(chk$pass), "/", nTopics, " topics reproduce S-LDSC\n", sep = "")
    R <- stats::cov2cor(Sigma)
    cat("mean off-diagonal Cor(e_k, e_m): ", round(mean(R[upper.tri(R)]), 4), "\n", sep = "")
  }
  if (!all(chk$pass))
    warning(sum(!chk$pass), " of ", nTopics, " topics fail the acceptance check; ",
            "verify baseline_prefix and frq_prefix match this run's S-LDSC call")

  list(Sigma = Sigma, check = chk, pass = all(chk$pass),
       enrichment = vapply(res, function(r) r$Enrichment[1], 0),
       se = vapply(res, function(r) r$Enrichment_std_error[1], 0),
       E_blocks = E_blocks)
}
