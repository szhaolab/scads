#' Cell-Type-Level Test for Heterogeneity in Disease Relevance
#'
#' Tests whether the cells within a cell type differ from one another in disease
#' relevance, giving one p-value per cell type. Summarising a cell type by its
#' mean score discards this question, which is often the biologically interesting
#' part.
#'
#' @section The statistic:
#' The cell score is linear in the topic enrichments, \eqn{s_i = w_i' E} with
#' \eqn{\sum_k w_{ik} = 1}, where the weights depend only on the cell's topic
#' loadings and the topic annotation sizes. For cell type \eqn{g}, with
#' \eqn{S_g} the covariance of \eqn{w} across its cells, the within-type variance
#' of the score is a quadratic form,
#' \deqn{T_g = \mathrm{Var}_g(s) = \hat E' S_g \hat E.}
#'
#' @section The null:
#' The null of *no heterogeneity* is that all topics are equally enriched,
#' \eqn{E = c \mathbf{1}} for some unknown \eqn{c}. This is not the null of no
#' enrichment: if every topic carries the same enrichment, every cell scores
#' \eqn{c} and the within-type variance is zero whatever the loadings. Because
#' the weights sum to one, \eqn{S_g \mathbf{1} = 0} exactly, so writing
#' \eqn{\hat E = c\mathbf{1} + \varepsilon} removes every term in \eqn{c} and
#' \deqn{T_g = \varepsilon' S_g \varepsilon \sim \sum_j \lambda_j \chi^2_1,}
#' with \eqn{\lambda} the eigenvalues of \eqn{\Sigma^{1/2} S_g \Sigma^{1/2}}. The
#' p-value is exact by Imhof inversion; no permutation is involved. That \eqn{c}
#' drops out matters in practice, since every topic annotation is open chromatin
#' and therefore genuinely enriched.
#'
#' @section Choice of Sigma:
#' Use [ldsc_jackknife_cov()]. The annotation-correlation approximation is *not*
#' adequate here: the statistic is orthogonal to the constant direction, which is
#' exactly where that approximation errs, and on null simulations it makes the
#' test conservative at every threshold. The test is also sensitive to the level
#' of the off-diagonal correlation, so an exchangeable stand-in is not a
#' substitute for the full matrix.
#'
#' @section Interpretation:
#' \eqn{T_g} and its null eigenvalues both scale linearly with \eqn{S_g}, so the
#' p-value is invariant to the overall size of the within-group loading spread.
#' Power depends on whether the *direction* in which cells vary within a group
#' aligns with the direction in which topic enrichments differ. A cell type's
#' score range therefore does not predict whether it will be called
#' heterogeneous, and a range-based summary will legitimately disagree with these
#' calls.
#'
#' Interpret heterogeneity only for cell types that pass the association test:
#' heterogeneity within a cell type carrying no overall signal is heterogeneity
#' in noise. The test also conditions on the topic loadings as fixed, ignoring
#' topic-model uncertainty, as the cell score variance does.
#'
#' @param topic_res Topic-model results with `Pmat` (peaks x topics) and `Lmat`
#'   (cells x topics), as returned by [run_fastTopics()].
#' @param enrichment Length-K vector of topic enrichment estimates.
#' @param Sigma `K x K` covariance of `enrichment`; the `Sigma` element of
#'   [ldsc_jackknife_cov()].
#' @param groups Character or factor vector of cell-type labels, one per cell,
#'   in the row order of `topic_res$Lmat`. `NA` labels are dropped.
#' @param min_cells Groups smaller than this are skipped. Default 20.
#' @param p_adjust_method Passed to [stats::p.adjust()]. Default `"BH"`.
#' @return A data frame with one row per cell type: `cell_type`, `n`, `stat`
#'   (\eqn{T_g}), `p`, `fdr`, sorted by p-value.
#' @seealso [ldsc_jackknife_cov()], [get_cs()]
#' @export
cell_type_heterogeneity <- function(topic_res, enrichment, Sigma, groups,
                                    min_cells = 20, p_adjust_method = "BH") {

  if (is.null(topic_res$Pmat) || is.null(topic_res$Lmat))
    stop("topic_res must contain 'Pmat' and 'Lmat'")
  l_ik <- topic_res$Lmat
  K <- ncol(l_ik)
  if (length(enrichment) != K) stop("enrichment must have length ", K)
  if (!is.matrix(Sigma) || any(dim(Sigma) != c(K, K)))
    stop("Sigma must be a ", K, " x ", K, " matrix")
  if (length(groups) != nrow(l_ik))
    stop("groups must have one entry per cell (", nrow(l_ik), ")")

  keep  <- !is.na(groups)
  l_ik  <- l_ik[keep, , drop = FALSE]
  groups <- as.character(groups)[keep]

  # weights sum to one within each cell, which is what makes S_g %*% 1 == 0
  a_k <- Matrix::colSums(topic_res$Pmat)
  W <- sweep(l_ik, 2, a_k, "*")
  W <- W / rowSums(W)

  gs <- names(which(table(groups) >= min_cells))
  if (!length(gs))
    stop("no group has at least ", min_cells, " cells")

  out <- do.call(rbind, lapply(gs, function(g) {
    Wg  <- W[groups == g, , drop = FALSE]
    S_g <- stats::cov(Wg)
    Tg  <- as.numeric(t(enrichment) %*% S_g %*% enrichment)
    p   <- .quadform_pvalue(Tg, S_g, Sigma)
    data.frame(cell_type = g, n = nrow(Wg), stat = Tg, p = p,
               stringsAsFactors = FALSE)
  }))
  out$fdr <- stats::p.adjust(out$p, method = p_adjust_method)
  out[order(out$p), ]
}

#' Upper-tail probability of a weighted sum of chi-square(1) variables
#'
#' Imhof inversion when CompQuadForm is available, with a Satterthwaite
#' moment-match fallback. The eigenvalues here are typically of order 1e-2, at
#' which scale `imhof` and `davies` return out-of-range values with `ifault` set;
#' rescaling by their sum fixes this and leaves the p-value unchanged.
#' @noRd
.quadform_pvalue <- function(Tg, S_g, Sigma) {
  if (!is.finite(Tg) || Tg <= 0) return(NA_real_)
  R   <- chol(Sigma + diag(1e-12, ncol(Sigma)))
  lam <- eigen(R %*% S_g %*% t(R), symmetric = TRUE, only.values = TRUE)$values
  lam <- lam[lam > max(lam) * 1e-10]
  if (!length(lam)) return(NA_real_)
  sc <- sum(lam)
  if (requireNamespace("CompQuadForm", quietly = TRUE)) {
    p <- tryCatch(CompQuadForm::imhof(Tg / sc, lambda = lam / sc)$Qq,
                  error = function(e) NA_real_)
    if (is.finite(p)) return(min(max(p, 0), 1))
  }
  m1 <- sum(lam); m2 <- 2 * sum(lam^2)
  stats::pchisq(Tg / (m2 / (2 * m1)), df = 2 * m1^2 / m2, lower.tail = FALSE)
}
