#' Resolvable-Error Rank FDR of the Cell Score
#'
#' Estimates how reliable the ranking of significant cells is. Calling the top
#' `k` significant cells by score, a call is false if the cell is below the
#' k-th score by more than the data can resolve: by more than one standard
#' deviation of the estimation noise of its score difference with the cells at
#' the cut-off. Near-ties within that noise are not counted as errors. The
#' estimated FDR at `k` is the average over the called cells of the probability
#' of such an error when the enrichments are drawn from `N(e_hat, Sigma)`.
#'
#' @section Definitions:
#' The cell score is `c_i = sum_k w_ik e_k`, with the same flooring rule as
#' [get_cs()]: a topic flagged in `small`, or whose (drawn) enrichment is
#' negative, contributes `e_k = 1`. Significant cells are those with one-sided
#' BH-adjusted p < `alpha` for `c_i > 1`; only they are ranked. The observed
#' order is by the score from `enrichment`.
#'
#' Margin: when the top `k` are called, called cell `i` gets
#' `m_i = sqrt(mean_j (w_i - w_j)' Sigma_f (w_i - w_j))`, the root-mean-square
#' estimation SD of its score difference with the cells `j` ranked within
#' 1% of the significant cells of position `k` in the observed order, where
#' `Sigma_f` is `Sigma` with the floored topics' rows and columns set to zero.
#' Cells with the same topic mix as the cells at the cut-off have a margin near
#' zero; cells from another part of the topic space have a larger one.
#'
#' Estimate: for each draw `e* ~ N(e_hat, Sigma)`, floored per draw, the drawn
#' scores are `s* = W e*` and `s*_(k)` is the k-th largest drawn score among
#' the significant cells. A called cell counts as false in the draw if
#' `s*_(k) - s*_i > m_i`. The FDR at `k` is the fraction of false calls,
#' averaged over draws.
#'
#' @section Interpretation:
#' The FDR is the expected share of called cells that are detectably worse than
#' the cut-off. It reflects the GWAS noise in the enrichment estimates; errors
#' in the fitted topic loadings are not modelled. The curve is usually highest
#' for the very top cells, where cells from different topics compete, and falls
#' toward zero as more cells are called, so it is best read at fixed fractions
#' (`fdr_top10`, or the `curve`) rather than as a cut-off rule. In simulations
#' (20 settings, about 480 runs) its average over runs agreed with
#' the actual share of resolvable errors within run-to-run noise, and that share
#' tracked how often a run put the wrong cell population on top. A single run
#' whose enrichment estimates confidently favour the wrong topic gets a low
#' FDR; such runs cannot be detected from the run itself.
#'
#' @param W `I x K` matrix of cell topic weights, rows summing to one
#'   (`w_ik` proportional to `l_ik` times the number of peaks in topic `k`).
#' @param enrichment Length-K vector of raw S-LDSC enrichment estimates.
#' @param Sigma `K x K` covariance of `enrichment`, e.g. the `Sigma` element of
#'   [ldsc_jackknife_cov()].
#' @param small Logical length-K vector, `TRUE` for topics floored for their
#'   annotation size (always contribute `e_k = 1`).
#' @param alpha BH level defining significant cells. Default 0.05.
#' @param pct Percentages of significant cells at which to evaluate the FDR.
#'   Default `1:50`.
#' @param n_draws Number of enrichment draws. Default 1000.
#' @param seed Optional random seed, for reproducibility.
#' @param chunk Draws processed at a time, to bound memory. Default 100.
#' @return A list with
#'   - `curve`: data frame with `pct` (percentage of significant cells called),
#'     `k` (number called) and `fdr` (estimated FDR)
#'   - `mean_fdr`: mean estimated FDR over the curve (top 1-50% by default)
#'   - `fdr_top10`: estimated FDR when calling the top 10% of significant cells
#'   - `n_sig`: number of significant cells
#'   - `sig`: logical length-I vector marking the significant cells
#'
#'   With fewer than two significant cells the FDR is undefined: `curve` is
#'   `NULL` and `mean_fdr`, `fdr_top10` are `NA`.
#' @seealso [get_cs()], which calls this and returns the result as `rank_fdr`.
#' @export
rank_fdr <- function(W, enrichment, Sigma, small = rep(FALSE, length(enrichment)),
                     alpha = 0.05, pct = 1:50, n_draws = 1000, seed = NULL,
                     chunk = 100) {
  K <- length(enrichment)
  if (ncol(W) != K) stop("W must have ", K, " columns, one per topic")
  if (!is.matrix(Sigma) || any(dim(Sigma) != c(K, K)))
    stop("Sigma must be a ", K, " x ", K, " matrix")
  if (length(small) != K) stop("small must have length ", K)

  fl <- function(e) ifelse(small | e < 0, 1, e)
  drop <- small | enrichment < 0
  Sf <- Sigma; Sf[drop, ] <- 0; Sf[, drop] <- 0

  score <- as.vector(W %*% fl(enrichment))
  z <- (score - 1) / sqrt(pmax(rowSums((W %*% Sf) * W), 1e-12))
  sig <- stats::p.adjust(stats::pnorm(-z), "BH") < alpha
  sig[is.na(sig)] <- FALSE
  n <- sum(sig)
  if (n < 2)
    return(list(curve = NULL, mean_fdr = NA_real_, fdr_top10 = NA_real_,
                n_sig = n, sig = sig))

  Ws <- W[sig, , drop = FALSE]
  o <- order(-score[sig])
  ks <- unique(pmax(1, round(pct / 100 * n)))

  # margin of each called cell at each cut: 1 SD of the estimation noise of its
  # score difference with the cells within 1% of the cut (observed order)
  h <- max(1, round(0.01 * n))
  CW <- Ws %*% Sf
  q <- rowSums(CW * Ws)
  margin <- lapply(ks, function(k) {
    cl <- o[seq_len(k)]
    wn <- o[max(1, k - h):min(n, k + h)]
    wbar <- colMeans(Ws[wn, , drop = FALSE])
    sqrt(pmax(q[cl] + mean(q[wn]) - 2 * as.vector(CW[cl, , drop = FALSE] %*% wbar), 0))
  })

  if (!is.null(seed)) set.seed(seed)
  E <- MASS::mvrnorm(n_draws, enrichment, Sigma)
  if (n_draws == 1) E <- matrix(E, nrow = 1)

  # per draw, the fraction of the observed top k below the k-th drawn score by
  # more than their margin
  Fd <- matrix(NA_real_, n_draws, length(ks))
  for (b in split(seq_len(n_draws), ceiling(seq_len(n_draws) / chunk))) {
    Eb <- E[b, , drop = FALSE]
    Eb[, small] <- 1
    Eb[Eb < 0] <- 1
    S <- Ws %*% t(Eb)
    bk <- matrix(apply(S, 2, function(x) sort(x, decreasing = TRUE)[ks]), nrow = length(ks))
    Fd[b, ] <- vapply(seq_along(ks), function(i) {
      cl <- o[seq_len(ks[i])]
      colMeans(sweep(-S[cl, , drop = FALSE], 2, bk[i, ], "+") > margin[[i]])
    }, numeric(length(b)))
  }
  fdr <- colMeans(Fd)
  pk <- 100 * ks / n
  list(curve = data.frame(pct = pk, k = ks, fdr = fdr),
       mean_fdr = mean(fdr),
       fdr_top10 = stats::approx(pk, fdr, 10, rule = 2)$y,
       n_sig = n, sig = sig)
}
