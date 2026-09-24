#' Fast closed-form (Laplace) log-fold-change statistics for de_analysis2
#'
#' Replaces the per-peak Metropolis MCMC used to quantify uncertainty of the
#' topic-specific log-rates with a Laplace / Fisher-information approximation.
#' The MCMC in \code{compute_lfc_stats2} samples \eqn{g = \log f} under a flat
#' prior on \eqn{f}; its posterior covariance is well approximated by the
#' inverse Poisson Fisher information in log space,
#' \deqn{I_{kl}(j) = f_{jk} f_{jl} \sum_i L_{ik} L_{il} / u_i, \quad u_i = \sum_k L_{ik} f_{jk}.}
#' For \code{lfc.stat = "vsnull"} the posterior mean is
#' \eqn{\log f_{jk} - \log f0_{jk}} and the z-score is that divided by the
#' posterior SD \eqn{\sqrt{(I^{-1})_{kk}}}, matching the definition used by
#' \code{fastTopics:::compute_zscores} at \code{conf.level = 0.68}.
#'
#' This is deterministic (no Monte Carlo noise) and typically 1-2 orders of
#' magnitude faster than the ns=1000 MCMC.
#'
#' @keywords internal
compute_lfc_stats_laplace <- function(X, F, L, f0, lfc.stat = "vsnull",
                                      nc = 1, nsplit = 100, ridge = 1e-10,
                                      verbose = TRUE) {
  if (!all(lfc.stat == "vsnull"))
    stop("compute_lfc_stats_laplace() currently supports only lfc.stat = \"vsnull\".")
  m <- nrow(F); k <- ncol(F)
  # per-block worker over a set of peak columns
  block <- function(js) {
    est  <- matrix(0, length(js), k)
    sdg  <- matrix(0, length(js), k)
    for (a in seq_along(js)) {
      j  <- js[a]
      fj <- F[j, ]
      u  <- as.numeric(L %*% fj)
      u[u < 1e-15] <- 1e-15
      A  <- L / sqrt(u)                 # n x k
      I  <- crossprod(A)                # k x k : sum_i L_ik L_il / u_i
      I  <- outer(fj, fj) * I           # Fisher info in g = log f
      diag(I) <- diag(I) + ridge * max(diag(I))
      Iinv <- tryCatch(solve(I), error = function(e) MASS::ginv(I))
      sdg[a, ] <- sqrt(pmax(diag(Iinv), 0))
      est[a, ] <- log(fj) - log(f0[j, ])
    }
    list(est = est, sdg = sdg)
  }
  cols <- parallel::splitIndices(m, min(m, nsplit))
  if (nc > 1) {
    if (verbose) op <- pbapply::pboptions(type = "txt", txt.width = 70)
    else         op <- pbapply::pboptions(type = NULL)
    ans <- pbapply::pblapply(cl = nc, cols, block)
    pbapply::pboptions(op)
  } else {
    ans <- lapply(cols, block)
  }
  est <- matrix(0, m, k); sdg <- matrix(0, m, k)
  dimnames(est) <- dimnames(F); dimnames(sdg) <- dimnames(F)
  for (i in seq_along(cols)) { est[cols[[i]], ] <- ans[[i]]$est
                               sdg[cols[[i]], ] <- ans[[i]]$sdg }
  postmean <- est
  lower <- postmean - sdg
  upper <- postmean + sdg
  z <- fastTopics:::compute_zscores(postmean, lower, upper)
  list(ar       = matrix(NA_real_, m, k, dimnames = dimnames(F)),
       est      = est / log(2),
       postmean = postmean / log(2),
       lower    = lower / log(2),
       upper    = upper / log(2),
       z        = z)
}
