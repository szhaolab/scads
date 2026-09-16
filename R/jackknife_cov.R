#' Covariance of S-LDSC Enrichment Estimates from Block-Jackknife Delete Values
#'
#' Reconstructs the K x K sampling covariance of the per-topic enrichment
#' estimates from the block-jackknife delete-one values that S-LDSC writes when
#' run with `--print-delete-vals`.
#'
#' This replaces the annotation-correlation approximation
#' `Cov(e_k, e_k') ~ w_k w_k' Cor(A_k, A_k')`, which captures only peak overlap.
#' Enrichment estimates from separate S-LDSC runs are correlated for reasons the
#' annotation correlation cannot see: they share the GWAS summary statistics, the
#' LD reference panel, and the baseline-LD covariates. Two topics with disjoint
#' peak sets therefore still have correlated estimates.
#'
#' Requires each per-topic S-LDSC run to have been invoked with
#' `--print-delete-vals`, which writes `<out>.part_delete`: an
#' `n_blocks x n_annot` matrix of delete-one per-category coefficient (tau)
#' values. The focal topic annotation is assumed to be the first column, matching
#' the convention used elsewhere in this package when reading `.results` files.
#'
#' The blocks must correspond across the K runs. They will if every run used the
#' same summary statistics and the same reference SNPs; the function checks that
#' the delete-value matrices all have the same number of rows, but cannot verify
#' block identity beyond that. Confirm the K logs report the same retained-SNP
#' count.
#'
#' @param ldsc_res_dir Directory containing `k*_output/results/<trait>.results`
#'   and the matching `.part_delete` files.
#' @param trait Trait name used in the S-LDSC output filenames.
#' @param nTopics Number of topics (K).
#' @param M_annot Optional numeric vector of per-category SNP counts (length
#'   `n_annot`) used to convert delete-one coefficients to per-category h2. If
#'   `NULL` (default), these are read from the `.results` file's `Prop._SNPs`
#'   column scaled by the total SNP count.
#' @return A `K x K` covariance matrix of the enrichment estimates, or `NULL` if
#'   the delete-value files are not present (callers should then fall back to the
#'   annotation-correlation approximation).
#' @export
ldsc_jackknife_cov <- function(ldsc_res_dir, trait, nTopics, M_annot = NULL) {

  delete_list <- vector("list", nTopics)
  prop_snps   <- numeric(nTopics)

  for (k in seq_len(nTopics)) {
    base_dir <- file.path(ldsc_res_dir, paste0("k", k, "_output"), "results")
    res_f    <- file.path(base_dir, paste0(trait, ".results"))
    del_f    <- file.path(base_dir, paste0(trait, ".part_delete"))

    if (!file.exists(del_f)) {
      warning("No .part_delete for topic k", k, " at ", del_f,
              ". Re-run S-LDSC with --print-delete-vals to use the jackknife ",
              "covariance; falling back to the annotation-correlation ",
              "approximation.")
      return(NULL)
    }
    if (!file.exists(res_f)) {
      warning("Missing ", res_f)
      return(NULL)
    }

    res <- read.table(res_f, header = TRUE, sep = "\t", check.names = FALSE)
    prop_snps[k] <- res$`Prop._SNPs`[1]

    # n_blocks x n_annot matrix of delete-one per-category coefficients
    delete_list[[k]] <- as.matrix(utils::read.table(del_f))
  }

  nb <- vapply(delete_list, nrow, integer(1))
  if (length(unique(nb)) != 1L) {
    warning("Delete-value matrices have differing block counts (",
            paste(nb, collapse = ", "), "); blocks do not correspond across ",
            "runs. Falling back to the annotation-correlation approximation.")
    return(NULL)
  }
  n_blocks <- nb[1]

  # Per-block enrichment for the focal annotation of each run.
  #   cat_(b)  = M * tau_(b)          (per-category h2)
  #   tot_(b)  = sum_j cat_j,(b)      (total h2)
  #   E_k,(b)  = (cat_focal,(b) / tot_(b)) / prop_snps_k
  E_blocks <- matrix(NA_real_, nrow = n_blocks, ncol = nTopics)
  for (k in seq_len(nTopics)) {
    dv <- delete_list[[k]]
    M_k <- if (is.null(M_annot)) rep(1, ncol(dv)) else M_annot
    if (length(M_k) != ncol(dv)) {
      warning("M_annot length (", length(M_k), ") does not match delete-value ",
              "columns (", ncol(dv), "); falling back.")
      return(NULL)
    }
    cat_b <- sweep(dv, 2, M_k, FUN = "*")
    tot_b <- rowSums(cat_b)
    E_blocks[, k] <- (cat_b[, 1] / tot_b) / prop_snps[k]
  }

  if (anyNA(E_blocks) || any(!is.finite(E_blocks))) {
    warning("Non-finite per-block enrichment values; falling back.")
    return(NULL)
  }

  # Standard block-jackknife covariance.
  Ebar  <- colMeans(E_blocks)
  Ecen  <- sweep(E_blocks, 2, Ebar, FUN = "-")
  Sigma <- ((n_blocks - 1) / n_blocks) * crossprod(Ecen)

  dimnames(Sigma) <- list(paste0("k", seq_len(nTopics)),
                          paste0("k", seq_len(nTopics)))
  Sigma
}
