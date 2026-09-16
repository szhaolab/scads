# scads 0.3.0

Revision of NCOMMS-26-026605-T. Changes to the cell-level variance and p-value.

## Breaking

* `get_cs()` now computes the cell score variance as the **full double sum**
  over ordered topic pairs, `Var(sum_k a_k e_k) = sum_k sum_k' a_k a_k' Cov(e_k, e_k')`.
  Previously the correlation matrix was masked to its upper triangle, which
  omits the factor of 2 on the off-diagonal terms and understates the variance
  whenever topic annotations are positively correlated. Cell-level z-scores
  decrease and p-values increase; the effect grows with the mean annotation
  correlation (roughly 17% at rho = 0.05, 34% at rho = 0.3 for K = 25).

* `get_cs()` returns `p_cell` of length I, with `NA` where the variance is zero
  or non-finite. Previously non-finite values were dropped, so the returned
  vector could be shorter than `cs` and `z_cell` and every subsequent p-value
  was attributed to the wrong cell.

## New

* `ldsc_jackknife_cov()` reconstructs the K x K covariance of the enrichment
  estimates from S-LDSC block-jackknife delete values (`--print-delete-vals`).
  Pass the result to `get_cs(Sigma = )` to replace the annotation-correlation
  approximation, which captures only peak overlap and ignores the correlation
  induced by the shared GWAS, LD reference and baseline-LD covariates.

* `get_cs(alternative = )` selects a two-sided (default, unchanged) or one-sided
  (`"greater"`) test. Under the two-sided test, significantly *depleted* cells
  also pass FDR.

* `get_cs()` returns `var_cell` and `Sigma_used`.

* `get_cs(min_prop_snps = )` exposes the annotation-size floor, previously
  hardcoded at 0.005. The flooring behaviour is now documented: floored topics
  remain in the weight denominator, diluting the score toward 1 rather than
  renormalising over retained topics.
