# scads 0.4.0

## New

* `cell_type_heterogeneity()` tests whether cells within a cell type differ in
  disease relevance, giving one p-value per cell type. The score is linear in the
  topic enrichments with weights summing to one, so the within-type variance is a
  quadratic form whose null does not depend on the unknown common enrichment
  level; the p-value is exact by Imhof inversion rather than permutation. See
  `vignette("HeterogeneityTest")`.

* `ldsc_jackknife_cov()` is rewritten and validated. It reconstructs the K x K
  covariance of the enrichment estimates from S-LDSC `--print-delete-vals`
  output, and now reproduces the reported `Enrichment` and `Enrichment_std_error`
  to four decimal places. Two corrections were needed: overlap-corrected
  per-category heritability (S-LDSC's `_overlap_output`) rather than the
  uncorrected ratio, and restriction to MAF 5-50% variants when `--frqfile-chr`
  is supplied. The function reports the acceptance check in `$check` and `$pass`,
  and no longer requires an opt-in flag.

  It takes `baseline_prefix` and `frq_prefix`, which must match the
  `--ref-ld-chr` and `--frqfile-chr` of the runs being reconstructed. A 25-topic
  run takes about three minutes.

* `get_cs()` gains `baseline_prefix` and `frq_prefix`. When both are supplied it
  builds the covariance with `ldsc_jackknife_cov()`, which is now the recommended
  route for the cell score variance. Supplying `Sigma` directly still works, and
  omitting all three falls back to the annotation-correlation approximation as
  before, so existing calls are unaffected.

## Notes

* `CompQuadForm` and `data.table` are new optional dependencies. Without
  `CompQuadForm` the quadratic-form p-value falls back to a Satterthwaite
  moment match; `data.table` is required by `ldsc_jackknife_cov()`.

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
