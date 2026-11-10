# scads 0.8.0

## Changed

* **The cell score uses the raw S-LDSC enrichments.** `get_cs()` no longer
  shrinks the topic enrichments with ASH before forming `cs`; the score is
  `sum_k w_ik e_k` from the reported `Enrichment`, floored as before. The
  `ash_res` element is removed from the output. On the revision simulations ASH
  shrinkage lowered the correlation with the truth, most where the causal topic
  had the largest standard error.
* `Sigma_used` from the annotation-correlation path is now the covariance of
  the reported enrichments for all topics (previously floored topics had zero
  rows), matching what the jackknife path returns.

## New

* `rank_fdr()` and a new `get_cs()` output, `rank_fdr`: the resolvable-error
  rank FDR of the significant cells (one-sided BH < `fdr_alpha`). Calling the
  top `k` significant cells by score, a call counts as false only if the cell
  is below the k-th score by more than one standard deviation of the estimation
  noise of its score difference with the cells at the cut-off (the cells within
  1% of position `k`); near-ties within that noise are not errors. The FDR is
  the average share of such calls when the enrichments are drawn from
  `N(e_hat, Sigma)`. It is returned as a curve over the top 1-50% of
  significant cells, with summaries `mean_fdr` and `fdr_top10`. The topics
  floored in the score (negative enrichment, or below the annotation-size
  threshold with `small_annot = "floor"`) are floored in the draws too.
  `get_cs()` gains `fdr_alpha` (0.05), `fdr_draws` (1000; 0 skips it) and
  `fdr_seed` (1); `cs_dat` gains the normalised cell topic weights `W` and the
  floored-topic indicator `drop`.
## Fixed

* `get_cs()` gains `chrs` (default `1:22`), passed to `ldsc_jackknife_cov()` and
  `annot_block_coverage()`, and `scads()` passes its `chrs` through. Previously
  a pipeline run on a subset of chromosomes (e.g. `scads(chrs = 10)`) failed in
  `get_cs()`, which read all 22 chromosomes for the covariance.

* On the revision simulations (20 settings, about 480 GWAS replicates) the
  rank FDR averaged over runs agrees with the actual share of resolvable errors
  within run-to-run noise, and that share tracks how often a run ranks the
  wrong cell population on top.

# scads 0.7.0

## Changed

* The annotation-size rule in `get_cs()` now follows Tashman et al. (2021,
  medRxiv 2021.03.13.21249938): a topic must overlap at least `min_pct_blocks`
  (default 1.7) percent of the genome's 0.01 cM blocks. This replaces the
  hard cut-off of 0.5% of SNPs. Block coverage reflects both an annotation's
  size and how clustered it is, and it is the metric that separated S-LDSC
  annotations with controlled type I error from those without. A topic with
  fewer than 0.5% of SNPs but spread across the genome is now kept; a larger
  but clustered topic can now be flagged.
* Topics below the threshold are always reported in a warning. They are
  floored to `e = 1` (`small_annot = "floor"`, default, as before) or kept
  (`small_annot = "warn"`). `scads()` passes both settings through.
* `get_cs()` returns `topic_qc`, one row per topic with `Prop._SNPs`,
  `pct_blocks` and the flooring decision.
* When `Sigma` is supplied (including the default jackknife covariance),
  floored topics now contribute no variance, as documented. Previously their
  rows of `Sigma` were kept, so cells loading on floored topics had an inflated
  variance and conservative p-values.
* `min_prop_snps` is deprecated. Supplying it restores the old rule. Without a
  `baseline_prefix` (or precomputed `block_coverage`) block coverage cannot be
  computed, and `get_cs()` falls back to the 0.5% rule with a message.

## New

* `annot_block_coverage()` computes, per topic, the percentage of 0.01 cM blocks
  overlapped, from the baseline annotation's genetic map and the thin topic
  annotations written by `run_sldsc()`.

# scads 0.6.0

## New

* `run_fastTopics()` and `scads()` gain `lfc.method = c("mcmc", "laplace")`.
  The default `"mcmc"` is unchanged. `"laplace"` replaces the per-peak
  random-walk MCMC in the DE step (the Step-1 runtime bottleneck) with a
  deterministic Laplace / Fisher-information approximation of the log-rate
  posterior for `lfc.stat = "vsnull"`. On a 452,000-peak, k=25 dataset the DE
  step drops from 3+ hours (did not finish) to ~8 minutes.

## When to use `"laplace"`

fastTopics uses MCMC because the log-rate posterior is skewed and
boundary-unstable at low counts, where a Gaussian approximation is poor. The
approximation is accurate when the effective per-topic count
`n_eff = F[peak, topic] * sum_i s_i L[i, topic]` is large. Benchmarks vs
`n_s = 1000` MCMC: topic annotations are 89-99% concordant (Jaccard) at
`k <= ~15` and `>= ~10,000` reads/cell, and ~84% concordant at `k = 25`
(disagreement on borderline-significant peaks). Recommended as an opt-in
speed-up for small-to-moderate `k`; keep the default `"mcmc"` at high `k` or low
depth when marginal calls matter. Benchmark scripts and full results are in
`inst/benchmarks/`.

# scads 0.5.0

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

* `get_cs()` now uses the jackknife covariance **by default**. It gains
  `baseline_prefix` and `frq_prefix`, which `scads()` passes through from its own
  arguments, so the covariance is built as part of the pipeline with no extra
  steps. Supplying `Sigma` directly still works.

  **Breaking:** calling `get_cs()` with none of `Sigma`, the reference prefixes,
  or `use_approximation = TRUE` is now an error rather than a silent fallback to
  the annotation-correlation approximation. Set `use_approximation = TRUE` to
  keep the previous behaviour.

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
