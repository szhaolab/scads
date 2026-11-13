# scads 0.8.0

* `get_cs()`: the cell score uses the raw S-LDSC enrichments; ASH shrinkage and
  the `ash_res` output are removed.
* New `rank_fdr()`, also returned by `get_cs()` as `rank_fdr`: resolvable-error
  rank FDR of the significant cells.
* `get_cs()` gains `fdr_alpha`, `fdr_draws` and `fdr_seed`; `cs_dat` gains `W`
  and `drop`.
* `Sigma_used` from the annotation-correlation approximation includes floored
  topics.
* `get_cs()` gains `chrs`, passed through by `scads()`; runs on a subset of
  chromosomes no longer fail.

# scads 0.7.0

* `get_cs()` judges topic size by 0.01 cM block coverage (`min_pct_blocks`,
  default 1.7%; Tashman et al. 2021) instead of 0.5% of SNPs.
* Small topics are reported in a warning and floored (`small_annot = "floor"`,
  default) or kept (`"warn"`); `scads()` passes both settings through.
* `get_cs()` returns `topic_qc`.
* Floored topics contribute no variance when `Sigma` is supplied.
* `min_prop_snps` is deprecated. Without `baseline_prefix` or `block_coverage`,
  `get_cs()` falls back to the 0.5% rule.
* New `annot_block_coverage()`.

# scads 0.6.0

* `run_fastTopics()` and `scads()` gain `lfc.method = c("mcmc", "laplace")`;
  `"laplace"` is a faster approximation of the differential accessibility step.
  Benchmarks in `inst/benchmarks/`.

# scads 0.5.0

* New `cell_type_heterogeneity()`: test of heterogeneity within each cell type
  (`vignette("HeterogeneityTest")`).
* `ldsc_jackknife_cov()` rewritten: overlap-corrected heritability, MAF 5-50%
  variants, reproduces the S-LDSC `Enrichment` and its SE; takes
  `baseline_prefix` and `frq_prefix`; returns `$check` and `$pass`.
* `get_cs()` uses the jackknife covariance by default and gains
  `baseline_prefix` and `frq_prefix`, passed through by `scads()`.
* **Breaking:** `get_cs()` without `Sigma`, the reference prefixes or
  `use_approximation = TRUE` is an error.
* `CompQuadForm` and `data.table` added to Suggests.

# scads 0.3.0

* **Breaking:** `get_cs()` cell-score variance includes both off-diagonal terms
  of each topic pair; z-scores decrease.
* **Breaking:** `get_cs()` returns `p_cell` of length I, `NA` where the variance
  is zero or non-finite.
* New `ldsc_jackknife_cov()`: covariance of the enrichment estimates from S-LDSC
  block-jackknife delete values, for `get_cs(Sigma = )`.
* `get_cs()` gains `alternative` and `min_prop_snps`, and returns `var_cell` and
  `Sigma_used`.
