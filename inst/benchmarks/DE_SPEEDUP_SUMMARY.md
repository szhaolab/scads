# Speeding up the DE step (`de_analysis2`) in scads

## The bottleneck

`run_fastTopics()` → `de_analysis2()` quantifies per-peak, per-topic
log-fold-change uncertainty by running a random-walk **Metropolis MCMC**
(`ns = 1000` samples) **for every peak** (`simulate_posterior_poisson_sparse_rcpp`).
On a real k=25 haematopoiesis run (452,004 peaks × 33,819 cells) this step takes
**3+ hours and timed out at a 3 h wall limit**. Only the resulting **z-scores**
flow downstream: `compute_topic_pvalues()` turns them into the binarized topic
annotation `Pmat` (FDR < 0.05), which is what S-LDSC consumes.


## Why fastTopics uses MCMC (not a Gaussian approximation)

fastTopics samples the posterior of the topic log-rate `g = log f` (flat prior on
`f`) by MCMC because that posterior is **not Gaussian in the low-count regime that
dominates single-cell data**:

- **Skew at low effective counts.** The Poisson log-likelihood for a rate is far
  from quadratic when few counts inform it, so the posterior of `log f` is
  right-skewed; a symmetric Laplace fit mis-sizes the interval.
- **Boundary breakdown.** When a peak has ~0 counts in a topic, the MLE `f -> 0`,
  `log f -> -Inf`, and the Fisher information is unstable, so the normal SE is
  undefined/inflated. MCMC with a proper prior stays finite and returns a valid
  *asymmetric* HPD interval (the z-code reads `postmean/(postmean-lower)` off one
  HPD bound for exactly this reason).

MCMC is therefore robust across the whole dynamic range; the price is speed. This
predicts the Laplace approximation is accurate only when the **effective per-topic
count** `n_jk = F_jk * sum_i s_i L_ik` is large enough for the log-rate posterior
to be ~symmetric -- i.e. high read depth, well-expressed peaks, and not-too-many
topics (more topics split the reads thinner).

## Approaches tried

1. **Reduce `ns`** (fewer MCMC samples). Exposed already via `control$ns`.
2. **Laplace / Fisher-information closed form** (`lfc.method = "laplace"`, new).
   The MCMC samples `g = log f` under a flat prior on `f`; that posterior is
   approximately Gaussian, so its mean and SD have closed forms:
   - mean: `log F[j,k] − log f0[j,k]` (already computed in the F-fitting step);
   - SD: from the Poisson Fisher information in log space,
     `I_j[k,l] = F_jk F_jl Σ_i L_ik L_il / u_i`, `u_i = Σ_k L_ik F_jk`;
     posterior SD of topic k = `sqrt((I_j⁻¹)_kk)`.
   - `z = mean / SD`, matching `fastTopics:::compute_zscores` at conf.level 0.68.
   One `crossprod` + one k×k solve per peak instead of a 1000-step chain.
   Deterministic (no Monte Carlo noise).

Both are behind a switch; **`lfc.method` defaults to `"mcmc"`, so current
behaviour is unchanged unless the user opts in.**

## Results

### Simulation (Fig2 sim, full 259,941 peaks, k=5)

Compared against the **existing** published run's stored ns=1000 z
(`run_fastTopics_res.rds$de_res$z`), using the same gc baseline:

| metric | value |
|--------|-------|
| Laplace DE wall-time | **0.60 min** (full 260k peaks) |
| z Spearman vs stored ns=1000 z | 0.929 |
| per-call FDR-significance agreement | **0.991** |
| Pmat Jaccard per topic | 0.963–0.998, **median 0.983** |
| total significant peaks (laplace / ref) | 0.837M / 0.846M (ratio 0.990) |

Annotations are **96–99.8% identical** to the published ns=1000 run — the
"minimally changed" bar is met at k=5.

A separate speed comparison of the three LFC methods on a 30k-peak subset
(mcmc ns=1000 vs ns=200 vs laplace) gave **laplace 44.7× faster** (5.5s vs
244s) with higher rank fidelity (Spearman 0.857) than ns=200 (0.758, 4.2×).
Reducing `ns` below ~500 makes z **magnitudes** unstable (HPD SD from few
samples is noisy); laplace is deterministic.

### Real data (eczema haematopoiesis run, full 452,004 peaks, k=25), vs the published ns=1000 z

| metric | value |
|--------|-------|
| Laplace DE wall-time | **8.4 min** (MCMC: 3 h+, did not finish) |
| z Spearman vs published z | 0.799 |
| per-call FDR-significance agreement | 0.956 |
| Pmat Jaccard per topic | min 0.642, **median 0.840**, mean 0.837 |
| total significant peaks (laplace / ref) | 2.98M / 2.77M (ratio 1.08) |

The Laplace annotation overlaps the published one at median Jaccard 0.84 and is
within 8% on total size, with 95.6% of significance calls identical. Agreement is near-perfect at k=5 (median Jaccard 0.98) and good at k=25
(median 0.84); the gap reflects greater topic collinearity at high k, which the
full k×k Fisher inverse partly but not fully absorbs.

### Real data: colon Fig6 (full 499,517 peaks, k=15), vs published ns=1000 z

| metric | value |
|--------|-------|
| Laplace DE wall-time | 9.2 min |
| z Spearman | 0.886 |
| per-call FDR agreement | 0.954 |
| Pmat Jaccard per topic | min 0.643, median 0.894, mean 0.872 |
| total sig-peak ratio (laplace/ref) | 0.895 |

### Read-depth sweep (Fig2 sim, 40k peaks, k=5; binomial thinning)

MCMC ns=1000 reference recomputed at each depth; Laplace vs that reference.

| thinning | median reads/cell | MCMC time | Laplace time | speedup | Jaccard median | FDR agree |
|----------|-------------------|-----------|--------------|---------|----------------|-----------|
| 1.00 | 19,631 | 267 s | 6 s | 49x | 0.991 | 0.995 |
| 0.30 | 5,898  | 93 s  | 4 s | 22x | 0.981 | 0.990 |
| 0.10 | 1,962  | 40 s  | 4 s | 9x  | 0.908 | 0.948 |

At k=5 the approximation holds up even at ~2k reads/cell (Jaccard 0.91). Note the
speedup shrinks at low depth: the sparse MCMC is cheaper when there are fewer
nonzeros, while Laplace time is roughly flat -- so Laplace helps most at high
depth / high nnz.

### When is the Gaussian approximation good?

Agreement binned by effective count `n_jk = F_jk * sum_i s_i L_ik`
(sig-agreement at z>4.5):

| n_eff bin | sim k=5 | colon k=15 | eczema k=25 |
|-----------|---------|------------|-------------|
| (0,1]     | 1.00    | 1.00       | 1.00  (both call non-sig) |
| (10,20]   | 1.00    | 0.99       | 0.97 |
| (20,50]   | 1.00    | 0.91       | 0.80 |
| (50,100]  | 0.85    | 0.84       | 0.85 |
| (100,Inf] | 1.00    | 0.98       | 0.93 |
| % of calls with n_eff>100 | 63% | 18% | 10% |

Disagreement concentrates in the **boundary range n_eff ~ 20-100**, where MCMC's
capture of posterior skew flips marginally-significant calls; it is small both far
below (both non-significant) and far above (both significant). The *fraction* of
peak-topics in that boundary range grows with k (topics split the reads), which is
why overall annotation concordance falls from ~0.98 (k=5) to ~0.89 (k=15) to
~0.84 (k=25) at comparable read depth.

**Practical guidance for `lfc.method = "laplace"`:**
- **Recommended** when k is small-to-moderate (k <= ~15) at typical depth
  (>= ~10k reads/cell): annotations 89-99% concordant, 9-49x faster.
- **Use with a caveat** at high k (>= ~25) or low depth: expect ~84% annotation
  concordance, with disagreement on borderline peaks (n_eff ~20-100). Prefer MCMC
  if those borderline calls are load-bearing.
- Read depth alone matters less than k: at k=5, thinning to 2k reads/cell still
  gave Jaccard 0.91.

## Recommendation

Adopt **`lfc.method = "laplace"`** as an opt-in fast path (keep MCMC the
default). It turns a 3 h+ step into ~8 min (~20×+) on real k=25 data and is
deterministic. Annotations differ from MCMC at the ~16% (1 − Jaccard) boundary
level; because S-LDSC enrichment aggregates over tens of thousands of peaks and
the annotation sizes match within 8%, the downstream enrichment change is
expected to be small.

**Remaining check before merge:** run per-topic S-LDSC on a laplace-derived vs
mcmc-derived annotation for one trait and confirm the enrichment estimates
agree within their standard errors. (Not yet run — the definitive
"minimally changed" test at the enrichment level.)

## Files
- `R/de_speedup.R` — `compute_lfc_stats_laplace()`
- `R/fastTopics_aux_func.R` — `de_analysis2(..., lfc.method=)` switch
- `R/run_fastTopics.R` — `run_fastTopics(..., lfc.method=)` passthrough
- `inst/benchmarks/de_speedup_benchmark_{sim,real}.R` — benchmark scripts
