# Speeding up the DE step (`de_analysis2`) in scads

## The bottleneck

`run_fastTopics()` → `de_analysis2()` quantifies per-peak, per-topic
log-fold-change uncertainty by running a random-walk **Metropolis MCMC**
(`ns = 1000` samples) **for every peak** (`simulate_posterior_poisson_sparse_rcpp`).
On a real k=25 haematopoiesis run (452,004 peaks × 33,819 cells) this step takes
**3+ hours and timed out at a 3 h wall limit**. Only the resulting **z-scores**
flow downstream: `compute_topic_pvalues()` turns them into the binarized topic
annotation `Pmat` (FDR < 0.05), which is what S-LDSC consumes.

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
