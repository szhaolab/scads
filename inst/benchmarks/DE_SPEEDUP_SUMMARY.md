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

### Simulation (Fig2 sim, 30,000-peak subset, k=5), vs mcmc ns=1000 reference

| method       | time   | speedup | z Spearman | sig-agree @z>2/4/6 |
|--------------|--------|---------|-----------|--------------------|
| mcmc ns=1000 | 244 s  | 1×      | (ref)     | —                  |
| mcmc ns=200  | 58 s   | 4.2×    | 0.758     | 0.949/0.946/0.924  |
| **laplace**  | 5.5 s  | **44.7×** | **0.857** | 0.948/0.941/0.919 |

Laplace dominates ns-reduction: ~10× faster **and** higher rank fidelity.
Reducing `ns` below ~500 also makes z **magnitudes** unstable (HPD SD from few
samples is noisy), while Laplace is deterministic.

### Real data (eczema haematopoiesis run, full 452,004 peaks, k=25), vs the published ns=1000 z

| metric | value |
|--------|-------|
| Laplace DE wall-time | **8.4 min** (MCMC: 3 h+, did not finish) |
| z Spearman vs published z | 0.799 |
| per-call FDR-significance agreement | 0.956 |
| Pmat Jaccard per topic | min 0.642, **median 0.840**, mean 0.837 |
| total significant peaks (laplace / ref) | 2.98M / 2.77M (ratio 1.08) |

The Laplace annotation overlaps the published one at median Jaccard 0.84 and is
within 8% on total size, with 95.6% of significance calls identical. The small
size bias is topic-dependent (slightly conservative on the k=5 sim, slightly
liberal on the k=25 real run).

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
