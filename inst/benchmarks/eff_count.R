# Derive when the Gaussian (laplace) approximation is reliable, from existing runs.
# The governing quantity is the effective per-topic count at a peak:
#   n_jk = F_jk * sum_i s_i L_ik   (expected counts attributable to topic k at peak j).
# Bin laplace-vs-MCMC z agreement by n_jk to find the reliability threshold.
suppressWarnings(suppressMessages({library(Matrix)}))
OUT <- "/dartfs/rc/lab/S/Szhao/simingz/scads-revision/de_speedup"
runs <- list(
  sim_k5 = list(dir="/dartfs/rc/lab/S/Szhao/liyang/SCSviaTM/analysis/Y_sim2/simu5_20250330_k3_final4",
                cnt="scATAC_SIM_k5_counts_matrix.rds", bench="bench_sim_vs_stored.rds"),
  eczema_k25 = list(dir="/dartfs/rc/lab/S/Szhao/liyang/SCSviaTM/analysis/part1/scads_20251030_eczema_gc",
                cnt="scATAC_Hema1_33K_counts_matrix.rds", bench="bench_real_eczema.rds"),
  colon_k15 = list(dir="/dartfs/rc/lab/S/Szhao/liyang/SCSviaTM/analysis/part1/scads_20250528_hg19_r2_k15_gc_new6",
                cnt="scATAC_Hickey_colon_epi_44K_counts_matrix.rds", bench="bench_colon.rds"))
brks <- c(0,1,2,5,10,20,50,100,Inf)
for (nm in names(runs)) {
  R <- runs[[nm]]; bf <- file.path(OUT, R$bench)
  if (!file.exists(bf)) { cat(nm, ": bench not ready, skip\n"); next }
  r <- readRDS(file.path(R$dir,"run_fastTopics_res.rds"))
  b <- readRDS(bf); zref <- b$zref; zl <- b$zl
  Fm <- r$de_res$F; L <- r$Lmat
  cm <- readRDS(file.path(R$dir,R$cnt)); if (nrow(cm)>ncol(cm)) cm <- Matrix::t(cm)
  s <- Matrix::rowSums(cm); rm(cm); gc()
  lsk <- as.numeric(Matrix::colSums(s * L))          # sum_i s_i L_ik, length k
  # effective count matrix n_jk = F_jk * lsk[k]  (F rows may be subset for eczema? no, full)
  if (nrow(Fm) != nrow(zref)) { Fm <- Fm[1:nrow(zref),,drop=FALSE] }
  neff <- sweep(Fm, 2, lsk, "*")
  # align: zref/zl are peaks x k, same order as Fm
  z1 <- as.numeric(zref); z2 <- as.numeric(zl); ne <- as.numeric(neff)
  ok <- is.finite(z1) & is.finite(z2) & is.finite(ne)
  z1<-z1[ok]; z2<-z2[ok]; ne<-ne[ok]
  # per-call significance agreement at genome-wide-ish z threshold 4.5, binned by n_eff
  bincut <- cut(ne, brks)
  agr <- tapply(seq_along(ne), bincut, function(ix) mean((z1[ix]>4.5)==(z2[ix]>4.5)))
  cnt <- table(bincut)
  sp  <- tapply(seq_along(ne), bincut, function(ix) if(length(ix)>5) cor(z1[ix],z2[ix],method="spearman") else NA)
  cat(sprintf("\n== %s (k=%d, median reads/cell %.0f) ==\n", nm, ncol(L), median(s)))
  cat("n_eff bin        frac_calls   sig-agree@4.5   zSpearman\n")
  for (i in seq_along(levels(bincut))) {
    lv<-levels(bincut)[i]
    cat(sprintf("  %-14s  %6.1f%%      %6.3f          %6.3f\n",
      lv, 100*cnt[lv]/sum(cnt), agr[lv], sp[lv])) }
}
cat("\nRule of thumb: agreement is high once n_eff (F_jk * sum_i s_i L_ik) exceeds the\nthreshold where the bins stabilize.\n")
