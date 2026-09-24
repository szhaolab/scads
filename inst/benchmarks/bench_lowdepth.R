# Low-read-depth test: binomially thin the Fig2 sim counts, then compare
# laplace vs mcmc ns=1000 (reference recomputed at the low depth). Isolates the
# DE step's sensitivity to read depth. 40k-peak subset for a tractable MCMC ref.
suppressWarnings(suppressMessages({library(Matrix); library(devtools)}))
SC <- "/dartfs/rc/lab/S/Szhao/simingz/scads-revision/de_speedup/scads"; devtools::load_all(SC, quiet=TRUE)
D <- "/dartfs/rc/lab/S/Szhao/liyang/SCSviaTM/analysis/Y_sim2/simu5_20250330_k3_final4"
OUT <- "/dartfs/rc/lab/S/Szhao/simingz/scads-revision/de_speedup"
r <- readRDS(file.path(D,"run_fastTopics_res.rds"))
cm <- readRDS(file.path(D,"scATAC_SIM_k5_counts_matrix.rds")); if (nrow(cm)>ncol(cm)) cm <- Matrix::t(cm)
fit <- r$fastTopics_fit; f0 <- r$baseline
rownames(cm) <- rownames(fit$L); colnames(cm) <- rownames(fit$F)
cat(sprintf("original median reads/cell: %.0f\n", median(Matrix::rowSums(cm))))
set.seed(1)
res <- list()
for (p in c(1.0, 0.30, 0.10)) {
  # binomial thinning of counts
  Xt <- cm; if (p < 1) { Xt@x <- rbinom(length(cm@x), cm@x, p) * 1.0; Xt <- drop0(Xt) }
  set.seed(42); idx <- sort(sample(ncol(Xt), 40000)); Xs <- Xt[, idx]; f0s <- as.matrix(f0)[idx,,drop=FALSE]
  fit_s <- fit; fit_s$F <- fit$F[idx,,drop=FALSE]
  colnames(Xs) <- rownames(fit_s$F); rownames(Xs) <- rownames(fit$L)
  s <- Matrix::rowSums(Xt)
  mdepth <- median(Matrix::rowSums(Xt))
  tm <- system.time(dm <- de_analysis2(fit_s, Xs, s=s, lfc.stat="vsnull", lfc.method="mcmc",
        shrink.method="none", control=list(ns=1000,nc=8,minval=1e-50), f0=f0s, verbose=FALSE))[3]
  tl <- system.time(dl <- de_analysis2(fit_s, Xs, s=s, lfc.stat="vsnull", lfc.method="laplace",
        shrink.method="none", control=list(nc=8,minval=1e-50), f0=f0s, verbose=FALSE))[3]
  bin<-function(z){z[is.na(z)]<- -Inf;p<-1-pnorm(z);apply(p,2,function(c) as.integer(p.adjust(c,"fdr",n=6e6)<0.05))}
  Bm<-bin(dm$z);Bl<-bin(dl$z); jac<-sapply(1:ncol(Bm),function(t){a<-Bm[,t];b<-Bl[,t];u<-sum(a|b);if(u==0)NA else sum(a&b)/u})
  z1<-as.numeric(dm$z);z2<-as.numeric(dl$z);ok<-is.finite(z1)&is.finite(z2)
  cat(sprintf("p=%.2f  median reads/cell %6.0f  mcmc %5.0fs laplace %4.0fs (%.0fx)  zSpear %.3f  FDRagree %.3f  Jac median %.3f  nsig m/l %d/%d\n",
      p, mdepth, tm, tl, tm/tl, cor(z1[ok],z2[ok],method="spearman"), mean(Bm==Bl),
      median(jac,na.rm=T), sum(Bm), sum(Bl)))
  res[[as.character(p)]] <- list(p=p, mdepth=mdepth, zm=dm$z, zl=dl$z, F=dl$est, jac=jac)
}
saveRDS(res, file.path(OUT,"bench_lowdepth.rds")); cat("saved bench_lowdepth.rds\n")
