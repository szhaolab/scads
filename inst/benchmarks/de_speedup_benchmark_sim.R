# Benchmark de_analysis2 LFC methods on real Fig2 sim data.
suppressWarnings(suppressMessages({library(Matrix); library(devtools)}))
SC <- "/dartfs/rc/lab/S/Szhao/simingz/scads-revision/de_speedup/scads"
devtools::load_all(SC, quiet=TRUE)
D <- "/dartfs/rc/lab/S/Szhao/liyang/SCSviaTM/analysis/Y_sim2/simu5_20250330_k3_final4"
OUT <- "/dartfs/rc/lab/S/Szhao/simingz/scads-revision/de_speedup"
r <- readRDS(file.path(D,"run_fastTopics_res.rds"))
cm <- readRDS(file.path(D,"scATAC_SIM_k5_counts_matrix.rds"))   # the matrix the fit was built on
if (nrow(cm) > ncol(cm)) cm <- Matrix::t(cm)                    # cells x peaks
fit <- r$fastTopics_fit
cat("counts (cells x peaks):", paste(dim(cm),collapse=" x "), "| fit$F:", paste(dim(fit$F),collapse=" x "), "\n")
stopifnot(ncol(cm) == nrow(fit$F))
f0 <- 1e-7                       # constant background shared by all methods (exact comparison)
s <- Matrix::rowSums(cm)
set.seed(42); idx <- sort(sample(ncol(cm), min(30000, ncol(cm))))
cms <- cm[, idx]
fit_sub <- fit; fit_sub$F <- fit$F[idx,,drop=FALSE]            # keep dims consistent; F is refit anyway
# same matrix/order as the fit was built on: force dimnames to match exactly
rownames(cms) <- rownames(fit$L)
colnames(cms) <- rownames(fit_sub$F)
run <- function(method, ns=1000){
  t <- system.time(de <- de_analysis2(fit_sub, cms, s=s, lfc.stat="vsnull", lfc.method=method,
        shrink.method="none", control=list(ns=ns, nc=8, minval=1e-50), f0=f0, verbose=FALSE))[3]
  list(z=de$z, t=as.numeric(t))
}
cat("\n== mcmc ns=1000 (reference) ==\n"); ref <- run("mcmc",1000)
cat("== mcmc ns=200 ==\n");             m200<- run("mcmc",200)
cat("== laplace ==\n");                 lap <- run("laplace")
bin <- function(z){ p<-1-pnorm(z); apply(p,2,function(c) (p.adjust(c,"fdr",n=6e6)<0.05)+0) }
Bref<-bin(ref$z)
cmp <- function(nm,o){ z1<-as.numeric(ref$z);z2<-as.numeric(o$z);ok<-is.finite(z1)&is.finite(z2)
  B<-bin(o$z); jac<-sapply(1:ncol(Bref),function(t){a<-Bref[,t];b<-B[,t];sum(a&b)/max(1,sum(a|b))})
  cat(sprintf("%-12s time %6.1fs (%.1fx)  zcorr %.4f  zSpear %.4f  Jac[%s]  nsig ref/this [%s]/[%s]\n",
    nm,o$t,ref$t/o$t,cor(z1[ok],z2[ok]),cor(z1[ok],z2[ok],method="spearman"),
    paste(sprintf("%.2f",jac),collapse=","),paste(colSums(Bref),collapse=","),paste(colSums(B),collapse=","))) }
cat("\n===== RESULTS (", length(idx), "peaks, k=", ncol(ref$z), ") =====\n")
cat(sprintf("%-12s time %6.1fs (ref)\n","mcmc ns1000",ref$t))
cmp("mcmc ns200",m200); cmp("laplace",lap)
saveRDS(list(ref=ref,m200=m200,lap=lap,idx=idx), file.path(OUT,"bench_fig2.rds"))
cat("\nsaved bench_fig2.rds\n")
