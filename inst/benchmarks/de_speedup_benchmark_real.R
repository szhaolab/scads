# Real-data test: run Laplace DE on a k=25 hema run and compare to the STORED
# ns=1000 MCMC z (the published reference) -- no MCMC rerun needed.
suppressWarnings(suppressMessages({library(Matrix); library(devtools)}))
SC <- "/dartfs/rc/lab/S/Szhao/simingz/scads-revision/de_speedup/scads"
devtools::load_all(SC, quiet=TRUE)
D <- "/dartfs/rc/lab/S/Szhao/liyang/SCSviaTM/analysis/part1/scads_20251030_eczema_gc"
OUT <- "/dartfs/rc/lab/S/Szhao/simingz/scads-revision/de_speedup"
r <- readRDS(file.path(D,"run_fastTopics_res.rds"))
zref <- r$de_res$z                       # stored ns=1000 MCMC reference (peaks x k)
cm <- readRDS(file.path(D,"scATAC_Hema1_33K_counts_matrix.rds"))
if (nrow(cm) > ncol(cm)) cm <- Matrix::t(cm)                 # cells x peaks
fit <- r$fastTopics_fit
cat("counts (cells x peaks):", paste(dim(cm),collapse=" x "), "| fit$F:", paste(dim(fit$F),collapse=" x "),
    "| zref:", paste(dim(zref),collapse=" x "), "| k:", ncol(fit$F), "\n")
stopifnot(ncol(cm)==nrow(fit$F))
rownames(cm) <- rownames(fit$L); colnames(cm) <- rownames(fit$F)
f0 <- r$baseline; s <- Matrix::rowSums(cm)
tl <- system.time(de_l <- de_analysis2(fit, cm, s=s, lfc.stat="vsnull", lfc.method="laplace",
      shrink.method="none", control=list(nc=8, minval=1e-50), f0=f0, verbose=TRUE))[3]
cat(sprintf("\nLaplace wall-time on full real data (%d peaks, k=%d): %.1f min\n",
    nrow(zref), ncol(zref), tl/60))
zl <- de_l$z
z1<-as.numeric(zref); z2<-as.numeric(zl); ok<-is.finite(z1)&is.finite(z2)
bin <- function(z){ p<-1-pnorm(z); apply(p,2,function(c) (p.adjust(c,"fdr",n=6e6)<0.05)+0) }
Br<-bin(zref); Bl<-bin(zl)
jac<-sapply(1:ncol(Br),function(t){a<-Br[,t];b<-Bl[,t];sum(a&b)/max(1,sum(a|b))})
cat(sprintf("z Pearson %.4f | z Spearman %.4f\n", cor(z1[ok],z2[ok]), cor(z1[ok],z2[ok],method="spearman")))
cat(sprintf("Pmat Jaccard per topic: min %.3f median %.3f mean %.3f\n", min(jac),median(jac),mean(jac)))
cat(sprintf("total nsig ref %d vs laplace %d (ratio %.3f)\n", sum(Br),sum(Bl),sum(Bl)/max(1,sum(Br))))
saveRDS(list(zref=zref, zl=zl, tl=tl, jac=jac, nsig_ref=colSums(Br), nsig_lap=colSums(Bl)),
        file.path(OUT,"bench_real_eczema.rds"))
cat("saved bench_real_eczema.rds\n")
