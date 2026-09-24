# Colon Fig6 (k=15): Laplace vs stored ns=1000 de_res$z, real baseline, full peaks.
suppressWarnings(suppressMessages({library(Matrix); library(devtools)}))
SC <- "/dartfs/rc/lab/S/Szhao/simingz/scads-revision/de_speedup/scads"; devtools::load_all(SC, quiet=TRUE)
D <- "/dartfs/rc/lab/S/Szhao/liyang/SCSviaTM/analysis/part1/scads_20250528_hg19_r2_k15_gc_new6"
OUT <- "/dartfs/rc/lab/S/Szhao/simingz/scads-revision/de_speedup"
r <- readRDS(file.path(D,"run_fastTopics_res.rds")); zref <- r$de_res$z
cm <- readRDS(file.path(D,"scATAC_Hickey_colon_epi_44K_counts_matrix.rds"))
if (nrow(cm) > ncol(cm)) cm <- Matrix::t(cm)
fit <- r$fastTopics_fit; f0 <- r$baseline
cat("counts:", paste(dim(cm),collapse=" x "), "| fit$F:", paste(dim(fit$F),collapse=" x "), "| k:", ncol(fit$F), "\n")
cat(sprintf("median reads/cell: %.0f\n", median(Matrix::rowSums(cm))))
stopifnot(ncol(cm)==nrow(fit$F))
rownames(cm) <- rownames(fit$L); colnames(cm) <- rownames(fit$F)
s <- Matrix::rowSums(cm)
tl <- system.time(de <- de_analysis2(fit, cm, s=s, lfc.stat="vsnull", lfc.method="laplace",
      shrink.method="none", control=list(nc=8, minval=1e-50), f0=f0, verbose=TRUE))[3]
zl <- de$z
cat(sprintf("\nLaplace wall-time (%d peaks, k=%d): %.2f min\n", nrow(zref), ncol(zref), tl/60))
bin<-function(z){z[is.na(z)]<- -Inf;p<-1-pnorm(z);apply(p,2,function(c) as.integer(p.adjust(c,"fdr",n=6e6)<0.05))}
Br<-bin(zref);Bl<-bin(zl); jac<-sapply(1:ncol(Br),function(t){a<-Br[,t];b<-Bl[,t];u<-sum(a|b);if(u==0)NA else sum(a&b)/u})
z1<-as.numeric(zref);z2<-as.numeric(zl);ok<-is.finite(z1)&is.finite(z2)
cat(sprintf("z Spearman %.4f | per-call FDR agreement %.4f\n", cor(z1[ok],z2[ok],method="spearman"), mean(Br==Bl)))
cat(sprintf("Pmat Jaccard: min %.3f median %.3f mean %.3f\n", min(jac,na.rm=T),median(jac,na.rm=T),mean(jac,na.rm=T)))
cat(sprintf("total nsig ref %d vs lap %d (ratio %.3f)\n", sum(Br),sum(Bl),sum(Bl)/sum(Br)))
saveRDS(list(zref=zref,zl=zl,tl=tl,jac=jac), file.path(OUT,"bench_colon.rds")); cat("saved\n")
