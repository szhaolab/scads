# Sim comparison against the EXISTING ns=1000 run (stored de_res$z) with the real
# gc baseline -- mirrors the real-data recipe. No MCMC rerun.
suppressWarnings(suppressMessages({library(Matrix); library(devtools)}))
SC <- "/dartfs/rc/lab/S/Szhao/simingz/scads-revision/de_speedup/scads"
devtools::load_all(SC, quiet=TRUE)
D <- "/dartfs/rc/lab/S/Szhao/liyang/SCSviaTM/analysis/Y_sim2/simu5_20250330_k3_final4"
OUT <- "/dartfs/rc/lab/S/Szhao/simingz/scads-revision/de_speedup"
r <- readRDS(file.path(D,"run_fastTopics_res.rds"))
zref <- r$de_res$z                              # stored ns=1000 MCMC z (real baseline)
cm <- readRDS(file.path(D,"scATAC_SIM_k5_counts_matrix.rds"))   # matrix the fit was built on
if (nrow(cm) > ncol(cm)) cm <- Matrix::t(cm)                    # cells x peaks
fit <- r$fastTopics_fit; f0 <- r$baseline
cat("counts:", paste(dim(cm),collapse=" x "), "| fit$F:", paste(dim(fit$F),collapse=" x "),
    "| zref:", paste(dim(zref),collapse=" x "), "| baseline:", paste(dim(as.matrix(f0)),collapse=" x "), "\n")
stopifnot(ncol(cm)==nrow(fit$F), nrow(as.matrix(f0))==nrow(fit$F))
rownames(cm) <- rownames(fit$L); colnames(cm) <- rownames(fit$F)   # same order, force names
s <- Matrix::rowSums(cm)
tl <- system.time(de_l <- de_analysis2(fit, cm, s=s, lfc.stat="vsnull", lfc.method="laplace",
      shrink.method="none", control=list(nc=8, minval=1e-50), f0=f0, verbose=TRUE))[3]
zl <- de_l$z
cat(sprintf("\nLaplace wall-time full sim (%d peaks, k=%d): %.2f min\n", nrow(zref), ncol(zref), tl/60))
bin<-function(z){ z[is.na(z)]<- -Inf; p<-1-pnorm(z); apply(p,2,function(c) as.integer(p.adjust(c,"fdr",n=6e6)<0.05)) }
Br<-bin(zref); Bl<-bin(zl)
jac<-sapply(1:ncol(Br),function(t){a<-Br[,t];b<-Bl[,t];u<-sum(a|b); if(u==0) NA else sum(a&b)/u})
z1<-as.numeric(zref);z2<-as.numeric(zl);ok<-is.finite(z1)&is.finite(z2)
cat(sprintf("z Spearman(finite) %.4f | per-call FDR agreement %.4f\n", cor(z1[ok],z2[ok],method="spearman"), mean(Br==Bl)))
cat(sprintf("Pmat Jaccard per topic: %s\n", paste(sprintf("%.3f",jac),collapse=",")))
cat(sprintf("Jaccard min %.3f median %.3f mean %.3f\n", min(jac,na.rm=T),median(jac,na.rm=T),mean(jac,na.rm=T)))
cat(sprintf("nsig ref: %s\n", paste(colSums(Br),collapse=",")))
cat(sprintf("nsig lap: %s\n", paste(colSums(Bl),collapse=",")))
cat(sprintf("total nsig ref %d vs lap %d (ratio %.3f)\n", sum(Br),sum(Bl),sum(Bl)/sum(Br)))
saveRDS(list(zref=zref,zl=zl,tl=tl,jac=jac,nsig_ref=colSums(Br),nsig_lap=colSums(Bl)),
        file.path(OUT,"bench_sim_vs_stored.rds"))
cat("saved bench_sim_vs_stored.rds\n")
