#' scads: A function for using scATAC-seq data and GWAS sumstats to calculate disease cell score
#'
#' @param count_matrix A peak-by-cell count matrix.
#' @param nTopics The number of topics (clusters) to fit. Default is 10.
#' @param n_s Number of simulations for fastTopics-DE 
#' @param n_c Number of cores to use for parallel
#' @param sumstats_dir Directory to location of GWAS sumstats
#' @param gwas_nsamps Sample size of GWAS
#' @param gwas_trait GWAS trait (e.g. IBD, CRC, etc.)
#' @param outdir Directory where outputs are saved
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A list containing the fitted topic model and cell scores.
#'
#' @examples
#' result <- scads(count_matrix, GWAS_ss, LD_reference, nTopics = 15)
#'
#' @export
#' 
scads <- function(count_matrix, nTopics = 10, n_s=1000, n_c=8, 
                  sumstats_dir, gwas_nsamps, gwas_trait, outdir, ...){
  
  # step0: create directories 
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # step1: run fastTopics
  cat("\n1. Run FastTopics\n")
  cat("\nStart time: ")
  print(Sys.time())
  cat("\nCount matrix dimensions", dim(count_matrix))
  out1 <- run_fastTopics(count_matrix, nTopics, n_s, n_c)
  saveRDS(out1, file.path(outdir, "run_fastTopics_res.rds"))
  # cat("\n Dimensions of topics_prob", dim(out1$Pmat))
  
  # step2: get topic annotations (bed files)
  cat("\n2. Get topic annotations\n")
  cat("\nStart time: ")
  print(Sys.time())
  beddir_list <- get_annot_beds(topics_res = out1, output_dir = outdir)

  # step3: run LDSC
  
  # Detect the number of physical cores
  num_cores <- parallel::detectCores(logical = FALSE)
  # Determine the number of tasks
  num_tasks <- length(beddir_list)
  # Set mc.cores to the minimum of num_cores and num_tasks
  mc_cores_to_use <- min(num_cores, num_tasks, 5) # 5 - set to 1 to see if this better for larger # of peaks
  
  cat("\nUsing", mc_cores_to_use, "cores for", num_tasks, "tasks.\n")
  cat("\nStart time: ")
  print(Sys.time())
  parallel::mclapply(1:length(beddir_list), function(i){
    cat("\nStart time for topic ", i)
    print(Sys.time())
    run_ldsc(polyfun_path = "/dartfs/rc/lab/S/Szhao/liyang/polyfun",
             sumstats_path = sumstats_dir,
             n = gwas_nsamps,
             trait = gwas_trait,
             onekg_path = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/1000G_EUR_Phase3_plink/1000G.EUR.QC",
             bed_dir = beddir_list[[i]],
             baseline_dir = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/Gazel_LD_baseline",
             weights_dir = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/weights",
             out_dir = beddir_list[[i]])
    cat("\nStop time for topic ", i)
    print(Sys.time())
  }, mc.cores = mc_cores_to_use)
  cat("\nStop time: ")
  print(Sys.time())
  
  # enrichment <- run_ldsc(polyfun_path = "/dartfs/rc/lab/S/Szhao/liyang/polyfun",
  #                        sumstats_path = sumstats_dir,
  #                        n = gwas_nsamps,
  #                        trait = gwas_trait,
  #                        onekg_path = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/1000G_EUR_Phase3_plink/1000G.EUR.QC",
  #                        bed_dir = beddir,
  #                        baseline_dir = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/Gazel_LD_baseline",
  #                        weights_dir = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/weights",
  #                        out_dir = outdir)

  # step4: calc the cell score
  cat("\n4. Calculate cell scores\n")
  cat("\nStart time: ")
  print(Sys.time())
  cs_res <- get_cs(topic_res = out1, 
               ldsc_res_dir = outdir, 
               trait = gwas_trait,
               nTopics = nTopics)
  print(summary(cs_res$cs))

  return(list(cs_res = cs_res, out1 = out1))
  
  cat("\nEnd time: ")
  print(Sys.time())
  
}





  
  