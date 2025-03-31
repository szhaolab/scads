#' scads: A function for using scATAC-seq data and GWAS sumstats to calculate disease cell score
#'
#' @param count_matrix A peak-by-cell count matrix.
#' @param nTopics The number of topics (clusters) to fit. Default is 10.
#' @param n_s Number of simulations for fastTopics-DE.
#' @param n_c Number of cores to use for parallel.
#' @param sumstats_dir Directory with GWAS sumstats.
#' @param gwas_nsamps Sample size of GWAS.
#' @param gwas_trait GWAS trait name (e.g., "SIM").
#' @param outdir Directory where outputs are saved.
#' @param polyfun_code_dir Directory to access polyFUN (LDSC) code scripts; polyFUN has python3 compatible code
#' @param ldsc_code_dir Directory to access LDSC code scripts
#' @param onekg_path Prefix to 1000G plink files in hg38 (e.g., "/path/1000G.EUR.hg38.").
#' @param baseline_dir Prefix to baseline v2.2 annotation (e.g. "/path/baselineLD_v2.2/baselineLD.").
#' @param frqfile_pref Prefix to the .frq or .afreq files (e.g. "/path/1000G.EUR.hg38.").
#' @param hm3_snps Path to the HapMap3 no‐MHC SNP list (e.g. "hm3_no_MHC.list.txt").
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A list containing the fitted topic model (`out1`) and the cell scores (`cs_res`).
#' @examples
#' result <- scads(count_matrix, nTopics = 15, sumstats_dir="/some/sumstats", gwas_nsamps=60000,
#'                 gwas_trait="SIM", outdir="results/", 
#'                 onekg_path="/some/path/1000G.EUR.hg38.",
#'                 baseline_dir="/some/path/baselineLD_v2.2/baselineLD.",
#'                 frqfile_pref="/some/path/1000G.EUR.hg38.",
#'                 hm3_snps="/some/path/hm3_no_MHC.list.txt")
#'
#' @export
scads <- function(count_matrix, 
                  nTopics = 10, 
                  n_s = 1000, 
                  n_c = 8, 
                  sumstats_dir, 
                  gwas_nsamps, 
                  gwas_trait, 
                  outdir,
                  polyfun_code_dir,
                  ldsc_code_dir,
                  onekg_path,
                  baseline_dir,
                  frqfile_pref,
                  hm3_snps,
                  ...){
  
  # 0) create directories 
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # 1) run fastTopics
  cat("\n1. Run FastTopics\n")
  cat("\nStart time: "); print(Sys.time())
  cat("\nCount matrix dimensions:", dim(count_matrix), "\n")
  out1 <- run_fastTopics(count_matrix, nTopics, n_s, n_c)
  saveRDS(out1, file.path(outdir, "run_fastTopics_res.rds"))
  # out1 <- readRDS(file.path(outdir, "run_fastTopics_res.rds"))
  
  # 2) get topic annotations (bed files)
  cat("\n2. Get topic annotations\n")
  cat("\nStart time: "); print(Sys.time())
  beddir_list <- get_annot_beds(topics_res = out1, output_dir = outdir)
  
  # 3) run LDSC
  cat("\n3. Run LDSC \n")
  cat("\nStart time: "); print(Sys.time())
  
  num_cores <- parallel::detectCores(logical = FALSE)
  num_tasks <- length(beddir_list)
  mc_cores_to_use <- min(num_cores, num_tasks, 5)
  
  cat("\nUsing", mc_cores_to_use, "cores for", num_tasks, "tasks.\n")
  
  parallel::mclapply(seq_along(beddir_list), function(i){
    cat("\nStart time for topic", i, ":", Sys.time(), "\n")
    
    # For each topic, call the updated run_ldsc
    run_sldsc(polyfun_path  = polyfun_code_dir,
              ldsc_path = ldsc_code_dir,
             sumstats_path = sumstats_dir,
             n             = gwas_nsamps,
             trait         = gwas_trait,
             onekg_path    = onekg_path,
             bed_dir       = beddir_list[[i]],
             baseline_dir  = baseline_dir,
             frqfile_pref  = frqfile_pref,
             hm3_snps      = hm3_snps,
             out_dir       = beddir_list[[i]])
    
    cat("\nStop time for topic", i, ":", Sys.time(), "\n")
  }, mc.cores = mc_cores_to_use)
  
  cat("\nStop time:", Sys.time(), "\n")
  
  # 4) calc the cell score
  cat("\n4. Calculate cell scores\n")
  cat("\nStart time:", Sys.time(), "\n")
  cs_res <- get_cs(topic_res  = out1,
                   ldsc_res_dir = outdir,
                   trait      = gwas_trait,
                   nTopics    = nTopics)
  print(summary(cs_res$cs))
  
  cat("\nEnd time:", Sys.time(), "\n")
  
  return(list(cs_res = cs_res, out1 = out1))
}


# scads <- function(count_matrix, nTopics = 10, n_s=1000, n_c=8, 
#                   sumstats_dir, gwas_nsamps, gwas_trait, outdir, ...){
#   
#   # step0: create directories 
#   if (!dir.exists(outdir)) {
#     dir.create(outdir, recursive = TRUE)
#   }
#   
#   # step1: run fastTopics
#   cat("\n1. Run FastTopics\n")
#   cat("\nStart time: ")
#   print(Sys.time())
#   cat("\nCount matrix dimensions", dim(count_matrix))
#   out1 <- run_fastTopics(count_matrix, nTopics, n_s, n_c)
#   saveRDS(out1, file.path(outdir, "run_fastTopics_res.rds"))
#   # cat("\n Dimensions of topics_prob", dim(out1$Pmat))
#   
#   # step2: get topic annotations (bed files)
#   cat("\n2. Get topic annotations\n")
#   cat("\nStart time: ")
#   print(Sys.time())
#   beddir_list <- get_annot_beds(topics_res = out1, output_dir = outdir)
# 
#   # step3: run LDSC
#   cat("\n3. Run LDSC \n")
#   cat("\nStart time: ")
#   print(Sys.time())
#   # Detect the number of physical cores
#   num_cores <- parallel::detectCores(logical = FALSE)
#   # Determine the number of tasks
#   num_tasks <- length(beddir_list)
#   # Set mc.cores to the minimum of num_cores and num_tasks
#   mc_cores_to_use <- min(num_cores, num_tasks, 5) # 5 - set to 1 to see if this better for larger # of peaks
#   
#   cat("\nUsing", mc_cores_to_use, "cores for", num_tasks, "tasks.\n")
#   cat("\nStart time: ")
#   print(Sys.time())
#   parallel::mclapply(1:length(beddir_list), function(i){
#     cat("\nStart time for topic ", i)
#     print(Sys.time())
#     run_ldsc(polyfun_path = "/dartfs/rc/lab/S/Szhao/liyang/polyfun",
#              sumstats_path = sumstats_dir,
#              n = gwas_nsamps,
#              trait = gwas_trait,
#              onekg_path = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/1000G_EUR_Phase3_plink/1000G.EUR.QC",
#              bed_dir = beddir_list[[i]],
#              baseline_dir = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/Gazel_LD_baseline",
#              # baseline_dir = "/dartfs/rc/lab/S/Szhao/liyang/enrichment/LDSCORE/baselineLF2.2.UKB",
#              weights_dir = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/weights",
#              # weights_dir = "/dartfs/rc/lab/S/Szhao/liyang/enrichment/LDSCORE/baselineLF2.2.UKB",
#              out_dir = beddir_list[[i]])
#     cat("\nStop time for topic ", i)
#     print(Sys.time())
#   }, mc.cores = mc_cores_to_use)
#   cat("\nStop time: ")
#   print(Sys.time())
#   
#   # enrichment <- run_ldsc(polyfun_path = "/dartfs/rc/lab/S/Szhao/liyang/polyfun",
#   #                        sumstats_path = sumstats_dir,
#   #                        n = gwas_nsamps,
#   #                        trait = gwas_trait,
#   #                        onekg_path = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/1000G_EUR_Phase3_plink/1000G.EUR.QC",
#   #                        bed_dir = beddir,
#   #                        baseline_dir = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/Gazel_LD_baseline",
#   #                        weights_dir = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/weights",
#   #                        out_dir = outdir)
# 
#   # step4: calc the cell score
#   cat("\n4. Calculate cell scores\n")
#   cat("\nStart time: ")
#   print(Sys.time())
#   cs_res <- get_cs(topic_res = out1, 
#                ldsc_res_dir = outdir, 
#                trait = gwas_trait,
#                nTopics = nTopics)
#   print(summary(cs_res$cs))
# 
#   return(list(cs_res = cs_res, out1 = out1))
#   
#   cat("\nEnd time: ")
#   print(Sys.time())
#   
# }
# 




  
  