#' scads: A function for using scATAC-seq data and GWAS sumstats to calculate disease cell score
#'
#' @param count_matrix A peak-by-cell count matrix.
#' @param nTopics The number of topics (clusters) to fit. Default is 10.
#' @param bed_dir 
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A list containing the fitted topic model and cell scores.
#'
#' @examples
#' result <- scads(count_matrix, GWAS_ss, LD_reference, nTopics = 15)
#'
#' @export
#' 
scads <- function(count_matrix, nTopics = 10, bed_dir, sumstats_dir, gwas_nsamps, gwas_trait, outdir, ...){
  
  # step0: create directories 
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # step1: run fastTopics
  cat("\n1. Run FastTopics\n")
  out1 <- run_fastTopics(count_matrix, nTopics)
  
  # step2: get topic annotations (bed files)
  cat("\n2. Get topic annotations\n")
  beddir <- get_annot_beds(topics_res = out1, output_dir = bed_dir)

  # step3: run LDSC
  cat("\n3. Perform LDSC\n")
  enrichment <- run_ldsc(polyfun_path = "/dartfs/rc/lab/S/Szhao/liyang/polyfun",
                         sumstats_path = sumstats_dir,
                         n = gwas_nsamps,
                         trait = gwas_trait,
                         rcode_path = "/dartfs/rc/lab/S/Szhao/liyang/enrichment/ARID1A_proj",
                         onekg_path = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/1000G_EUR_Phase3_plink/1000G.EUR.QC",
                         bed_dir = beddir,
                         baseline_dir = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/Gazel_LD_baseline",
                         weights_dir = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/weights",
                         out_dir = outdir)

  # step4: calc the cell score
  cat("\n4. Calculate cell scores\n")
  cs <- get_cs(topic_res = out1, 
               ldsc_res = outdir)

  return(list(cs = cs, out1 = out1))
  
}





  
  