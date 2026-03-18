#' scads: A package for using scATAC-seq data and GWAS sumstats to calculate disease cell score
#'
#' @param count_matrix A peak-by-cell count matrix.
#' @param nTopics The number of topics (clusters) to fit. Default is 10.
#' @param n_s Number of simulations for fastTopics-DE.
#' @param n_c Number of cores to use for parallelization.
#' @param baseline_method Method for calculating baseline f_kj0 (Default is gc; other options: constant, average, estimate)
#' @param bl_celltype Use only if baseline_method == "estimate", otherwise NULL
#' @param bl_celltype_peak_file Use only if baseline_method == "estimate", otherwise NULL
#' @param sumstats_dir Directory containing GWAS sumstats (file format: TRAIT.sumstats.txt.gz).
#' @param gwas_nsamps Sample size of GWAS.
#' @param gwas_trait GWAS trait name (e.g., "SIM", "IBD").
#' @param outdir Directory where outputs are saved.
#' @param polyfun_code_dir Directory to access polyFUN (LDSC) code scripts; must be python3-compatible.
#' @param ldsc_code_dir Directory to access LDSC code scripts.
#' @param onekg_path Prefix to 1000G plink files in hg19 (/LDSCORE/zenodo/1000G_EUR_Phase3_plink/1000G.EUR.QC.). (use hg38 version if input data is hg38).
#' @param baseline_dir Prefix to 1000G baseline v2.2 annotation in hg19 (e.g. "/LDSCORE/zenodo/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD."). (use hg38 version if input data is hg38).
#' @param frqfile_pref Prefix to the 1000G  .frq or .afreq files in hg19 (e.g. "/LDSCORE/zenodo/1000G_Phase3_frq/1000G.EUR.QC."). (use hg38 version if input data is hg38).
#' @param hm3_snps Path to the HapMap3 no‐MHC SNP list in hg19 (e.g. "/LDSCORE/zenodo/hm3_no_MHC.list.txt").(use hg38 version if input data is hg38).
#' @param weights_pref Path to the 1000G weights in hg19 (e.g. "/LDSCORE/zenodo/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.").(use hg38 version if input data is hg38).
#' @param continuous_topic_annot Logical. If \code{TRUE}, create \emph{continuous} topic annotations 
#'   (via \code{get_continuous_annot_beds()}) and run \code{run_sldsc_cont()}. 
#'   If \code{FALSE}, use the \emph{binary} annotation approach (\code{get_annot_beds()} and \code{run_sldsc()}). 
#'   Default is \code{FALSE}.
#' @param chrs Number of chromosomes for running S-LDSC (Default is chr1-22)
#' @param genome Genome build (Default is hg19)
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A list containing the fitted topic model (\code{out1}) and the cell scores (\code{cs_res}) results.
#'
#' @examples
#' result <- scads(count_matrix, nTopics = 15, sumstats_dir="/some/sumstats", gwas_nsamps=60000,
#'                 gwas_trait="SIM", outdir="results/", 
#'                 polyfun_code_dir = "/some/path/polyfun/", ldsc_code_dir = "/some/path/ldsc/",
#'                 onekg_path="/some/path/1000G.EUR.QC.",
#'                 baseline_dir="/some/path/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.",
#'                 frqfile_pref="/some/path/1000G_Phase3_frq/1000G.EUR.QC.",
#'                 hm3_snps="/some/path/hm3_no_MHC.list.txt", 
#'                 weights_pref="/some/path/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.")
#'
#' @export
scads <- function(count_matrix, 
                  nTopics = 10, 
                  n_s = 1000, 
                  n_c = 8, 
                  baseline_method = "gc", # c("constant", "average", "estimate"),
                  bl_celltype = NULL, 
                  bl_celltype_peak_file = NULL,
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
                  weights_pref,
                  continuous_topic_annot = FALSE,
                  chrs = 1:22,
                  genome = "hg19",
                  ...) {
  
  # 0) create directories 
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # 1) run fastTopics
  cat("\n1. Run FastTopics\n")
  cat("\nStart time: "); print(Sys.time())
  cat("\nCount matrix dimensions:", dim(count_matrix), "\n")
  
  cat("\nTransposing matrix for fastTopics")
  count_matrix_t <- Matrix::t(count_matrix)
  print(dim(count_matrix))
  cat("\nTransposed count matrix dimensions", dim(count_matrix_t), "\n")
  rm(count_matrix)
  
  # Ensure that count_matrix has cells as rows and peaks as columns
  out1 <- run_fastTopics(count_matrix_t, nTopics, n_s, n_c,
                         baseline_method, bl_celltype, bl_celltype_peak_file,
                         outdir = outdir, genome = genome)
  saveRDS(out1, file.path(outdir, "run_fastTopics_res.rds"))
  # out1 <- readRDS(file.path(outdir, "run_fastTopics_res.rds"))
  
  # 2) get topic annotations (bed files)
  cat("\n2. Get topic annotations\n")
  cat("\nStart time: "); print(Sys.time())

  if (!continuous_topic_annot) {
    # Binary (peak-based) annotation approach
    beddir_list <- get_annot_beds(topics_res = out1, output_dir = outdir)
  } else {
    # Continuous annotation approach
    beddir_list <- get_continuous_annot_beds(topics_res = out1, output_dir = outdir)
  }

  # 3) run LDSC
  cat("\n3. Run LDSC\n")
  cat("\nStart time: "); print(Sys.time())

  num_cores <- parallel::detectCores(logical = FALSE)
  num_tasks <- length(beddir_list)
  mc_cores_to_use <- min(num_cores, num_tasks, 5)

  cat("\nUsing", mc_cores_to_use, "cores for", num_tasks, "tasks.\n")

  # sldsc_results <- 
  parallel::mclapply(seq_along(beddir_list), function(i) {

    #tryCatch({
    cat("\nStart time for topic", i, ":", Sys.time(), "\n")
    # If using continuous annotations, call run_sldsc_cont; else run_sldsc
    if (!continuous_topic_annot) {
      run_sldsc(
        chrs          = chrs,
        polyfun_path  = polyfun_code_dir,
        ldsc_path     = ldsc_code_dir,
        sumstats_path = sumstats_dir,
        n             = gwas_nsamps,
        trait         = gwas_trait,
        onekg_path    = onekg_path,
        bed_dir       = beddir_list[[i]],
        baseline_dir  = baseline_dir,
        frqfile_pref  = frqfile_pref,
        hm3_snps      = hm3_snps,
        weights_pref  = weights_pref,
        out_dir       = beddir_list[[i]]
      )
    } else {
      run_sldsc_cont(
        polyfun_path  = polyfun_code_dir,
        ldsc_path     = ldsc_code_dir,
        sumstats_path = sumstats_dir,
        n             = gwas_nsamps,
        trait         = gwas_trait,
        onekg_path    = onekg_path,
        bed_dir       = beddir_list[[i]],
        baseline_dir  = baseline_dir,
        frqfile_pref  = frqfile_pref,
        hm3_snps      = hm3_snps,
        out_dir       = beddir_list[[i]]
      )
    }
    # }, error = function(e) e)

    cat("\nStop time for topic", i, ":", Sys.time(), "\n")
  }, mc.cores = mc_cores_to_use)

  # print(sldsc_results[sapply(sldsc_results, inherits, "error")])

  cat("\nStop time:", Sys.time(), "\n")

  # 4) calc the cell score
  cat("\n4. Calculate cell scores\n")
  cat("\nStart time:", Sys.time(), "\n")
  cs_res <- get_cs(topic_res  = out1,
                   ldsc_res_dir = outdir,
                   trait      = gwas_trait,
                   nTopics    = nTopics)
  cat("\nCell score summary: ", summary(cs), "\n")

  cat("\nEnd time:", Sys.time(), "\n")

  return(list(cs_res = cs_res, out1 = out1))
}
