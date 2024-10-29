#' Run LDSC Pipeline
#'
#'
#' This function runs the LDSC pipeline using PolyFun for enrichment analysis.
#' 
#' @importFrom reticulate use_virtualenv
#'
#' @param polyfun_path A character string specifying the directory to the PolyFun package.
#' @param sumstats_path A character string specifying the directory with the cleaned summary statistics file.
#' @param n Number of samples in the GWAS.
#' @param trait GWAS trait (e.g. SIM for simulated)
#' @param rcode_path A character string specifying the directory to Rscripts for running LDSC
#' @param onekg_path A character string specifying the directory to 1000G genotypes data
#' @param bed_dir A character string specifying the directory for the BED files.
#' @param baseline_dir A character string specifying the directory to baseline annotations for LDSC
#' @param weights_dir A character string specifying the directory to weights provided for LDSC
#' @param out_dir A character string specifying the output directory for LDSC enrichment results.
#' 
#' @return None. LDSC results are written to the specified output directory.
#' @export
#' 
run_ldsc <- function(polyfun_path, 
                     sumstats_path, 
                     n, 
                     trait,
                     rcode_path, 
                     onekg_path, 
                     bed_dir, 
                     baseline_dir, 
                     weights_dir, 
                     out_dir) { 
  
  # # Use the virtual environment
  # reticulate::use_virtualenv("scads_env", required = TRUE)
  # 
  # # Check if required packages are available
  # required_packages <- c("pandas", "numpy", "scipy", "logging", "tqdm", "pandas_plink")
  # missing_packages <- required_packages[!sapply(required_packages, reticulate::py_module_available)]
  # if (length(missing_packages) > 0) {
  #   stop(paste("The following required Python packages are missing in the 'scads_env' environment:", 
  #              paste(missing_packages, collapse = ", "), 
  #              "\nPlease run setup_python_env() to install them."))
  # }
  
  # Ensure the output directory exists
  if (!dir.exists(paste0(out_dir, "/annotations"))) {
    dir.create(paste0(out_dir, "/annotations"), recursive = TRUE)
  }
  if (!dir.exists(paste0(out_dir, "/results"))) {
    dir.create(paste0(out_dir, "/results"), recursive = TRUE)
  }
  if (!dir.exists(paste0(out_dir, sprintf("/annotations/%s", trait)))) {
    dir.create(paste0(out_dir, sprintf("/annotations/%s", trait)), recursive = TRUE)
  }
  
  # 1. Munge sumstats
  cmd1 <- paste(
    "python3", 
    file.path(polyfun_path, "munge_polyfun_sumstats.py"),
    "--sumstats", file.path(sumstats_path, sprintf('%s_sumstats.txt.gz', trait)),
    "--n", n,
    "--out", file.path(out_dir, sprintf('%s_munged_sumstats.parquet', trait))
  ) 
  print(cmd1)
  #system(cmd1)
  
  
  # Get significant loci - for fine-mapping - works but not necessary for now
  # sig_loci <- get_signif_loci(sumstats_path, 
  #                             ldblocks_path, 
  #                             pval_threshold = 5e-8, 
  #                             max_loci = 150, 
  #                             outfile = sig_loci_output) 
  
  # 2. Create annotations 
  for (chr in 1:22){
    cmd2 <- paste(
      "Rscript", 
      file.path(rcode_path, "create_annotations_LY.R"),
      paste0(onekg_path, sprintf(".%s.bim",chr)), 
      file.path(bed_dir),
      file.path(baseline_dir, sprintf("baseline_MAF_LD.%s.annot.gz", chr)), 
      file.path(out_dir, sprintf('annotations/%s/%s.%s.annot.gz', trait, trait, chr)),
      file.path(out_dir, sprintf('annotations/%s/%s.%s.l2.M', trait, trait, chr))
    )
    print(cmd2)
    #system(cmd2)
    
  }
  
  # 3. Compute LD score 
  for (chr in 1:22){
    cmd3 <- paste(
      "python3", 
      file.path(polyfun_path, "compute_ldscores.py"),
      "--bfile", paste0(onekg_path, sprintf(".%s",chr)),
      "--annot", file.path(out_dir, sprintf('annotations/%s/%s.%s.annot.gz', trait, trait, chr)),
      "--out", file.path(out_dir, sprintf('annotations/%s/%s.%s.l2.ldscore.parquet', trait, trait, chr))
    ) 
    print(cmd3)
    #system(cmd3)
    
  }
  
  # 4. Calculate enrichment 
  cmd4 <- paste(
    "python3", 
    file.path(polyfun_path, "ldsc.py"),
    "--h2", file.path(out_dir, sprintf('%s_munged_sumstats.parquet', trait)),
    "--ref-ld-chr", file.path(out_dir, sprintf('annotations/%s/%s.', trait, trait)),
    "--w-ld-chr", file.path(weights_dir, 'weights.'),
    "--not-M-5-50",
    "--out", file.path(out_dir, sprintf('results/%s_enrichment', trait)),
    "--overlap-annot && sed -i 's/_0//g'", file.path(out_dir, sprintf('results/%s_enrichment.results', trait))
  ) 
  print(cmd4)
  #system(cmd4)
  
  
}