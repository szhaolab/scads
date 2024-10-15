#' Run LDSC Pipeline
#'
#'
#' This function runs the LDSC pipeline using PolyFun for enrichment analysis.
#' 
#' @import reticulate
#'
#' @param polyfun_path Path to the PolyFun package.
#' @param sumstats_path Path to the cleaned summary statistics file.
#' @param n Number of samples in the GWAS.
#' @param mss_output Munged sumstats file.
#' @param ldblocks_path Path to file with LD blocks.
#' @param sig_loci_output Signficant loci output file.
#' 
#' @return enrichment
#' @export
#' 
run_ldsc <- function(polyfun_path, sumstats_path, n, mss_output, 
                     ldblocks_path, sig_loci_output = NULL) {
  
  # Use the virtual environment
  reticulate::use_virtualenv("scads_env", required = TRUE)
  
  # Check if required packages are available
  required_packages <- c("pandas", "numpy", "scipy", "logging", "tqdm", "pandas_plink")
  missing_packages <- required_packages[!sapply(required_packages, reticulate::py_module_available)]
  if (length(missing_packages) > 0) {
    stop(paste("The following required Python packages are missing in the 'scads_env' environment:", 
               paste(missing_packages, collapse = ", "), 
               "\nPlease run setup_python_env() to install them."))
  }
  
  # 1. Munge sumstats
  cmd1 <- paste(
    "python3", 
    file.path(polyfun_path, "munge_polyfun_sumstats.py"),
    "--sumstats", sumstats_path,
    "--n", n,
    "--out", mss_output
  ) 
  system(cmd1)
  
  # Get significant loci - for fine-mapping - works but not necessary for now
  # sig_loci <- get_signif_loci(sumstats_path, 
  #                             ldblocks_path, 
  #                             pval_threshold = 5e-8, 
  #                             max_loci = 150, 
  #                             outfile = sig_loci_output) 
  
  # 2. Create annotations 
  
  
}

