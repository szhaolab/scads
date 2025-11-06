#' Run S-LDSC Pipeline (hg38) with a Custom BED + Baseline v2.2
#'
#' This function runs the S-LDSC pipeline using PolyFun scripts for:
#'   1) Munging summary statistics ("munge_polyfun_sumstats.py"),
#'   2) Generating annotation files from a user-supplied BED (hg38) via "make_annot.py",
#'   3) Computing LD scores with "ldsc.py --l2",
#'   4) Performing final heritability estimation ("ldsc.py --h2")
#'      by combining the custom annotation with baseline v2.2 annotations.
#'
#' @param polyfun_path A character string specifying the directory containing the PolyFun scripts 
#'   (e.g.,\code{ldsc.py}, \code{munge_polyfun_sumstats.py}).
#' @param ldsc_path A character string specifying the directory containing the ldsc script ( (e.g., \code{make_annot.py})
#' @param sumstats_path A character string specifying the directory with the cleaned summary 
#'   statistics file \code{<trait>_sumstats.txt.gz}.
#' @param n Integer. The number of samples in the GWAS.
#' @param trait A character string representing the GWAS trait (e.g., \code{"SIM"}).
#' @param onekg_path A character string specifying the prefix to 1000 Genomes Plink files in hg38,
#'   e.g. \code{"/path/1000G.EUR.hg38."}. Each chromosome file is then 
#'   \code{1000G.EUR.hg38.1.bed/.bim/.fam}, etc.
#' @param bed_dir A character string specifying the directory containing the user-supplied BED file(s)
#'   in hg38. Only the first \code{.bed} found is used to create the annotation.
#' @param baseline_dir A character string specifying the prefix to baseline v2.2 annotations for LDSC,
#'   e.g. \code{"/path/baselineLD_v2.2/baselineLD."}.
#' @param frqfile_pref A character string specifying the prefix to allele frequency files in hg38,
#'   used by LDSC (e.g. \code{"/path/1000G.EUR.hg38."}).
#' @param hm3_snps A character string specifying the path to the HapMap3 SNP file 
#'   (e.g. \code{"hm3_no_MHC.list.txt"}).
#' @param out_dir A character string specifying the output directory for LDSC results. 
#'   The function will create:
#'   \itemize{
#'     \item \code{"annotations/<trait>"} for annotation files
#'     \item \code{"results"} for final LDSC logs/results
#'   }
#'
#' @details
#' **Step-by-step**:
#' \enumerate{
#'   \item Munge summary stats \code{<trait>_sumstats.txt.gz} into a parquet with 
#'         \code{munge_polyfun_sumstats.py}.
#'   \item For each chromosome \code{1..22}:
#'         \itemize{
#'           \item Create annotation (\code{make_annot.py}) from the user .bed plus the .bim for that chromosome.
#'           \item Compute LD scores (\code{ldsc.py --l2}) for the annotation.
#'         }
#'   \item Finally, run \code{ldsc.py --h2} referencing both the newly created annotation 
#'         (\code{out_dir/annotations/<trait>/<trait>.<chr>}) and \code{baseline_dir}, 
#'         providing \code{frqfile_pref} for the allele frequency files and the \code{hm3_snps} 
#'         for SNP filtering.
#' }
#'
#' @return No object is returned. All intermediate files (\code{.annot.gz}, \code{.ldscore.gz}, etc.)
#'   are written to \code{out_dir/annotations/<trait>}, and the final LDSC results (\code{.log} / 
#'   \code{.results}) to \code{out_dir/results}.
#'
#' @examples
#' \dontrun{
#' run_sldsc(
#'   polyfun_path = "/path/polyfun",
#'   sumstats_path= "/path/sumstats",
#'   n            = 60000,
#'   trait        = "SIM",
#'   onekg_path   = "/path/1000G.EUR.hg38.",
#'   bed_dir      = "/path/my_bedfiles",
#'   baseline_dir = "/path/baselineLD_v2.2/baselineLD.",
#'   frqfile_pref = "/path/1000G.EUR.hg38.",
#'   hm3_snps     = "/path/hm3_no_MHC.list.txt",
#'   out_dir      = "/path/k1_annotations"
#' )
#' }
#'
#' @export
run_sldsc <- function(polyfun_path,
                      ldsc_path,
                      sumstats_path,
                      n,
                      trait,
                      onekg_path,
                      bed_dir,
                      baseline_dir,
                      frqfile_pref,
                      hm3_snps,
                      weights_pref,
                      out_dir) {
  
  # 1) Create subfolders
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  ann_dir <- file.path(out_dir, "annotations", trait)  # e.g. "out_dir/annotations/SIM"
  if (!dir.exists(ann_dir)) {
    dir.create(ann_dir, recursive = TRUE)
  }
  res_dir <- file.path(out_dir, "results")
  if (!dir.exists(res_dir)) {
    dir.create(res_dir, recursive = TRUE)
  }
  
  # 2) Munge sumstats
  sumstats_gz <- list.files(
    path = sumstats_path,
    pattern = paste0("^", trait, ".*sumstats.*\\.txt\\.gz$"),
    full.names = TRUE
  )
  
  # Check number of matches
  if (length(sumstats_gz) == 0) {
    stop(paste0("No matching sumstats file found for trait: ", trait))
  } else if (length(sumstats_gz) > 1) {
    stop(paste0("Multiple matching files found for trait: ", trait, 
                "\nFiles:\n", paste(sumstats_gz, collapse = "\n")))
  } else {
    message(paste0("Found summary stats file: ", basename(sumstats_gz)))
  }
  
  munged_out  <- file.path(out_dir, sprintf("%s_munged_sumstats.parquet", trait))
  cmd_munge <- paste(
    "python3", file.path(polyfun_path, "munge_polyfun_sumstats.py"),
    "--sumstats", sumstats_gz,
    "--n", n,
    "--out", munged_out
  )
  message("\n[run_sldsc] Step 1: Munging sumstats:\n", cmd_munge)
  system(cmd_munge)
  
  # 3) Build annotation + compute LD scores for each chromosome
  bed_files <- list.files(unlist(bed_dir), pattern = "\\.bed$", full.names = TRUE)
  if (length(bed_files) < 1) {
    stop("No BED file (.bed) found in: ", bed_dir)
  }
  message("[run_sldsc] Found BED file(s): ", paste(bed_files, collapse = ", "))
  user_bed <- bed_files[1]
  
  for (chr in 10) {
    # (a) make_annot
    out_annot <- file.path(ann_dir, sprintf("%s.%d.annot.gz", trait, chr))
    bim_file  <- sprintf("%s%d.bim", onekg_path, chr)
    
    cmd_annot <- paste(
      "python3", file.path(ldsc_path, "make_annot.py"),
      "--bed-file", user_bed,
      "--bimfile", bim_file,
      "--annot-file", out_annot
    )
    message("\n[run_sldsc] Step 2a: make_annot:\n", cmd_annot)
    system(cmd_annot)
    
    # (b) ldsc.py --l2
    out_prefix <- file.path(ann_dir, sprintf("%s.%d", trait, chr))
    bfile_chr  <- sprintf("%s%d", onekg_path, chr)
    
    cmd_l2 <- paste(
      "python3", file.path(polyfun_path, "ldsc.py"),
      "--l2",
      "--bfile", bfile_chr,
      "--print-snps", hm3_snps,
      "--ld-wind-cm 1",
      "--annot", out_annot,
      "--thin-annot",
      "--out", out_prefix
    )
    message("\n[run_sldsc] Step 2b: ldsc --l2:\n", cmd_l2)
    system(cmd_l2)
    message("\nFinish step 2b: ", cmd_l2)
  }
  
  message("\nBegin step 3")
  # 4) Final S-LDSC referencing new annotation + baseline v2.2
  final_out <- file.path(res_dir, trait)
  message("\n Creating new file: ", final_out)
  
  # Count how many annotation files exist
  annot_files <- list.files(ann_dir, pattern = "\\.annot\\.gz$", full.names = TRUE)
  n_annot <- length(annot_files)
  message("\n Number of chromosomes with annotation: ", n_annot)
  
  if (n_annot == 1) {
    # for one chromosome 
    cmd_h2 <- paste(
      "python3", file.path(polyfun_path, "ldsc.py"),
      "--h2", munged_out,
      "--ref-ld", paste0(
        file.path(ann_dir, paste0(trait, ".", chr)), ",",
        paste0(baseline_dir,chr)
      ),
      "--frqfile", paste0(frqfile_pref, chr),
      "--w-ld", paste0(weights_pref, chr),
      "--overlap-annot",
      "--print-coefficients",
      "--print-delete-vals",
      "--out", paste0(final_out, "_chr", chr)
    )
    
  } else if (n_annot == 22) {
    # for all chromosomes (1-22)
    cmd_h2 <- paste(
      "python3", file.path(polyfun_path, "ldsc.py"),
      "--h2", munged_out,
      "--ref-ld-chr", paste0(file.path(ann_dir, paste0(trait, ".")), ",", baseline_dir),
      "--frqfile-chr", frqfile_pref,
      "--w-ld-chr", weights_pref,
      "--overlap-annot",
      "--print-coefficients",
      "--print-delete-vals",
      "--out", final_out
    )
    
  } else {
    # Invalid number of annotation files
    stop("Please check number of CHR — S-LDSC only takes in either one chromosome or 22 chromosomes.")
  }
  
  message("\n[run_sldsc] Step 3: ldsc.py --h2:\n", cmd_h2)
  system(cmd_h2)
  
  message("\n[run_sldsc] Completed S-LDSC steps. Output in: ", res_dir)
}

# ref-ld-chr => "out_dir/annotations/trait/trait.," + baseline_dir
# We also pick up the weights from the parent of the parent of baseline_dir (since baseline_dir might be .../baselineLD_v2.2/baselineLD.)
# Adjust if your weights directory is stored differently.
# weights_pref <- file.path(dirname(dirname(baseline_dir)), "weights/weights.hm3_noMHC.") #hg38
# weights_pref <- file.path(dirname(dirname(baseline_dir)), "1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.") #hg19

