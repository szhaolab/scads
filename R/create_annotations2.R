#' Create SNP Annotations Based on BED Files and Baseline Annotations
#'
#' This function annotates SNPs based on overlap with BED files and merges them with baseline annotations.
#' It outputs the combined annotations and computes the number of SNPs in each annotation.
#'
#' @param bim_file A data frame or character string specifying the path to the BIM file. The BIM file should have columns:
#'   - \code{SNP}: SNP identifier
#'   - \code{CHR}: Chromosome number
#'   - \code{BP}: Base pair position
#'   - \code{A1}: Allele 1
#'   - \code{A2}: Allele 2
#' @param bed_files A character string specifying the directory containing BED files, or a vector of BED file paths. Default is \code{NULL}.
#' @param baseline_file A data frame or character string specifying the path to the baseline annotations file.
#' @param outfile Optional. A character string specifying the path to save the combined annotations. Default is \code{NULL}.
#' @param Mfile Optional. A character string specifying the path to save the number of SNPs per annotation. Default is \code{NULL}.
#' @return A data frame containing the combined annotations.
#' @export
#' @importFrom GenomicRanges GRanges 
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom vroom vroom vroom_write
#' @importFrom plyranges mutate
#' @importFrom dplyr mutate inner_join
#' @importFrom rtracklayer import
#' @importFrom data.table fread
#' @importFrom stringr str_detect str_c str_replace
#' @importFrom tools file_path_sans_ext
create_annotations2 <- function(bim_file, bed_files = NULL, baseline_file, outfile = NULL, Mfile = NULL) {
  
  # Load necessary libraries
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(IRanges)
    library(vroom)
    library(plyranges)
    library(dplyr)
    library(rtracklayer)
    library(data.table)
    library(stringr)
    library(tools)
  })
  
  # Read BIM file    
  if (is.character(bim_file)) {
    annot <- vroom::vroom(bim_file, delim = '\t', col_names = FALSE)
    colnames(annot) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
    annot <- annot[, c("SNP", "CHR", "BP", "A1", "A2")]
  } else if (is.data.frame(bim_file)) {
    annot <- bim_file
    # Check required columns
    required_cols <- c("SNP", "CHR", "BP", "A1", "A2")
    if (!all(required_cols %in% colnames(annot))) {
      stop("bim_file data frame must contain columns: SNP, CHR, BP, A1, A2")
    }
  } else {
    stop("bim_file must be a data frame or a file path")
  }
  
  # Check for NAs introduced by coercion
  if (any(is.na(annot$BP))) {
    stop("Some BP values could not be coerced to integers. Please check the BIM file for formatting issues.")
  }
  
  # Get BED file paths
  if (is.null(bed_files)) {
    bed_files_list <- NULL
  } else if (is.character(bed_files)) {
    # If bed_files is a directory
    if (length(bed_files) == 1 && dir.exists(bed_files)) {
      bed_files_list <- list.files(bed_files, pattern = '\\.bed$', full.names = TRUE)
    } else {
      # Assume bed_files is a vector of file paths
      bed_files_list <- bed_files
    }
  } else {
    stop("bed_files must be NULL, a directory path, or a vector of BED file paths")
  }
  
  # Process BED files if provided
  if (!is.null(bed_files_list) && length(bed_files_list) > 0) {
    cat("Importing BED files and creating annotations...\n")
    
    # Use the name of the first BED file as the annotation name.
    annot_name <- file_path_sans_ext(basename(bed_files_list[1]))
    
    # Initialize a list to store GRanges objects from each BED file
    bed_gr_list <- lapply(bed_files_list, function(f) {
      # Here we assume you have a function read_preprocess_bed() that reads a BED file and returns a GRanges object.
      gr <- tryCatch({
        read_preprocess_bed(f)
      }, error = function(e) {
        stop(paste("Error processing BED file:", f, ":", e$message))
      })
      return(gr)
    })
    
    # Combine all GRanges into one unified GRanges object
    combined_bed_gr <- do.call(c, bed_gr_list)
    
    # Convert SNP positions to GRanges
    snp_gr <- GRanges(seqnames = annot$CHR, 
                      ranges = IRanges(start = annot$BP, end = annot$BP),
                      SNP = annot$SNP)
    
    # Find overlaps between SNPs and the combined BED regions
    overlaps <- findOverlaps(snp_gr, combined_bed_gr, ignore.strand = TRUE)
    
    # Create a new binary column with the annotation name: 1 if SNP overlaps any BED region, 0 otherwise
    annot[[annot_name]] <- 0
    if (length(overlaps) > 0) {
      overlapping_indices <- unique(queryHits(overlaps))
      annot[[annot_name]][overlapping_indices] <- 1
    } else {
      warning("No overlaps found between SNPs and BED annotations.")
    }
    
    # Create the complement column, named <annot_name>_not
    complement_name <- paste0(annot_name, "_not")
    annot[[complement_name]] <- 1 - annot[[annot_name]]
    cat("Annotations", annot_name, "and its complement", complement_name, "created.\n")
    
  } else {
    # If no BED files provided, add default weights column
    annot$weights <- 1
    cat("No BED files provided. Added weights column.\n")
  }
  
  # Write combined annotations to outfile if specified
  if (!is.null(outfile)) {
    cat("Writing combined annotations to outfile...\n")
    vroom::vroom_write(annot, outfile)
  }
  
  # Compute M values and write to Mfile if specified
  if (!is.null(Mfile)) {
    cat("Computing and writing M values to Mfile...\n")
    # Identify annotation columns (all columns except SNP, CHR, BP, A1, and A2)
    annotation_cols <- setdiff(colnames(annot), c("SNP", "CHR", "BP", "A1", "A2"))
    M_values <- colSums(annot[, annotation_cols, drop = FALSE])
    writeLines(paste0(M_values, collapse = ' '), con = Mfile)
  }
  
  cat("Annotation process completed successfully.\n")
  return(annot)
}
