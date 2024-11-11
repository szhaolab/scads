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
#' 
create_annotations <- function(bim_file, bed_files = NULL, baseline_file, outfile = NULL, Mfile = NULL) {
  
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
  
  # Define annotator function
  annotator <- function(bim_file, annotations = NULL){
    
    if (is.null(annotations)){
      bim_file$weights <- 1
      return(bim_file)
    }
    
    snpRanges <- GenomicRanges::GRanges(seqnames = bim_file$CHR, 
                                        ranges = IRanges::IRanges(start = bim_file$BP, end = bim_file$BP))
    snpRanges <- plyranges::mutate(snpRanges, SNP = bim_file$SNP)
    
    for(f in annotations){
      
      name <- basename(f)
      name <- strsplit(name, split='[.]')[[1]][1]
      curr <- rtracklayer::import(f, format = 'bed')
      subdf <- IRanges::subsetByOverlaps(snpRanges, curr)
      snpsIn <- unique(subdf$SNP)
      bim_file <- dplyr::mutate(bim_file, !!name := ifelse(SNP %in% snpsIn, 1, 0))
    }
    return(bim_file)
  }
  
  # Annotate SNPs
  base_annot <- annotator(annot, bed_files_list)
  
  # Read baseline annotations
  if (is.character(baseline_file)) {
    baseline_annots <- vroom::vroom(baseline_file)
  } else if (is.data.frame(baseline_file)) {
    baseline_annots <- baseline_file
  } else {
    stop("baseline_file must be a data frame or a file path")
  }
  
  # Merge with baseline annotations
  base_annot <- dplyr::inner_join(base_annot, baseline_annots, by = c('SNP','CHR','BP','A1','A2'))
  
  # Write combined annotations to outfile if specified
  if (!is.null(outfile)) {
    vroom::vroom_write(base_annot, outfile)
  }
  
  # Compute M values and write to Mfile if specified
  if (!is.null(Mfile)) {
    M_values <- unname(colSums(base_annot[, c(6:ncol(base_annot))]))
    writeLines(paste0(M_values, collapse = ' '), con = Mfile)
  }
  
  # return(base_annot)
}
