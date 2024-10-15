#' Get Significant Loci from Summary Statistics and LD Blocks
#'
#' This function identifies significant loci based on summary statistics and LD blocks.
#'
#' @param sumstats_path Directory to a data frame containing summary statistics. Must include columns:
#'   - `chr`: Chromosome number.
#'   - `pos`: Base pair position.
#'   - `snp`: SNP identifier.
#'   - `pval`: P-value of association.
#' @param ldblocks_path Directory to a data frame containing LD blocks. Must include columns:
#'   - `X1`: Chromosome number.
#'   - `X2`: Start position of the LD block.
#'   - `X3`: End position of the LD block.
#'   - `X4`: Locus identifier.
#' @param pval_threshold A numeric value specifying the p-value threshold for significance. Default is \code{5e-8}.
#' @param max_loci An integer specifying the maximum number of loci to return. Default is \code{150}.
#' @param outfile An optional character string specifying the output file path. If provided, the significant loci data frame will be written to this file.
#' @return A data frame of significant loci.
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom plyranges mutate join_overlap_inner
#' @importFrom dplyr inner_join group_by summarise filter arrange desc slice
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom vroom vroom
#' 
get_signif_loci <- function(sumstats_path, ldblocks_path, pval_threshold = 5e-8, max_loci = 150, outfile = NULL) {
  
  # Helper function to make GRanges objects
  make_ranges <- function(seqname, start, end){
    return(GenomicRanges::GRanges(seqnames = seqname, ranges = IRanges::IRanges(start = start, end = end)))
  }
  
  # Function to assign SNPs to loci
  assign_locus_snp <- function(cleaned_sumstats, ldBed){
    
    # Drop NA values and filter p-values greater than zero
    cleaned_sumstats <- cleaned_sumstats %>% tidyr::drop_na()
    cleaned_sumstats <- cleaned_sumstats[cleaned_sumstats$pval > 0, ]
    
    # Create GRanges for LD blocks
    ldRanges <- make_ranges(ldBed$X1, ldBed$X2, ldBed$X3)
    ldRanges <- ldRanges %>% plyranges::mutate(locus = ldBed$X4)
    
    # Create GRanges for SNPs
    snpRanges <- GenomicRanges::GRanges(seqnames = cleaned_sumstats$chr, 
                                        ranges = IRanges::IRanges(start = cleaned_sumstats$pos, 
                                                                  end   = cleaned_sumstats$pos,
                                                                  names = cleaned_sumstats$snp))
    snpRanges <- snpRanges %>% plyranges::mutate(snp = names(snpRanges))
    
    # Find overlapping SNPs and LD blocks
    snp_ld_overlap <- plyranges::join_overlap_inner(snpRanges, ldRanges)
    snp_ld_block <- tibble::as_tibble(mcols(snp_ld_overlap))
    
    # Remove duplicates
    snp_ld_block <- snp_ld_block[!duplicated(snp_ld_block$snp), ]
    
    # Merge with summary statistics
    cleaned_annot_sumstats <- dplyr::inner_join(cleaned_sumstats, snp_ld_block, by = 'snp')
    
    return(cleaned_annot_sumstats)
  }
  
  # Assign SNPs to loci
  sumstats <- vroom::vroom(sumstats_path, col_names = T)
  ldblocks <- vroom::vroom(ldblocks_path, col_names = F)
  cleaned_sumstats <- assign_locus_snp(sumstats, ldblocks)
  
  # Identify significant loci
  signif_loci <- cleaned_sumstats %>%
    dplyr::group_by(locus) %>%
    dplyr::summarise(n_signif = sum(pval <= pval_threshold)) %>%
    dplyr::filter(n_signif > 0)
  
  # Limit to top loci if necessary
  if (nrow(signif_loci) > max_loci) {
    signif_loci <- signif_loci %>%
      dplyr::arrange(dplyr::desc(n_signif)) %>%
      dplyr::slice(1:max_loci)
  }
  
  top_loci <- signif_loci$locus
  
  # Get LD blocks for significant loci
  sig_loci_df <- ldblocks[ldblocks$X4 %in% top_loci, ]
  
  # Write to output file if specified
  if (!is.null(outfile)) {
    vroom::vroom_write(sig_loci_df, outfile, col_names = F)
  }
  
  return(sig_loci_df)
}
