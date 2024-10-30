#' Generate BED Files from Topic Annotations
#'
#' This function generates BED files for each topic based on the provided topic annotations.
#'
#' @param topics_res A list containing the factor matrices \code{Fmat} and \code{Lmat}, 
#' differential expression results \code{de_res}, and p-values matrix \code{p_jk}.
#'    - topics_annot A matrix or data frame containing topic annotations (peaks x topics).
#'      Rows correspond to peaks with names in the format "chr:start-end".
#'      Columns correspond to topics.
#'    - topics_prob A matrix probability of peak is in a topic
#' @param output_dir A character string specifying the output directory for the BED files.
#' @param fuzzy Logical value indicating whether to create fuzzy annotations (default: FALSE).
#' @param cutoff Numeric value specifying the cutoff for binary annotations when `fuzzy = FALSE` (default: 0.5).
#' @return None. BED files are written to the specified output directory.
#' @importFrom stringr str_extract
#' @export
#' 
get_annot_beds <- function(topics_res, output_dir, fuzzy = FALSE, cutoff = 0.5) {
  
  # Obtain data from topics_res object from run_fastTopics step
  topics_annot <- topics_res$Fmat
  topics_prob <- topics_res$Pmat
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Ensure row names contain peak ranges in "chr:start-end" format
  peak_ranges <- rownames(topics_annot)
  if (is.null(peak_ranges)) {
    stop("Row names of topics_annot must contain peak ranges in 'chr:start-end' format.")
  }
  
  # Extract seqnames, start, and end from peak ranges
  seqnames <- stringr::str_extract(peak_ranges, "^[^:]+")
  positions <- stringr::str_extract(peak_ranges, "(?<=:).+")
  start <- as.numeric(stringr::str_extract(positions, "^[^-]+"))
  end <- as.numeric(stringr::str_extract(positions, "(?<=-).+"))
  
  # Create a data frame with seqnames, start, end
  bed_df <- data.frame(seqnames = gsub("chr", "", seqnames), 
                       start = start, 
                       end = end)
  
  # Add topic names
  topics <- paste0('k', 1:ncol(topics_prob))
  colnames(topics_prob) <- topics
  
  # If fuzzy annotations are requested
  if (fuzzy) {
    # Implement fuzzy annotation logic here
  } else {
    
    # Binary annotations based on cutoff
    topics_annot_bin <- ifelse(topics_prob >= cutoff, 1, 0)
    
    # For each topic, create a BED file
    for (topic in topics) {
      # Subset to peaks where the topic annotation is 1
      indices <- which(topics_annot_bin[, topic] == 1)
      if (length(indices) > 0) {
        
        df_subset <- bed_df[indices, ]
        
        # Write the BED file
        write.table(df_subset,
                    file = file.path(output_dir, paste0(topic, "_annotations.bed")),
                    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      }
    }
  }
  
  return(output_dir)
}



