#' Generate BED Files from Topic Annotations
#'
#' This function generates BED files for each topic based on the provided topic annotations.
#' Each BED file is saved into its own folder, and a list of output folder paths is returned.
#'
#' @param topics_res A list containing the factor matrices \code{Fmat} and \code{Lmat}, 
#' differential expression results \code{de_res}, and p-values matrix \code{p_jk}.
#'    - \code{Fmat}: A matrix containing factor loadings (peaks x topics).
#'    - \code{Pmat}: A matrix of p-values or probabilities (peaks x topics).
#'      Rows correspond to peaks with names in the format "chr:start-end".
#'      Columns correspond to topics.
#' @param output_dir A character string specifying the base output directory for the BED files.
#' @param fuzzy Logical value indicating whether to create fuzzy annotations (default: FALSE).
#' @param cutoff Numeric value specifying the cutoff for binary annotations when `fuzzy = FALSE` (default: 0.5).
#'
#' @return A list of output folder paths (`bed_dir_list`), one for each topic.
#' @importFrom stringr str_extract
#' @export
#' 
get_annot_beds <- function(topics_res, output_dir, fuzzy = FALSE, cutoff = 0.9) {
  library(stringr)
  
  topics_annot <- topics_res$Fmat
  topics_prob  <- topics_res$Pmat
  peak_ranges  <- rownames(topics_annot)
  head(peak_ranges)
  if (is.null(peak_ranges)) stop("Row names must contain peak ranges")
  
  # Detect format and split accordingly:
  if (all(grepl("^[^_]+_[0-9]+_[0-9]+$", peak_ranges))) {
    # underscore format: chr1_10234_10734
    parts <- strsplit(peak_ranges, "_", fixed = TRUE)
    seqnames <- vapply(parts, `[`, 1, FUN.VALUE = "")
    start    <- as.integer(vapply(parts, `[`, 2, FUN.VALUE = ""))
    end      <- as.integer(vapply(parts, `[`, 3, FUN.VALUE = ""))
    
  } else if (all(grepl("^[^:]+:[0-9]+-[0-9]+$", peak_ranges))) {
    # colon‐dash format: chr1:10234-10734
    seqnames <- stringr::str_extract(peak_ranges, "^[^:]+")
    positions <- stringr::str_extract(peak_ranges, "(?<=:).+")
    start    <- as.integer(stringr::str_extract(positions, "^[0-9]+"))
    end      <- as.integer(stringr::str_extract(positions, "[0-9]+$"))
    
  } else if (all(grepl("^[^\\-]+-[0-9]+-[0-9]+$", peak_ranges))) {
    # dash format: chr1-10234-10734
    parts <- strsplit(peak_ranges, "-", fixed = TRUE)
    seqnames <- vapply(parts, `[`, 1, FUN.VALUE = "")
    start    <- as.integer(vapply(parts, `[`, 2, FUN.VALUE = ""))
    end      <- as.integer(vapply(parts, `[`, 3, FUN.VALUE = ""))
    
  } else {
    stop("Row names must be one of: 'chr_start_end', 'chr:start-end', or 'chr-start-end'")
  }
  
  # Build bed_df
  bed_df <- data.frame(
    seqnames = ifelse(grepl("^chr", seqnames), seqnames, paste0("chr", seqnames)),
    start    = start,
    end      = end,
    stringsAsFactors = FALSE
  )
  
  # Add topic names
  topics <- paste0('k', 1:ncol(topics_prob))
  colnames(topics_prob) <- topics
  
  # Initialize list to store output folder paths
  bed_dir_list <- list()
  
  # If fuzzy annotations are requested
  if (fuzzy) {
    # Implement fuzzy annotation logic here
    stop("Fuzzy annotations are not yet implemented.")
  } else {
    
    # Binary annotations based on cutoff
    topics_annot_bin <- ifelse(topics_prob >= cutoff, 1, 0)
    
    # For each topic, create a BED file in its own folder
    for (topic in topics) {
      # Subset to peaks where the topic annotation is 1
      indices <- which(topics_annot_bin[, topic] == 1)
      if (length(indices) > 0) {
        
        df_subset <- bed_df[indices, ]
        
        # Create output folder for this topic
        output_folder <- file.path(output_dir, paste0(topic, "_output"))
        if (!dir.exists(output_folder)) {
          dir.create(output_folder, recursive = TRUE)
        }
        
        # Save the output folder path to the list
        bed_dir_list[[topic]] <- output_folder
        
        # Write the BED file into the topic's folder
        write.table(df_subset,
                    file = file.path(output_folder, paste0(topic, "_annotations.bed")),
                    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      } else {
        # If no peaks, still create the folder
        output_folder <- file.path(output_dir, paste0(topic, "_output"))
        if (!dir.exists(output_folder)) {
          dir.create(output_folder, recursive = TRUE)
        }
        bed_dir_list[[topic]] <- output_folder
        # Optionally, write an empty BED file or skip writing
        warning(paste("No peaks passed the cutoff for topic", topic, "- Empty BED file created."))
      }
    }
  }
  
  return(bed_dir_list)
}
