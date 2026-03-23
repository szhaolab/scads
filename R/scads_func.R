
#' Estimate “background” accessibility rate from a reference cell‐type peak list
#'
#' @param count_matrix A peaks-by-cells matrix (dgCMatrix or dense).
#' @param cell_type_cells list of cells of the specific cell type
#' @param cell_type_peaks Path to a BED file of that cell type’s peaks (three cols: chr, start, end).
#' @return A single numeric: the estimated baseline accessibility probability.
#' @examples
#' \dontrun{
#' get_average_bg(cm, "CD8", "cd8_peaks.bed")
#' }
get_average_bg <- function(count_matrix,
                           cell_type_cells,
                           cell_type_peaks) {

  cat("\nEstimating baseline using cell type\n")
  ## 1a) Convert CT peaks into a GRanges object
  CTpeaks <- fread(cell_type_peaks)
  gr_ct <- GRanges(
    seqnames = paste0("chr", CTpeaks$V1),
    ranges   = IRanges(start = CTpeaks$V2, end = CTpeaks$V3)
  )

  ## 1b) Convert all peaks in count marix into a GRanges object
  cells <- readRDS(cell_type_cells)
  print(head(cells))
  cm <- count_matrix[cells,]
  print(dim(cm))
  print(class(cm))
  # if (!is.matrix(cm)) {
  #   cm <- as.matrix(cm)
  # }
  # print(class(cm))
  peak_ids <- colnames(cm)
  print(head(peak_ids))
  ## Normalize separators: replace “:”, “-” or any mix with “_”
  peak_ids_norm <- gsub("[:\\-]", "_", peak_ids)
  # print(head(peak_ids_norm))
  ## Split into chromosome / start / end
  parts <- strsplit(peak_ids_norm, "_", fixed = TRUE)
  print(head(parts))
  # bind into a 3-column matrix
  mat <- do.call(rbind, parts)
  print(dim(mat))
  ## Turn into numeric where appropriate
  chr   <- mat[,1]
  start <- as.integer(mat[,2])
  end   <- as.integer(mat[,3])
  
  print("Start")
  print(head(start))
  print(sum(is.na(start)))
  print("End")
  print(head(end))
  print(sum(is.na(end)))

  ## Build your GRanges
  gr_peaks <- GenomicRanges::GRanges(seqnames = chr,
                      ranges   = IRanges::IRanges(start = start, end = end))

  ## 2) Find which peaks overlap any CT peaks
  is_overlap <- overlapsAny(gr_peaks, gr_ct, ignore.strand=TRUE)
  print(table(is_overlap))

  ## 3) Subset the counts‐matrix accordingly and compare total reads.

  # 3b) Sum *across cells* for each peak
  peak_totals <- colSums(cm)

  # 3c) Now split “overlapping” vs. “non‐overlapping” peaks
  overlapped_peaks     <- which(is_overlap)
  non_overlapped_peaks <- which(!is_overlap)

  sum_overlap_reads     <- sum( peak_totals[overlapped_peaks] )
  sum_non_overlap_reads <- sum( peak_totals[non_overlapped_peaks] )

  dt <- data.table::data.table(
    category    = c("overlap", "non_overlap"),
    n_peaks     = c(length(overlapped_peaks), length(non_overlapped_peaks)),
    total_reads = c(sum_overlap_reads, sum_non_overlap_reads)
  )
  print(dt)

  n_cells <- nrow(count_matrix)
  reads_per_cell <- (sum_overlap_reads + sum_non_overlap_reads)/n_cells # reads per cell
  baseline <- sum_non_overlap_reads/length(non_overlapped_peaks)/n_cells/reads_per_cell # 4.762821e-07
  cat("\nBaseline: ", baseline)
  
  return(baseline)
}



# Helper function to read and preprocess BED files
read_preprocess_bed <- function(f) {
  # Read BED file using fread
  bed_df <- tryCatch({
    fread(f, header = FALSE, sep = "\t", data.table = FALSE, fill = TRUE)
  }, error = function(e) {
    stop(paste("Error reading BED file:", f, ":", e$message))
  })
  
  # Ensure BED file has at least 3 columns: chrom, start, end
  if (ncol(bed_df) < 3) {
    stop(paste("BED file", f, "does not have at least 3 columns."))
  }
  
  # Rename columns for clarity
  colnames(bed_df)[1:3] <- c("chrom", "start", "end")
  
  # Convert start and end to integers (handle scientific notation)
  bed_df$start <- as.integer(as.numeric(bed_df$start))
  bed_df$end <- as.integer(as.numeric(bed_df$end))
  
  # Check for NAs introduced by coercion
  if (any(is.na(bed_df$start)) || any(is.na(bed_df$end))) {
    problematic_rows <- which(is.na(bed_df$start) | is.na(bed_df$end))
    stop(paste("Non-integer values found in BED file:", f, "at rows:", paste(problematic_rows, collapse = ", ")))
  }
  
  # Create GRanges object
  gr <- GRanges(seqnames = bed_df$chrom,
                ranges = IRanges(start = bed_df$start, end = bed_df$end))
  
  return(gr)
}

# Defensive helper function for formatting peak names
normalize_peak_names <- function(nms) {
  # 1. Strip whitespace
  nms <- trimws(nms)
  
  # 2. If format is "chr1_12345_67890" → "chr1:12345-67890"
  nms <- gsub("^(chr[^_]+)_([0-9]+)_([0-9]+)$", "\\1:\\2-\\3", nms)
  
  # 3. If format is "chr1:12345:67890" → "chr1:12345-67890"
  nms <- sub("^([^:]+:[0-9]+):([0-9]+)", "\\1-\\2", nms)
  
  # 4. Check what's still malformed after fixes
  valid_pattern <- "^[^:]+:[0-9]+-[0-9]+(:[+\\-\\*])?$"
  bad <- !grepl(valid_pattern, nms)
  if (any(bad)) {
    warning(sprintf(
      "%d rowname(s) could not be parsed as genomic ranges and will be dropped:\n%s",
      sum(bad),
      paste(head(nms[bad], 5), collapse = "\n")
    ))
  }
  
  nms
}
