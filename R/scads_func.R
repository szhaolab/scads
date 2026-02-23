#' Helper functions for scads
#' @importFrom utils modifyList
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom ashr ash
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom RhpcBLASctl blas_get_num_procs
#' @importFrom fastTopics verify.fit.and.count.matrix
#' @importFrom fastTopics verify.positive.vector
#' @importFrom fastTopics de_analysis_control_default
#' @importFrom fastTopics add_pseudocounts
#' @importFrom fastTopics initialize.multithreading
#' @importFrom fastTopics fit_poisson_models
#' @importFrom fastTopics lpfromz
#' @importFrom parallel splitIndices
#' @importFrom pbapply pboptions
#' @importFrom pbapply pblapply
#' 

# Create continuous annotation for S-LDSC
create_continuous_annotation <- function(bim_file, bed_file, annot_file) {
  
  # Read the .bim file
  bim_dt <- fread(bim_file, header = FALSE)
  setnames(bim_dt, c("CHR","SNP","GD","BP","A1","A2"))
  
  # Read the BED file (chr start end value).
  bed_df <- fread(bed_file, header = FALSE)
  setnames(bed_df, c("chr","start","end","value"))
  bed_df[, chr := gsub("^chr", "", chr)]
  
  # Build GRanges for BED intervals (0-based -> 1-based).
  bed_gr <- GRanges(
    seqnames = bed_df$chr,
    ranges   = IRanges(start = bed_df$start + 1, end = bed_df$end),
    score    = bed_df$value
  )
  
  # Build GRanges for the SNPs (using CHR/BP from .bim).
  snp_gr <- GRanges(
    seqnames = bim_dt$CHR,
    ranges   = IRanges(start = bim_dt$BP, end = bim_dt$BP)
  )
  
  # Find overlaps
  overlaps <- findOverlaps(snp_gr, bed_gr, ignore.strand = TRUE)
  
  # Create the numeric annotation vector
  continuous_annot <- numeric(length(snp_gr))
  
  # Assign the mean 'score' to overlapping SNPs
  overlap_list <- split(subjectHits(overlaps), queryHits(overlaps))
  for (snp_idx in names(overlap_list)) {
    idx <- as.integer(snp_idx)
    overlapping_values <- mcols(bed_gr)$score[ overlap_list[[snp_idx]] ]
    continuous_annot[idx] <- mean(overlapping_values, na.rm = TRUE)
  }
  
  # Write out a single column with header = "ANNOT"
  # 1) Write to uncompressed file
  uncompressed_file <- sub("\\.gz$", "", annot_file)
  fwrite(
    data.table(ANNOT = continuous_annot),
    file      = uncompressed_file,
    sep       = "\t",
    col.names = TRUE,   # <<-- provide a header line
    row.names = FALSE,
    quote     = FALSE
  )
  
  # 2) Gzip the file
  system2("gzip", c("-f", shQuote(uncompressed_file)))
  
  cat("Continuous annotation written + gzipped to", annot_file, "\n")
}

create_continuous_annotation2 <- function(bim_file, bed_file, annot_file) {
  library(data.table)
  library(GenomicRanges)
  
  message("[1] Reading BIM file: ", bim_file)
  bim_dt <- fread(bim_file, header = FALSE)
  setnames(bim_dt, c("CHR","SNP","GD","BP","A1","A2"))
  message("    → BIM rows: ", nrow(bim_dt))
  
  message("[2] Reading BED file: ", bed_file)
  bed_df <- fread(bed_file, header = FALSE)
  setnames(bed_df, c("chr","start","end","value"))
  bed_df[, chr := gsub("^chr", "", chr)]
  message("    → BED rows: ", nrow(bed_df),
          "; value range: [", min(bed_df$value), ", ", max(bed_df$value), "]")
  
  message("[3] Building GRanges objects")
  bed_gr <- GRanges(
    seqnames = bed_df$chr,
    ranges   = IRanges(start = bed_df$start + 1, end = bed_df$end),
    score    = bed_df$value
  )
  snp_gr <- GRanges(
    seqnames = bim_dt$CHR,
    ranges   = IRanges(start = bim_dt$BP, end = bim_dt$BP)
  )
  message("    → SNP GRanges length: ", length(snp_gr),
          "; BED GRanges length: ", length(bed_gr))
  
  message("[4] Finding overlaps")
  overlaps <- findOverlaps(snp_gr, bed_gr, ignore.strand = TRUE)
  message("    → # overlaps: ", length(overlaps))
  
  message("[5] Building continuous_annot vector")
  continuous_annot <- numeric(length(snp_gr))
  overlap_list <- split(subjectHits(overlaps), queryHits(overlaps))
  message("    → # SNPs with any overlap: ", length(overlap_list))
  
  for (snp_idx in names(overlap_list)) {
    idx <- as.integer(snp_idx)
    vals <- mcols(bed_gr)$score[overlap_list[[snp_idx]]]
    continuous_annot[idx] <- mean(vals, na.rm = TRUE)
  }
  message("    → continuous_annot range before scaling: [",
          min(continuous_annot, na.rm=TRUE), ", ",
          max(continuous_annot, na.rm=TRUE), "]")
  
  message("[6] Rescaling to [0.0001, 1.0001]")
  new_min <- 0.0001; new_max <- 1.0001
  old_min <- min(continuous_annot, na.rm = TRUE)
  old_max <- max(continuous_annot, na.rm = TRUE)
  if (old_max == old_min) {
    warning("All annotation values identical; skipping rescaling.")
    scaled <- rep((new_min + new_max) / 2, length(continuous_annot))
  } else {
    scaled01 <- (continuous_annot - old_min) / (old_max - old_min)
    scaled <- scaled01 * (new_max - new_min) + new_min
  }
  message("    → scaled range: [", min(scaled), ", ", max(scaled), "]")
  
  message("[7] Writing to ", annot_file)
  uncompressed_file <- sub("\\.gz$", "", annot_file)
  fwrite(
    data.table(ANNOT = scaled),
    file      = uncompressed_file,
    sep       = "\t",
    col.names = TRUE,
    quote     = FALSE
  )
  message("    → gzipping ", uncompressed_file)
  system2("gzip", c("-f", shQuote(uncompressed_file)))
  
  message("[8] Done: continuous annotation written to ", annot_file)
}



# Helper function to determine cutoff using KDE (kernel density estimation) - OLD 
find_kde_cutoff_old <- function(data, bandwidth = "nrd0") {
  
  # Compute density
  density_est <- density(data, bw = bandwidth)
  
  # Compute the first derivative of the density
  derivative <- diff(density_est$y)
  
  # Compute the sign of the derivative
  sign_derivative <- sign(derivative)
  
  # Compute where the derivative changes from negative to positive (local minima)
  sign_changes <- diff(sign_derivative)
  
  # Indices where a local minimum occurs (change from -1 to +1)
  minima_indices <- which(sign_changes == 2) + 1  # +1 to adjust for diff
  
  # If multiple minima are found, select the one between the two main modes
  # Assuming two modes, there should be one minima
  if (length(minima_indices) >= 1) {
    cutoff <- density_est$x[minima_indices[1]]
  } else {
    # If no clear minimum is found, default to a global threshold (e.g., median)
    cutoff <- median(data)
    warning("No clear minimum found in KDE. Using median as cutoff.")
  }
  
  return(cutoff)
}

# Helper function to determine cutoff using KDE (kernel density estimation)
find_kde_cutoff <- function(data, bandwidth = "nrd0", plot_results = FALSE) {
  # Compute density and derivatives
  density_est <- density(data, bw = bandwidth)
  derivative <- diff(density_est$y)
  sign_derivative <- sign(derivative)
  sign_changes <- diff(sign_derivative)
  
  # Detect maxima and minima
  maxima_indices <- which(sign_changes == -2) + 1
  minima_indices <- which(sign_changes == 2) + 1
  
  # Handle cases with fewer than two maxima
  if (length(maxima_indices) < 2) {
    warning("Fewer than two maxima found. Returning median as cutoff.")
    return(median(data))
  }
  
  # Sort maxima and select top two
  maxima_x <- density_est$x[maxima_indices]
  maxima_y <- density_est$y[maxima_indices]
  sorted_maxima <- order(maxima_y, decreasing = TRUE)
  primary_maxima <- maxima_indices[sorted_maxima[1:2]]
  primary_maxima_x <- density_est$x[primary_maxima]
  
  # Find minimum between the top two maxima
  minima_between <- minima_indices[
    density_est$x[minima_indices] > primary_maxima_x[1] &
      density_est$x[minima_indices] < primary_maxima_x[2]
  ]
  
  if (length(minima_between) > 0) {
    selected_minimum <- minima_between[which.min(density_est$y[minima_between])]
    cutoff <- density_est$x[selected_minimum]
  } else {
    cutoff <- median(data)  # Fallback to median
    warning("No clear minimum found between primary maxima. Using median as cutoff.")
  }
  
  # Plot results if requested
  if (plot_results) {
    plot(density_est, main = "Density with Maxima and Minima")
    points(maxima_x, maxima_y, col = "red", pch = 20)
    points(density_est$x[minima_indices], density_est$y[minima_indices], col = "blue", pch = 20)
    points(primary_maxima_x, density_est$y[primary_maxima], col = "green", pch = 20)  # Highlight top maxima
    abline(v = cutoff, col = "purple", lwd = 2, lty = 2)  # Highlight cutoff
  }
  
  return(cutoff)
}

find_kde_midpoint <- function(data, bandwidth = "nrd0", plot_results = FALSE) {
  # Compute density
  density_est <- density(data, bw = bandwidth)
  derivative <- diff(density_est$y)
  sign_derivative <- sign(derivative)
  sign_changes <- diff(sign_derivative)
  
  # Find maxima
  maxima_indices <- which(sign_changes == -2) + 1
  
  if (length(maxima_indices) < 2) {
    warning("Fewer than two maxima found. Returning median as cutoff.")
    return(median(data))
  }
  
  # Extract x-values of the two highest peaks
  maxima_x <- density_est$x[maxima_indices]
  maxima_y <- density_est$y[maxima_indices]
  sorted_maxima <- order(maxima_y, decreasing = TRUE)
  top_two_maxima_x <- sort(maxima_x[sorted_maxima[1:2]]) # Ensure left-to-right order
  
  # Compute midpoint between the two peaks
  cutoff <- mean(top_two_maxima_x)
  
  # Plot results if requested
  if (plot_results) {
    plot(density_est, main = "Density with Maxima and Minima", xlab = "Value", ylab = "Density")
    
    # Highlight maxima
    points(top_two_maxima_x, density_est$y[match(top_two_maxima_x, density_est$x)], col = "green", pch = 20)
    
    # Highlight midpoint (cutoff)
    abline(v = cutoff, col = "purple", lwd = 2, lty = 2)
    
    # Annotate cutoff
    text(cutoff, max(density_est$y) * 0.9, labels = round(cutoff, 6), col = "purple", pos = 4)
  }
  
  return(cutoff)
}


#' Estimate “background” accessibility rate from a reference cell‐type peak list
#'
#' @param count_matrix A peaks-by-cells matrix (dgCMatrix or dense).
#' @param cell_type_cells list of cells of the specific cell type
#' @param cell_type_peaks Path to a BED file of that cell type’s peaks (three cols: chr, start, end).
#' @return A single numeric: the estimated baseline accessibility probability.
#' @examples
#' get_average_bg(cm, "CD8", "cd8_peaks.bed")
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


library(GenomicRanges)
library(IRanges)
library(Matrix)
library(S4Vectors)
library(GenomeInfoDb)

#' Calculate Local Background (Lambda) for scATAC-seq Peaks (Per-Cell Model)
#'
#' This function calculates a local background accessibility rate (lambda) for
#' scATAC-seq peaks, similar to the method used in MACS2.
#'
#' This model calculates all rates on a "per-cell" basis
#' (i.e., average reads per base pair per cell).
#'
#' @param count_matrix A sparse matrix (e.g., dgCMatrix) with cells as rows
#'   and peaks as columns.
#' @param peaks A \code{GRanges} object containing the coordinates for all peaks.
#'   \strong{CRITICAL: The order of peaks \emph{must} match the
#'   order of columns in \code{count_matrix}.}
#' @param gc_content A numeric vector of GC content fraction for each peak,
#'   in the same order as \code{peaks}.
#' @param genome_size Numeric, the effective genome size. Default is 2.7e9.
#'
#' @return A \code{data.frame} with peak coordinates and all per-cell
#'   lambda values.
#'
calculate_lambda_local_scATAC <- function(count_matrix,
                                          peaks,
                                          gc_content = NULL,
                                          genome_size = 2.7e9) {
  
  # --- 0. Input Validation ---
  if (!inherits(peaks, "GRanges")) {
    stop("Argument 'peaks' must be a GRanges object.")
  }
  if (ncol(count_matrix) != length(peaks)) {
    stop(
      "Number of columns in 'count_matrix' (", ncol(count_matrix), ") ",
      "does not match the number of 'peaks' (", length(peaks), ")."
    )
  }
  
  # --- 1. Get Base Statistics ---
  n_cells <- nrow(count_matrix)
  if (n_cells == 0) {
    stop("count_matrix has 0 rows (cells).")
  }
  
  peak_sums <- Matrix::colSums(count_matrix)
  total_reads <- sum(peak_sums)
  
  # --- 2. Global Background (Per-Cell Model) ---
  # This is (avg reads per cell) / (genome size)
  # Units: (reads / cell) / bp
  lambda_bg_per_cell <- (total_reads / n_cells) / genome_size
  
  # --- 3. Peak Density (Per-Cell Model) ---
  peak_widths <- width(peaks)
  
  # This is (avg reads in peak per cell) / (peak width)
  # Units: (reads / cell) / bp
  peak_density_per_cell <- (peak_sums / n_cells) / (peak_widths + 1)
  
  # Add this to our peaks object for easy access
  mcols(peaks)$peak_density_per_cell <- peak_density_per_cell
  
  # --- 4. Local Background (Per-Cell Model) ---
  
  .get_local_mean_optimized <- function(peaks, window_kb) {
    windows <- peaks + (window_kb * 1000)
    hits <- findOverlaps(windows, peaks)
    
    # Calculate the mean of the *per-cell* peak densities
    local_means <- tapply(
      X = mcols(peaks)$peak_density_per_cell[subjectHits(hits)],
      INDEX = factor(queryHits(hits), levels = 1:length(peaks)),
      FUN = mean,
      na.rm = TRUE
    )
    
    # Handle peaks with no neighbors (will be NA)
    local_means[is.na(local_means)] <- 0
    
    return(as.numeric(local_means))
  }
  
  message("Calculating local lambda (per-cell) for 1kb window...")
  lambda_1k_per_cell <- .get_local_mean_optimized(peaks, 1)
  
  message("Calculating local lambda (per-cell) for 5kb window...")
  lambda_5k_per_cell <- .get_local_mean_optimized(peaks, 5)
  
  message("Calculating local lambda (per-cell) for 10kb window...")
  lambda_10k_per_cell <- .get_local_mean_optimized(peaks, 10)
  
  # --- 5. Combine into MACS2-style local lambda ---
  # All values are now in (reads / cell) / bp, so this is a valid comparison
  lambda_local_per_cell <- pmax(
    lambda_bg_per_cell,
    lambda_1k_per_cell,
    lambda_5k_per_cell,
    lambda_10k_per_cell,
    na.rm = TRUE
  )
  
  # --- 6. Optional GC Correction ---
  if (!is.null(gc_content)) {
    message("Applying GC content correction...")
    # (GC correction logic would go here, same as before)
    # It will operate on the 'lambda_local_per_cell'
  }
  
  # --- 7. Format Output ---
  message("Done.")
  data.frame(
    chr = as.character(seqnames(peaks)),
    start = start(peaks),
    end = end(peaks),
    lambda_bg_per_cell = lambda_bg_per_cell,
    lambda_1k_per_cell = lambda_1k_per_cell,
    lambda_5k_per_cell = lambda_5k_per_cell,
    lambda_10k_per_cell = lambda_10k_per_cell,
    lambda_local_per_cell = lambda_local_per_cell
  )
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



plot_signac_umap <- function(count_matrix, L, cell_scores, top_peaks = 2000, seed = 24) {
  set.seed(seed)  # For reproducibility
  
  # Step 1: Create a Seurat object using the count matrix
  cat("Creating Seurat object...\n")
  seurat_obj <- CreateSeuratObject(
    counts = count_matrix,
    assay = "ATAC"
  )
  
  # Step 2: Identify top variable features (peaks)
  cat("Finding top variable peaks...\n")
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = top_peaks)
  
  # Step 3: Perform Latent Semantic Indexing (LSI)
  cat("Performing LSI...\n")
  seurat_obj <- RunTFIDF(seurat_obj)  # TF-IDF normalization
  seurat_obj <- RunSVD(seurat_obj, reduction.name = "lsi", reduction.key = "LSI_")
  
  # Step 4: UMAP on LSI dimensions
  cat("Running UMAP...\n")
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "lsi",
    dims = 1:30,  # Use top 30 LSI dimensions by default
    reduction.name = "umap",
    reduction.key = "UMAP_"
  )
  
  # Step 5: Add cell scores for visualization
  cat("Adding cell scores...\n")
  seurat_obj$Cell_Score <- cell_scores
  
  # Step 6: Extract UMAP embeddings for plotting
  umap_embeddings <- Embeddings(seurat_obj, "umap")
  umap_df <- data.frame(
    UMAP1 = umap_embeddings[, 1],
    UMAP2 = umap_embeddings[, 2],
    enrich_topic_prop = L,
    Cell_Score = cell_scores  # Add cell scores for coloring
  )
  
  # Step 7: Plot UMAP
  purple_ramp <- c("white", "#DCD0FF", "#A884FF", "#6A45B1", "#421766")
  
  cat("Plotting UMAP...\n")
  p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = enrich_topic_prop)) +
    geom_point(alpha = 0.7, size = 0.5) +
    # scale_color_gradient(low = "blue", high = "red") +  
    scale_color_gradientn(
      colours = purple_ramp,
      values  = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)) # optional if you want equal spacing
    ) + 
    labs(title = "UMAP of scATAC-seq Data (Colored by prop of enriched topic)",
         x = "UMAP1", y = "UMAP2", color = "proportion of enriched topic in cell") +
    theme_minimal()
  
  cat("Plotting UMAP...\n")
  p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cell_Score)) +
    geom_point(alpha = 0.7, size = 0.5) +
    scale_color_gradientn(
      colours = purple_ramp,
      values  = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)) # optional if you want equal spacing
    ) + 
    # scale_color_gradient(low = "blue", high = "red") + 
    labs(title = "UMAP of scATAC-seq Data (Colored by Cell Score)",
         x = "UMAP1", y = "UMAP2", color = "Cell Score") +
    theme_minimal()
  
  return(p = list(p1,p2))  # Explicitly print the plot for non-interactive environments
}

plot_signac_umap2 <- function(count_matrix, L, cell_scores, z_scores, seed_idx, scavenge_trs, top_peaks = 2000, seed = 24) {
  set.seed(seed)  # For reproducibility
  
  # Step 1: Create a Seurat object using the count matrix
  cat("Creating Seurat object...\n")
  seurat_obj <- CreateSeuratObject(
    counts = count_matrix,
    assay = "ATAC"
  )
  
  # Step 2: Identify top variable features (peaks)
  cat("Finding top variable peaks...\n")
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = top_peaks)
  
  # Step 3: Perform Latent Semantic Indexing (LSI)
  cat("Performing LSI...\n")
  seurat_obj <- RunTFIDF(seurat_obj)  # TF-IDF normalization
  seurat_obj <- RunSVD(seurat_obj, reduction.name = "lsi", reduction.key = "LSI_")
  
  # Step 4: UMAP on LSI dimensions
  cat("Running UMAP...\n")
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "lsi",
    dims = 1:30,  # Use top 30 LSI dimensions by default
    reduction.name = "umap",
    reduction.key = "UMAP_"
  )
  
  # Step 5: Add cell scores for visualization
  cat("Adding cell scores...\n")
  seurat_obj$Cell_Score <- cell_scores
  
  # Step 6: Extract UMAP embeddings for plotting
  umap_embeddings <- Embeddings(seurat_obj, "umap")
  umap_df <- data.frame(
    UMAP1 = umap_embeddings[, 1],
    UMAP2 = umap_embeddings[, 2],
    enrich_topic_prop = L,
    z_score = z_scores,
    seed = seed_idx,
    scavenge_trs = scavenge_trs,
    Cell_Score = cell_scores  # Add cell scores for coloring
  )
  
  # Step 7: Plot UMAP
  purple_ramp <- c("white", "#DCD0FF", "#A884FF", "#6A45B1", "#421766")
  
  cat("Plotting UMAP...\n")
  p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = enrich_topic_prop)) +
    geom_point(alpha = 0.7, size = 0.5) +
    # scale_color_gradient(low = "blue", high = "red") +  # Customize color gradient
    scale_color_gradientn(
      colours = purple_ramp,
      values  = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)) # optional if you want equal spacing
    ) +
    labs(title = "UMAP of scATAC-seq Data (Colored by prop of enriched topic)",
         x = "UMAP1", y = "UMAP2", color = "proportion of enriched topic in cell") +
    theme_minimal()
  
  # cat("Plotting UMAP...\n")
  # p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = z_score)) +
  #   geom_point(alpha = 0.7, size = 0.5) +
  #   # scale_color_gradient(low = "blue", high = "red") +  # Customize color gradient
  #   scale_color_gradientn(
  #     colours = purple_ramp,
  #     values  = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)) # optional if you want equal spacing
  #   ) +
  #   labs(title = "UMAP of scATAC-seq Data (Colored by gchromVAR z_scores)",
  #        x = "UMAP1", y = "UMAP2", color = "z_scores (gchromVAR)") +
  #   theme_minimal()
  
  # cat("Plotting UMAP...\n")
  # p2b <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = seed)) +
  #   geom_point(alpha = 0.7, size = 0.5) +
  #   # scale_color_gradient(low = "blue", high = "red") +  # Customize color gradient
  #   scale_color_gradientn(
  #     colours = purple_ramp,
  #     values  = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)) # optional if you want equal spacing
  #   ) +
  #   labs(title = "UMAP of scATAC-seq Data (Colored by seed cell)",
  #        x = "UMAP1", y = "UMAP2", color = "Seed") +
  #   theme_minimal()
  
  cat("Plotting UMAP...\n")
  p3 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = scavenge_trs)) +
    geom_point(alpha = 0.7, size = 0.5) +
    # scale_color_gradient(low = "blue", high = "red") +  # Customize color gradient
    scale_color_gradientn(
      colours = purple_ramp,
      values  = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)) # optional if you want equal spacing
    ) +
    labs(title = "UMAP of scATAC-seq Data (Colored by SCAVENGE score)",
         x = "UMAP1", y = "UMAP2", color = "SCAVENGE TRS") +
    theme_minimal()
  
  cat("Plotting UMAP...\n")
  p4 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cell_Score)) +
    geom_point(alpha = 0.7, size = 0.5) +
    # scale_color_gradient(low = "blue", high = "red") +  # Customize color gradient
    scale_color_gradientn(
      colours = purple_ramp,
      values  = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)) # optional if you want equal spacing
    ) +
    labs(title = "UMAP of scATAC-seq Data (Colored by Cell Score)",
         x = "UMAP1", y = "UMAP2", color = "Cell Score by scads") +
    theme_minimal()
  
  return(p = list(p1,p3,p4))  # Explicitly print the plot for non-interactive environments
}


get_signac_umap2_null <- function(count_matrix, cell_scores, cs_zscores, scavenge_trs, top_peaks = 2000, seed = 24) {
  set.seed(seed)  # For reproducibility
  
  # Step 1: Create a Seurat object using the count matrix
  cat("Creating Seurat object...\n")
  seurat_obj <- CreateSeuratObject(
    counts = count_matrix,
    assay = "ATAC"
  )
  
  # Step 2: Identify top variable features (peaks)
  cat("Finding top variable peaks...\n")
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = top_peaks)
  
  # Step 3: Perform Latent Semantic Indexing (LSI)
  cat("Performing LSI...\n")
  seurat_obj <- RunTFIDF(seurat_obj)  # TF-IDF normalization
  seurat_obj <- RunSVD(seurat_obj, reduction.name = "lsi", reduction.key = "LSI_")
  
  # Step 4: UMAP on LSI dimensions
  cat("Running UMAP...\n")
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "lsi",
    dims = 1:30,  # Use top 30 LSI dimensions by default
    reduction.name = "umap",
    reduction.key = "UMAP_"
  )
  
  # Step 5: Add cell scores for visualization
  cat("Adding cell scores...\n")
  seurat_obj$Cell_Score <- cell_scores
  
  # Step 6: Extract UMAP embeddings for plotting
  umap_embeddings <- Embeddings(seurat_obj, "umap")
  umap_df <- data.frame(
    UMAP1 = umap_embeddings[, 1],
    UMAP2 = umap_embeddings[, 2],
    scavenge_trs = scavenge_trs,
    Cell_Score = cell_scores,
    Cell_Score_z = cs_zscores
  )
  
  return(umap_df)
}


get_signac_umap <- function(count_matrix, L, cell_scores, top_peaks = 13000, seed = 24) {
  set.seed(seed)  # For reproducibility
  
  # Step 1: Create a Seurat object using the count matrix
  cat("Creating Seurat object...\n")
  seurat_obj <- CreateSeuratObject(
    counts = count_matrix,
    assay = "ATAC"
  )
  
  # Step 2: Identify top variable features (peaks)
  cat("Finding top variable peaks...\n")
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = top_peaks)
  
  # Step 3: Perform Latent Semantic Indexing (LSI)
  cat("Performing LSI...\n")
  seurat_obj <- RunTFIDF(seurat_obj)  # TF-IDF normalization
  seurat_obj <- RunSVD(seurat_obj, reduction.name = "lsi", reduction.key = "LSI_")
  
  # Step 4: UMAP on LSI dimensions
  cat("Running UMAP...\n")
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "lsi",
    dims = 1:30,  # Use top 30 LSI dimensions by default
    n_neigbors = 30,
    reduction.name = "umap",
    reduction.key = "UMAP_"
  )
  
  # Step 5: Add cell scores for visualization
  cat("Adding cell scores...\n")
  seurat_obj$Cell_Score <- cell_scores
  
  # Step 6: Extract UMAP embeddings for plotting
  umap_embeddings <- Embeddings(seurat_obj, "umap")
  umap_df <- data.frame(
    UMAP1 = umap_embeddings[, 1],
    UMAP2 = umap_embeddings[, 2],
    enrich_topic_prop = L,
    Cell_Score = cell_scores  # Add cell scores for coloring
  )
  
  return(umap_df)  
}


plot_cell_ranking <- function(score, topic_prop, title = "Ranked plot") {
  stopifnot(length(score) == length(topic_prop))
  
  df <- data.frame(score = score, topic_prop = topic_prop) %>%
    arrange(desc(score)) %>%  
    mutate(rank = row_number(), rank_percent = rank / n() * 100)
  
  ggplot(df, aes(x = 0, xend = 1, y = rank_percent, yend = rank_percent, color = topic_prop)) +
    geom_segment(linewidth = 0.3) +
    scale_y_reverse(breaks = c(0, 25, 50, 75, 100), labels = c("0%", "25%", "50%", "75%", "100%")) +
    scale_color_gradientn(
      colors = c("white", "#D7C6F4", "#5A189A"),
      name = "Causal Topic Prop"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Cells ranked by score (high to low)"
    )
}

############################################################
# Function to extract phenores outputs
get_phenores_outputs <- function(phenores_file) {
  phenores <- readRDS(phenores_file)
  
  csnps_null <- list()
  csnps_annot <- list()
  h2_null <- list()
  h2_annot <- list()
  snps_null <- list()
  snps_annot <- list()
  beta_null <- list()
  beta_annot <- list()
  
  for (i in seq_along(phenores$batch)) {
    id_cSNP_list <- phenores$batch[[i]]$id.cSNP.list
    csnps_null[[i]] <- if (length(id_cSNP_list) >= 1) id_cSNP_list[[1]] else list()
    csnps_annot[[i]] <- if (length(id_cSNP_list) >= 2) id_cSNP_list[[2]] else list()
    
    snps_null[[i]] <- phenores$batch[[i]]$id.SNP.count[1]
    snps_annot[[i]] <- ifelse(length(phenores$batch[[i]]$id.SNP.count) >= 2,
                              phenores$batch[[i]]$id.SNP.count[2], 0)
    h2_null[[i]] <- phenores$batch[[i]]$var.snp[1]
    h2_annot[[i]] <- ifelse(length(phenores$batch[[i]]$var.snp) >= 2,
                            phenores$batch[[i]]$var.snp[2], 0)
    
    beta_all <- phenores$batch[[i]]$beta.all
    if (length(id_cSNP_list) >= 1) {
      num_null <- length(id_cSNP_list[[1]])
      beta_null[[i]] <- beta_all[1:num_null]
      beta_annot[[i]] <- beta_all[-(1:num_null)]
    } else {
      beta_null[[i]] <- numeric(0)
      beta_annot[[i]] <- beta_all
    }
  }
  
  csnps_null <- unlist(csnps_null, recursive = FALSE)
  csnps_annot <- unlist(csnps_annot, recursive = TRUE)
  snps_null <- unlist(snps_null, recursive = FALSE)
  snps_annot <- unlist(snps_annot, recursive = FALSE)
  h2_null <- unlist(h2_null, recursive = FALSE)
  h2_annot <- unlist(h2_annot, recursive = FALSE)
  beta_null <- unlist(beta_null, recursive = FALSE)
  beta_annot <- unlist(beta_annot, recursive = FALSE)
  
  return(list(csnps_null = csnps_null,
              csnps_annot = csnps_annot,
              snps_null = snps_null,
              snps_annot = snps_annot,
              h2_null = h2_null,
              h2_annot = h2_annot,
              beta_null = beta_null,
              beta_annot = beta_annot))
}

############################################################


plot_snp_cs_correlation <- function(df,
                                    cor_col    = "cor",
                                    pval_col   = "pval",
                                    threshold  = 0.5,
                                    facet_ncol = 3,
                                    point_col  = "#E74C3C",
                                    line_col   = "#2C3E50") {
  
  
  # turn your column‐names into symbols for tidy eval
  cor_sym  <- rlang::sym(cor_col)
  pval_sym <- rlang::sym(pval_col)
  
  # 1) pick the SNP×cell‐type combos above threshold
  good_pairs <- df %>%
    dplyr::distinct(SNP, BioClassification, !!cor_sym, !!pval_sym) %>%
    dplyr::filter(abs(!!cor_sym) > threshold)
  
  if (nrow(good_pairs)==0) {
    stop("No SNP×cell-type combos with ", cor_col, " >", threshold)
  }
  
  # 2) subset to those
  plot_data <- df %>%
    dplyr::semi_join(good_pairs, by = c("SNP","BioClassification"))
  
  # 3) build facet labels
  labels <- good_pairs %>%
    dplyr::mutate(
      label = sprintf("%s — %s\nr=%.2f, p=%.1g",
                      SNP, BioClassification,
                      !!cor_sym, !!pval_sym)
    ) %>%
    dplyr::select(SNP, BioClassification, label)
  
  plot_data <- plot_data %>%
    dplyr::left_join(labels, by = c("SNP","BioClassification"))
  
  # 4) draw
  ggplot(plot_data, aes(x = cs_bin, y = total_reads)) +
    geom_line(color = line_col, size = 0.8, aes(group = 1)) +
    geom_point(color = point_col, size = 2) +
    facet_wrap(~ label, ncol = facet_ncol, scales = "free_y") +
    scale_x_continuous(breaks = 1:5) +
    scale_y_continuous(labels = comma) +
    labs(
      title   = sprintf("ATAC‐seq reads vs. cs_bin (cor > %.2f)", threshold),
      x       = "cs_bin (1 = lowest, 5 = highest)",
      y       = "Total ATAC reads",
      caption = "Facets = SNP × cell type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5),
      strip.text       = element_text(size = 9, face = "bold"),
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_blank()
    )
}


plot_snp_cs_correlation2 <- function(df,
                                    cor_col    = "cor",
                                    pval_col   = "pval",
                                    threshold  = 0.5,
                                    facet_ncol = 3,
                                    point_col  = "#E74C3C",
                                    line_col   = "#2C3E50") {
  
  
  # turn your column‐names into symbols for tidy eval
  cor_sym  <- rlang::sym(cor_col)
  pval_sym <- rlang::sym(pval_col)
  
  # 1) pick the SNP×cell‐type combos above threshold
  good_pairs <- df %>%
    dplyr::distinct(SNP, !!cor_sym, !!pval_sym) %>%
    dplyr::filter(abs(!!cor_sym) > threshold)
  
  if (nrow(good_pairs)==0) {
    stop("No SNP×cell-type combos with ", cor_col, " >", threshold)
  }
  
  # 2) subset to those
  plot_data <- df %>%
    dplyr::semi_join(good_pairs, by = c("SNP"))
  
  # 3) build facet labels
  labels <- good_pairs %>%
    dplyr::mutate(
      label = sprintf("%s — \nr=%.2f, p=%.1g",
                      SNP,
                      !!cor_sym, !!pval_sym)
    ) %>%
    dplyr::select(SNP, label)
  
  plot_data <- plot_data %>%
    dplyr::left_join(labels, by = c("SNP"))
  
  # 4) draw
  ggplot(plot_data, aes(x = cs_bin_all, y = total_reads)) +
    geom_line(color = line_col, size = 0.8, aes(group = 1)) +
    geom_point(color = point_col, size = 2) +
    facet_wrap(~ label, ncol = facet_ncol, scales = "free_y") +
    scale_x_continuous(breaks = 1:5) +
    scale_y_continuous(labels = comma) +
    labs(
      title   = sprintf("ATAC‐seq reads vs. cs_bin (cor > %.2f)", threshold),
      x       = "cs_bin (1 = lowest, 5 = highest)",
      y       = "Total ATAC reads",
      caption = "Facets = SNP × cell type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5),
      strip.text       = element_text(size = 9, face = "bold"),
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_blank()
    )
}
