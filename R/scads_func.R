#' Helper functions for scads
#' @importFrom utils modifyList
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom ashr ash
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom RhpcBLASctl blas_get_num_procs
#' 
#' 
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
  
  # ----- NEW: Min–max rescale into [0.0001, 1.0001] ----- #
  new_min <- 0.0001
  new_max <- 1.0001
  
  old_min <- min(continuous_annot, na.rm = TRUE)
  old_max <- max(continuous_annot, na.rm = TRUE)
  
  if (old_max == old_min) {
    warning("All annotation values are identical; skipping rescaling.")
    scaled <- rep((new_min + new_max) / 2, length(continuous_annot))
  } else {
    # first bring into [0, 1]
    scaled01 <- (continuous_annot - old_min) / (old_max - old_min)
    # then stretch into [new_min, new_max]
    scaled  <- scaled01 * (new_max - new_min) + new_min
  }
  
  # Write out a single column with header = "ANNOT"
  uncompressed_file <- sub("\\.gz$", "", annot_file)
  fwrite(
    data.table(ANNOT = scaled),
    file      = uncompressed_file,
    sep       = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote     = FALSE
  )
  
  # Gzip the file
  system2("gzip", c("-f", shQuote(uncompressed_file)))
  
  cat("Continuous annotation (min–max scaled) written + gzipped to", annot_file, "\n")
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


# Modified de_analysis function from fastTopics
de_analysis2 <- function (fit, X, s = rowSums(X), pseudocount = 0.01,
                          fit.method = c("scd","em","mu","ccd","glm"),
                          shrink.method = c("ash","none"), lfc.stat = "le",
                          control = list(), verbose = TRUE, ...) {
  
  # CHECK AND PROCESS INPUTS
  # ------------------------
  # Check and process input argument "fit".
  if (is.matrix(fit)) {
    L <- fit
    if (any((rowSums(L) - 1) > 1e-15))
      warning("\"fit\" is a matrix but may not be topic proportions matrix; ",
              "rowSums(fit) should be a vector of all ones")
    m <- ncol(X)
    k <- ncol(L)
    F <- matrix(1,m,k)
    rownames(F) <- colnames(X)
    colnames(F) <- colnames(L)
    fit <- init_poisson_nmf(X,F = F,L = L)
    fit <- poisson2multinom(fit)
  }
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\", or a matrix of topic proportions")
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  
  is.sparse.matrix <- function(x) is(x, 'sparseMatrix') # added
  
  # Check and process input argument "X".
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  fastTopics:::verify.fit.and.count.matrix(X,fit)
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"
  
  # Get the number of rows (n) and columns (m) in the counts matrix, and
  # the number of groups or topics (k).
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(fit$F)
  
  # Check input argument "s".
  fastTopics:::verify.positive.vector(s)
  if (length(s) != n)
    stop("Input argument \"s\" should be a vector of positive numbers, ",
         "in which length(s) = nrow(X)")
  
  # Check input argument "pseudocount".
  if (any(pseudocount <= 0))
    stop("Input argument \"pseudocount\" should be a positive number")
  
  # Process input arguments "fit.method" and "shrink.method".
  fit.method <- match.arg(fit.method)
  shrink.method <- match.arg(shrink.method)
  
  # Check and process input argument "control".
  control <- modifyList(fastTopics:::de_analysis_control_default(),control,keep.null = TRUE)
  if (control$nc > 1 & .Platform$OS.type == "windows")
    stop("Multithreading is not available on Windows; try again with ",
         "control$nc = 1")
  
  # Check and process input argument "lfc.stat".
  if (!(all(lfc.stat == "vsnull") | all(lfc.stat == "le"))) {
    if (!(any(lfc.stat == 1:k) | any(lfc.stat == colnames(fit$F))))
      stop("Input argument \"lfc.stat\" should be either \"vsnull\", \"le\" ",
           "a number between 1 and k, where k is the number of topics, or ",
           "a name of a topic (column of fit$F)")
    if (is.character(lfc.stat))
      lfc.stat <- match(lfc.stat,colnames(fit$F))
  }
  
  # FIT NULL MODELS - *****MODIFIED*****
  # ---------------
  # Compute the MLE f0 in the "null" model x ~ Poisson(u), with u =
  # s*f0. (This calculation is performed for each column of the counts
  # matrix. It must be done before adding pseudocounts to the data.)
  # Equivalently, this is the maximum-likelihood estimate of the
  # binomial probability in the "null" Binomial model x ~ Binom(s,p0).
  
  # background = f0
  # cat("Use f0 provided:", background, "\n")
  # f0 <- c(rep(background, times=ncol(X)))
  # names(f0) <- rownames(fit$F)
  
  # if (is.null(f0)) {
  #   # If f0 is not provided, use the default single cutoff
  #   background = 1/(3E9/500) # *100
  #   f0 <- c(rep(background, times=ncol(X)))
  #   names(f0) <- rownames(fit$F)
  #   cat("No f0 cutoffs provided. Using default background:", background, "\n")
  # } else {
  #   # Ensure f0 is a vector with length equal to number of topics
  #   if (!is.numeric(f0) | length(f0) != k) {
  #     stop("Input \"f0\" should be a numeric vector with length equal to the number of topics.")
  #   }
  #   names(f0) <- colnames(fit$F)
  #   cat("Using provided f0 values for each topic.\n")
  #   print(f0)
  # }
  
  # SET UP DATA FOR FITTING POISSON MODELS
  # --------------------------------------
  # Ensure that none of the topic proportions are exactly zero or
  # exactly one.
  L <- fit$L
  L <- pmax(L,control$minval)
  L <- pmin(L,1 - control$minval)
  
  # Add "pseudocounts" to the data, and get the Poisson NMF loadings
  # matrix. From this point on, we will fit Poisson glm models x ~
  # Poisson (u), u = sum(L*f), where x is a column of X.
  out <- fastTopics:::add_pseudocounts(X,s*L,pseudocount)
  X <- out$X
  L <- out$L
  
  # FIT POISSON MODELS
  # ------------------
  # For each column j of the counts matrix, compute MLEs of the
  # parameters in the Poisson glm, x ~ Poisson(u), in which the
  # Poisson rates are u = sum(L*f), and f = F[j,].
  if (verbose)
    cat(sprintf("Fitting %d Poisson models with k=%d using method=\"%s\".\n",
                m,k,fit.method))
  nc <- fastTopics:::initialize.multithreading(control$nc,verbose)
  F <- fastTopics:::fit_poisson_models(X,L,fit.method,control$eps,control$numiter,
                                       control$tol,control$nc)
  print(head(F))
  print(summary(F))
  F <- pmax(F,control$minval)
  dimnames(F) <- dimnames(fit$F)
  
  # calc f0 
  global_cutoff <- 10^find_kde_midpoint(log10(c(F)))
  cat("Use f0 calculated from refitted F:", global_cutoff, "\n")
  f0 <- c(rep(global_cutoff, times=ncol(X)))
  names(f0) <- rownames(fit$F)
  
  # COMPUTE LOG-FOLD CHANGE STATISTICS
  # ----------------------------------
  # Perform MCMC to simulate the posterior distribution of the LFC
  # statistics, then compute key posterior quantities from the
  # simulated Monte Carlo samples.
  if (verbose) {
    cat("Computing log-fold change statistics from ")
    cat(sprintf("%d Poisson models with k=%d.\n",m,k))
  }
  ns <- control$ns
  D <- matrix(rnorm(ns*k),ns,k)
  U <- matrix(runif(ns*k),ns,k)
  M <- matrix(sample(k,ns*k,replace = TRUE),ns,k) - 1
  ncb <- RhpcBLASctl::blas_get_num_procs()
  RhpcBLASctl::blas_set_num_threads(control$nc.blas)
  if (nc == 1)
    out <- fastTopics:::compute_lfc_stats(X,F,L,f0,D,U,M,lfc.stat,control$conf.level,
                                          control$rw,control$eps,verbose)
  else {
    out <- fastTopics:::compute_lfc_stats_multicore(X,F,L,f0,D,U,M,lfc.stat,
                                                    control$conf.level,control$rw,
                                                    control$eps,control$nc,control$nsplit,
                                                    verbose)
  }
  blas_set_num_threads(ncb)
  if (any(out$ar == 0))
    warning("One or more MCMC simulations yielded acceptance rates of zero; ",
            "consider increasing the number of Monte Carlo samples ",
            "(control$ns) or modifying the noise level of the random-walk ",
            "proposal distribution (control$rw) to improve the acceptance ",
            "rates")
  
  # STABILIZE ESTIMATES USING ADAPTIVE SHRINKAGE
  # --------------------------------------------
  # If requested, use adaptive shrinkage to stabilize the log-fold
  # change estimates. Here we need to carefully edge cases such as
  # se's of zero.
  if (shrink.method == "ash") {
    if (verbose)
      cat("Stabilizing posterior log-fold change estimates using adaptive",
          "shrinkage.\n")
    se <- with(out,postmean/z)
    se[out$z == 0] <- as.numeric(NA)
    se[out$postmean == 0] <- 0
    res          <- shrink_estimates(out$postmean,se,...)
    out$postmean <- res$b
    out$z        <- res$z
    out$lfsr     <- res$lfsr
    out$lpval    <- as.numeric(NA)
    out$svalue   <- res$svalue
    out$positive_prob <- res$positive_prob
    out$negative_prob <- res$negative_prob
    out$ash      <- res$ash
    dimnames(out$lfsr)   <- dimnames(F)
    dimnames(out$svalue) <- dimnames(F)
  } else {
    
    # Compute the -log10 two-tailed p-values computed from the z-scores.
    out$lpval  <- -fastTopics:::lpfromz(out$z)
    out$svalue <- as.numeric(NA)
    out$lfsr   <- as.numeric(NA)
  }
  
  # Return the Poisson model MLEs (F), the log-fold change statistics
  # (est, postmean, lower, upper, z, lpval) and local false sign
  # rates (lfsr), and the relative rates under the "null" model (f0).
  out$F  <- F
  out$f0 <- f0
  class(out) <- c("topic_model_de_analysis","list")
  return(out)
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

get_signac_umap <- function(count_matrix, L, cell_scores, top_peaks = 2000, seed = 24) {
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
