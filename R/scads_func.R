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

