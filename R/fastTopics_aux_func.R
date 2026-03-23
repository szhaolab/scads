#' Helper functions for scads
#' @importFrom utils modifyList
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom ashr ash
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom RhpcBLASctl blas_get_num_procs
#' @importFrom fastTopics de_analysis_control_default
#' @importFrom parallel splitIndices
#' @importFrom pbapply pboptions
#' @importFrom pbapply pblapply
#' 

# Use compute_lfc_stats function from this instead of fastTopics
compute_lfc_stats2 <- #LY
  function (X, F, L, f0,
            D = matrix(rnorm(1000*k),1000,ncol(F)),
            U = matrix(runif(1000*k),1000,ncol(F)),
            M = matrix(sample(k,1000*k,replace=TRUE),1000,k)-1,
            lfc.stat = "le", conf.level = 0.68, rw = 0.3, e = 1e-15,
            verbose = TRUE) {
    
    # Get the number of counts matrix columns (m) and the number of
    # topics (k).
    m <- nrow(F)
    k <- ncol(F)
    
    # Allocate storage for the outputs.
    est      <- matrix(0,m,k)
    postmean <- matrix(0,m,k)
    lower    <- matrix(0,m,k)
    upper    <- matrix(0,m,k)
    ar       <- matrix(0,m,k)
    dimnames(est)      <- dimnames(F)
    dimnames(postmean) <- dimnames(F)
    dimnames(lower)    <- dimnames(F)
    dimnames(upper)    <- dimnames(F)
    dimnames(ar)       <- dimnames(F)
    
    # Fill in the outputs, row by row. The core computation is performed
    # by compute_lfc_helper.
    ls <- colSums(L)
    if (verbose)
      pb <- progress_bar$new(total = m)
    for (j in 1:m) {
      if (verbose)
        pb$tick()
      out <- compute_lfc_stats_helper2(j,X,F,L,D,U,M,ls,f0,lfc.stat, #LY
                                       conf.level,rw,e) 
      est[j,]      <- out$dat["est",]
      postmean[j,] <- out$dat["postmean",]
      lower[j,]    <- out$dat["lower",]
      upper[j,]    <- out$dat["upper",]
      ar[j,]       <- out$ar
    }
    
    # Compute the z-scores, then output the LFC point estimates (est),
    # the posterior means (postmean), lower and upper limits of the HPD
    # intervals, the z-scores (z), and MCMC acceptance rates (ar).
    return(list(ar       = ar,
                est      = est/log(2),
                postmean = postmean/log(2),
                lower    = lower/log(2),
                upper    = upper/log(2),
                z        = fastTopics:::compute_zscores(postmean,lower,upper)))
  }


compute_lfc_vsnull2 <- function (f, f0, samples, conf.level) {
  
  est0     <- log(f0) #SZ
  est      <- log(f) - est0 #SZ
  # samples  <- samples - est0 #SZ
  samples <- sweep(samples, MARGIN = 2, STATS = as.numeric(est0), FUN = "-") #LY
  postmean <- colMeans(samples)
  ans      <- fastTopics:::compute_hpd_intervals(samples,conf.level)
  return(rbind(est      = est,
               postmean = postmean,
               lower    = ans$lower,
               upper    = ans$upper))
}


# Modified compute_lfc_stats_helper function from fastTopics
compute_lfc_stats_helper2 <- function (j, X, F, L, D, U, M, ls, f0, #LY
                                       lfc.stat, conf.level, rw, e) {
  k <- ncol(F)
  if (fastTopics:::is.sparse.matrix(X)) {
    dat <- fastTopics:::get.nonzeros(X,j)
    out <- fastTopics:::simulate_posterior_poisson_sparse_rcpp(dat$x,L[dat$i,],ls,
                                                               F[j,],D,U,M,rw,e)
  } else
    out <- fastTopics:::simulate_posterior_poisson_rcpp(X[,j],L,F[j,],D,U,M,rw,e)
  if (lfc.stat == "vsnull") {
    dat <- compute_lfc_vsnull2(F[j,], f0[j,], out$samples, conf.level) #SZ
  } else if (lfc.stat == "le") {
    dat <- fastTopics:::compute_lfc_le(F[j,],out$samples,conf.level)
  } else {
    dat <- fastTopics:::compute_lfc_pairwise(F[j,],lfc.stat,out$samples,conf.level)}
  return(list(ar = out$ar,dat = dat))
}

# Modified compute_lfc_stats_multicore function from fastTopics 
# Add this temporary debugging version to see what's actually failing:
compute_lfc_stats_multicore_debug <- function (X, F, L, f0, D, U, M, lfc.stat,
                                               conf.level, rw, e, nc, nsplit = 100,
                                               verbose = TRUE) {
  
  m <- nrow(F)
  k <- ncol(F)
  
  nsplit <- min(m, nsplit)
  cols   <- parallel::splitIndices(m, nsplit)
  dat    <- vector("list", nsplit)
  
  for (i in 1:nsplit) {
    j <- cols[[i]]
    dat[[i]] <- list(
      X = X[, j, drop = FALSE],
      F = F[j, , drop = FALSE],
      f0 = f0[j, , drop = FALSE]
    )
  }
  
  # Test function with error capturing
  parlapplyf <- function(dat, L, D, U, M, lfc.stat, conf.level, rw, e) {
    tryCatch({
      compute_lfc_stats2(dat$X, dat$F, L, dat$f0, D, U, M, lfc.stat, 
                         conf.level, rw, e, verbose = FALSE)
    }, error = function(err) {
      # Return the error so we can see what went wrong
      list(error = TRUE, message = as.character(err), 
           traceback = as.character(sys.calls()))
    })
  }
  
  if (verbose)
    op <- pbapply::pboptions(type = "txt", txt.width = 70)
  else
    op <- pbapply::pboptions(type = NULL)
  
  ans <- pbapply::pblapply(cl = nc, dat, parlapplyf, L, D, U, M, 
                           lfc.stat, conf.level, rw, e)
  
  pbapply::pboptions(op)
  
  # Check for errors in results
  for (i in 1:length(ans)) {
    if (is.list(ans[[i]]) && !is.null(ans[[i]]$error) && ans[[i]]$error) {
      cat("Error in worker", i, ":\n")
      cat("Message:", ans[[i]]$message, "\n")
      cat("Traceback:", paste(ans[[i]]$traceback, collapse = "\n"), "\n\n")
    }
  }
  
  return(ans)
}

compute_lfc_stats_multicore2 <- function (X, F, L, f0, D, U, M, lfc.stat, #LY
                                          conf.level, rw, e, nc, nsplit = 100,
                                          verbose = TRUE) {
  
  # Get the number of counts matrix columns (m) and the number of
  # topics (k).
  m <- nrow(F)
  k <- ncol(F)
  
  # Split the data.
  nsplit <- min(m,nsplit)
  cols   <- parallel::splitIndices(m,nsplit)
  dat    <- vector("list",nsplit)
  for (i in 1:nsplit) {
    j        <- cols[[i]]
    # dat[[i]] <- list(X = X[,j,drop = FALSE],F = F[j,,drop = FALSE],f0 = f0[j,]) #SZ 
    dat[[i]] <- list(
      X = X[, j, drop = FALSE], 
      F = F[j, , drop = FALSE], 
      f0 = f0[j, , drop = FALSE] 
    ) # LY (12/22/25)
  }
  
  # Distribute the calculations using pblapply.
  parlapplyf <- function (dat, L, D, U, M, lfc.stat, conf.level, rw, e)
    compute_lfc_stats2(dat$X,dat$F,L,dat$f0,D,U,M,lfc.stat,conf.level,rw,e, #LY
                       verbose = FALSE)
  if (verbose)
    op <- pbapply::pboptions(type = "txt",txt.width = 70)
  else
    op <- pbapply::pboptions(type = NULL)
  ans <- pbapply::pblapply(cl = nc,dat,parlapplyf,L,D,U,M,lfc.stat,conf.level,rw,e)
  pbapply::pboptions(op)
  
  # Combine the individual compute_lfc_stats outputs, and output the
  # combined result.
  out <- list(ar       = matrix(0,m,k),
              est      = matrix(0,m,k),
              postmean = matrix(0,m,k),
              lower    = matrix(0,m,k),
              upper    = matrix(0,m,k),
              z        = matrix(0,m,k))
  dimnames(out$ar)       <- dimnames(F)
  dimnames(out$est)      <- dimnames(F)
  dimnames(out$postmean) <- dimnames(F)
  dimnames(out$lower)    <- dimnames(F)
  dimnames(out$upper)    <- dimnames(F)
  dimnames(out$z)        <- dimnames(F)
  for (i in 1:nsplit) {
    j <- cols[[i]]
    out$ar[j,]       <- ans[[i]]$ar
    out$est[j,]      <- ans[[i]]$est
    out$postmean[j,] <- ans[[i]]$postmean
    out$lower[j,]    <- ans[[i]]$lower
    out$upper[j,]    <- ans[[i]]$upper
    out$z[j,]        <- ans[[i]]$z
  }
  return(out)
}

# Modified de_analysis function from fastTopics
de_analysis2 <- function (fit, X, s = rowSums(X), pseudocount = 0.01,
                          fit.method = c("scd","em","mu","ccd","glm"),
                          shrink.method = c("ash","none"), lfc.stat = "le",
                          control = list(), verbose = TRUE, f0 = NULL, ...) { #LY
  
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
  
  # Check f0 input # LY----------------------------------------
  if (!is.null(f0)) {
    if (length(f0) == 1) {
      # Case 1: user provided a single background value
      background <- f0
      cat("Use f0 provided as single background:", background, "\n")
      f0 <- rep(background, times = ncol(X)) # vector 
      # names(f0) <- colnames(X)
      f0 <- matrix(rep(f0, times = k), ncol = k, byrow = F) # matrix   
      cat("f0: ", dim(f0), "\n")
      rownames(f0) <- colnames(X)
    } else if (length(f0) == ncol(X)) {
      # Case 2: user provided a vector matching number of columns
      cat("Use f0 provided as vector of length", length(f0), "\n")
      names(f0) <- colnames(X)
    } else if (is.matrix(f0) || is.data.frame(f0)) {
      cat("Use f0 provided as a matrix", dim(f0), "\n", class(f0), "\n")
      f0 <- as.matrix(f0)
      rownames(f0) <- colnames(X)
    } else {
      cat("f0 provided is unacceptable, please check")
    }
  } else {
    cat("No f0 provided, using default behavior.\n")
  }
  # LY, end----------------------------------------
  
  # # calc f0 
  # global_cutoff <- 10^find_kde_midpoint(log10(c(F)))
  # cat("Use f0 calculated from refitted F:", global_cutoff, "\n")
  # f0 <- c(rep(global_cutoff, times=ncol(X)))
  # names(f0) <- rownames(fit$F)
  
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
  if (verbose) {
    cat("F matrix summary:\n")
    print(summary(F))
  }
  F <- pmax(F,control$minval)
  dimnames(F) <- dimnames(fit$F)
  
  
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
    out <- compute_lfc_stats2(X,F,L,f0,D,U,M,lfc.stat,control$conf.level, #LY
                              control$rw,control$eps,verbose)
  else {
    out <- compute_lfc_stats_multicore2(X,F,L,f0,D,U,M,lfc.stat, #LY
                                        control$conf.level,control$rw,
                                        control$eps,control$nc,control$nsplit,
                                        verbose)
  }
  RhpcBLASctl::blas_set_num_threads(ncb) # LY
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

