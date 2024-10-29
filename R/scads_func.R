#' Helper functions for scads
#' 
#' modified de_analysis function from fastTopics
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
  background = 100/(3E9/500)
  f0 <- c(rep(background, times=ncol(X))) # colSums(X)/sum(s)
  # f0 <- c(rep(0, times=ncol(count_matrix)))
  # f0 <- ifelse(peak_assignments==1, 9.268e-05, 
  #              ifelse(peak_assignments==2, 9.722e-05, 
  #                     ifelse(peak_assignments==3, 9.292e-05, 
  #                            ifelse(peak_assignments==4, 9.648e-05, 9.263e-05))))
  names(f0) <- rownames(fit$F)
  
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
  ncb <- blas_get_num_procs()
  blas_set_num_threads(control$nc.blas)
  if (nc == 1)
    out <- fastTopics:::compute_lfc_stats(X,F,L,f0,D,U,M,lfc.stat,control$conf.level,
                             control$rw,control$eps,verbose)
  else {
    out <- compute_lfc_stats_multicore(X,F,L,f0,D,U,M,lfc.stat,
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