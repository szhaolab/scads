#' Compute Single-Cell Level Scores (c_i)
#'
#' This function computes the single-cell level scores (\eqn{c_i}) based on the fastTopics model outputs
#' and other required parameters.
#'
#' @param topic_res Results from run_fastTopics
#'  - Pmat: Matrix of probabilities of peaks in topics (peaks x topics), obtained from the F matrix of fastTopics.
#'  - Lmat: Matrix of topic proportions for each cell (cells x topics), obtained from the L matrix of fastTopics.
#' @param ldsc_res_dir Directory containing folders for each topic with enrichment results from run_ldsc
#'  - tau_Ck: Vector of per-topic heritability contributions (\eqn{\tau_{C_k}}) of length K (number of topics).
#' @param trait Trait name used in LDSC results files.
#' @param nTopics Number of topics in the model.
#' @return A list containing:
#'   - cs: A vector of single-cell level scores (\eqn{c_i}) of length I (number of cells).
#'   - ldsc_res_table: A data frame of LDSC results for each topic.
#' @details
#' The function computes the scores based on the corrected approximation.
#' @export
get_cs <- function(topic_res, ldsc_res_dir, trait, nTopics) {

  # Extract matrices from topic_res
  p_jk <- topic_res$Pmat  # J x K matrix
  l_ik <- topic_res$Lmat  # I x K matrix

  # Initialize an empty list to store per-snp heritability (tau) values
  tau <- list()
  tau_display <- list()

  # Loop over each topic
  for (k in 1:nTopics) {
    # Construct the folder name for the current topic
    topic_folder <- paste0("k", k, "_output")

    # Construct the full path to the .results file
    results_file <- file.path(ldsc_res_dir, topic_folder, "results", sprintf("%s_enrichment.results", trait))

    # Check if the file exists
    if (file.exists(results_file)) {
      # Read the .results file
      res <- read.table(
        results_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      # Ensure that the data has at least one row
      if (nrow(res) >= 1) {

        # Extract the tau value (per-snp h2 for annot) for the topic annotation
        # Store the enrichment value in the tau list
        tau[[k]] <- res$Prop._h2[1]

        # Store the LDSC results for display
        if (is.null(tau_display)) {
          tau_display <- res[1, ]
        } else {
          tau_display <- rbind(tau_display, res[1, ])
        }

      } else {
        warning(paste("File", results_file, "does not have at least one row. Skipping."))
        tau[[k]] <- NA  # Assign NA if there are not enough rows
      }
    } else {
      warning(paste("File", results_file, "does not exist. Skipping."))
      tau[[k]] <- NA  # Assign NA if the file does not exist
    }
  }

  # Convert tau to a numeric vector
  tau_Ck <- unlist(tau)

  # Print the tau vector
  print(tau_display)

  # Number of peaks (J), topics (K), and cells (I)
  J <- nrow(p_jk)
  K <- ncol(p_jk)
  I <- nrow(l_ik)

  # Check dimensions
  if (ncol(l_ik) != K) {
    stop("Number of topics in l_ik and p_jk must be the same.")
  }
  if (length(tau_Ck) != K) {
    stop("Length of tau_Ck must be equal to the number of topics (K).")
  }

  # Since v_j is the number of variants and each peak corresponds to one variant,
  # we set v_j to a vector of ones.
  v_j <- rep(1, J)

  # get enrichment per topic
  e_Ck <- c(tau_display$Enrichment)

  # Step 1: Compute variants in open peaks
  v_j_p_jk <- p_jk * v_j  # Matrix of m by k

  # Step 2: Weight enrichment by variants in open peaks
  v_j_p_jk_e_Ck <- e_Ck * t(v_j_p_jk) # Matrix of k by m

  # Step 3: Compute E(M_i)
  E_Mi <- rowSums(l_ik %*% v_j_p_jk_e_Ck)  # Vector of length I

  # Step 4: Compute E(N_i)
  E_Ni <- rowSums(l_ik %*% t(v_j_p_jk))  # Vector of length I

  # Step 5:
  S_k_Cov1 <- colSums(p_jk * tau_Ck * v_j)
  First_Term_Cov_Mi_Ni <- l_ik %*% S_k_Cov1
  S_kk_Cov2 <- t(p_jk * tau_Ck * v_j) %*% p_jk
  Second_Term_Cov_Mi_Ni <- rowSums((l_ik %*% S_kk_Cov2) * l_ik)
  Cov_Mi_Ni <- First_Term_Cov_Mi_Ni - Second_Term_Cov_Mi_Ni

  # p_ij <- l_ik %*% t(p_jk)
  # Cov_Mi_Ni <- h_j * v_j * p_ij * (1 - p_ij)

  # Step 6:
  S_k_N <- colSums(p_jk * (v_j^2))
  First_Term_Var_Ni <- l_ik %*% S_k_N
  S_kk_Var_Ni <- t(p_jk * (v_j^2)) %*% p_jk
  Second_Term_Var_Ni <- rowSums((l_ik %*% S_kk_Var_Ni) * l_ik)
  Var_Ni <- First_Term_Var_Ni - Second_Term_Var_Ni

  # Step 7:
  cs <- (E_Mi / E_Ni) * (1 - (Cov_Mi_Ni / (E_Mi * E_Ni)) + (Var_Ni / (E_Ni^2)))

  # Return the single-cell scores and LDSC results table
  return(list(cs = cs,
              ldsc_res_table = tau_display,
              cs_dat = list(J = J, K =K, I=I, v_j=v_j,
                         E_Mi=E_Mi, E_Ni=E_Ni, Cov_Mi_Ni=Cov_Mi_Ni, Var_Ni=Var_Ni)))
}





#' Compute Single-Cell Level Scores (c_i) and Adjust for Small Annotations
#'
#' This function computes the single-cell level scores (c_i) based on the fastTopics model outputs
#' and other required parameters. Additionally, it extracts peak names, computes peak coordinates,
#' determines the total region covered by peaks per topic, and replaces enrichment values when the
#' total size is <3% of the genome.
#'
#' @param topic_res Results from run_fastTopics 
#'  - Pmat: Matrix of probabilities of peaks in topics (peaks x topics), obtained from the F matrix of fastTopics.
#'  - Lmat: Matrix of topic proportions for each cell (cells x topics), obtained from the L matrix of fastTopics.
#' @param ldsc_res_dir Directory containing folders for each topic with enrichment results from run_ldsc 
#'  - tau_Ck: Vector of per-topic heritability contributions (\eqn{\tau_{C_k}}) of length K (number of topics).
#'  - e_Ck: Vector of per-topic enrichment values of length K (number of topics).
#' @param trait Trait name used in LDSC results files.
#' @param nTopics Number of topics in the model.
#' @return A list containing:
#'   - cs: A vector of single-cell level scores (\eqn{c_i}) of length I (number of cells).
#'   - ldsc_res_table: A data frame of LDSC results for each topic.
#'   - peak_regions: A list of peak regions covered per topic.
#' @details
#' The function computes the scores based on the corrected approximation and adjusts for small annotations.
#' @export
get_cs <- function(topic_res, ldsc_res_dir, trait, nTopics) {
  
  # Extract matrices from topic_res
  p_jk <- as.matrix((topic_res$Pmat > 0.5) + 0)
  l_ik <- topic_res$Lmat  # I x K matrix
  
  # Extract peak names and compute coordinates
  peak_names <- rownames(topic_res$Fmat)
  peak_ranges <- do.call(rbind, strsplit(peak_names, "[:-]"))
  colnames(peak_ranges) <- c("chr", "start", "end")
  peak_ranges <- as.data.frame(peak_ranges, stringsAsFactors = FALSE)
  peak_ranges$start <- as.numeric(peak_ranges$start)
  peak_ranges$end <- as.numeric(peak_ranges$end)
  peak_ranges$length <- peak_ranges$end - peak_ranges$start
  
  # Initialize lists
  tau <- list()
  e_Ck <- list()
  tau_display <- list()
  peak_regions <- list()
  total_annot_region <- list()
  
  # Genome size threshold (3% of 3B = 90M base pairs)
  genome_threshold <- 3E9 * 0.0003
  
  # Loop over each topic
  for (k in 1:nTopics) {
    # Get peaks assigned to the topic
    topic_peaks <- peak_ranges[p_jk[, k] > 0, ]
    total_region_size <- sum(topic_peaks$length, na.rm = TRUE)
    total_annot_region[[k]] <- total_region_size
    
    # Store peak regions per topic
    peak_regions[[paste0("Topic_", k)]] <- topic_peaks
    
    # Construct file path to enrichment results
    topic_folder <- paste0("k", k, "_output")
    results_file <- file.path(ldsc_res_dir, topic_folder, "results", sprintf("%s_enrichment.results", trait))
    
    # Default enrichment value
    tau_k <- NA
    e_Ck_k <- NA
    
    # Read enrichment results file if it exists
    if (file.exists(results_file)) {
      res <- read.table(
        results_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      if (nrow(res) >= 1) {
        tau_k <- res$Prop._h2[1]  # Extract per-topic heritability contributions
        e_Ck_k <- res$Enrichment[1]  # Extract enrichment for each topic
        tau_display <- rbind(tau_display, res[1, ])
      }
    }
    
    tau[[k]] <- tau_k
    e_Ck[[k]] <- e_Ck_k
  }
  
  # Loop over each topic
  for (k in 1:nTopics) {
    # Adjust for small annotation size by replacing with the smallest reliable enrichment value
    if (!is.na(e_Ck[[k]]) && total_annot_region[[k]] < genome_threshold) {
      message(sprintf("Topic %d has small annotation size (%d bp). Applying baseline enrichment.", k, total_region_size))
      non_negative_enrichments <- unlist(e_Ck)[!is.na(unlist(e_Ck)) & unlist(e_Ck) >= 0]
      e_Ck[[k]] <- min(non_negative_enrichments, na.rm = TRUE)
    }
  }
  
  # Convert lists to numeric vectors
  tau_Ck <- unlist(tau)
  e_Ck <- unlist(e_Ck)
  cat("Enrichment for cell score calculation: ", e_Ck)
  
  # Compute single-cell scores
  v_j <- rep(1, nrow(p_jk))
  
  # Compute enrichment-weighted variant presence in peaks
  v_j_p_jk <- p_jk * v_j  
  v_j_p_jk_e_Ck <- e_Ck * t(v_j_p_jk) 
  
  # Compute expected number of variants in a cell
  E_Mi <- rowSums(l_ik %*% v_j_p_jk_e_Ck)  
  E_Ni <- rowSums(l_ik %*% t(v_j_p_jk))  
  
  # Compute covariance between expected number of variants and peak presence
  S_k_Cov1 <- colSums(p_jk * tau_Ck * v_j)
  First_Term_Cov_Mi_Ni <- l_ik %*% S_k_Cov1
  S_kk_Cov2 <- t(p_jk * tau_Ck * v_j) %*% p_jk
  Second_Term_Cov_Mi_Ni <- rowSums((l_ik %*% S_kk_Cov2) * l_ik)
  Cov_Mi_Ni <- First_Term_Cov_Mi_Ni - Second_Term_Cov_Mi_Ni
  
  # Compute variance of the expected number of variants per cell
  S_k_N <- colSums(p_jk * (v_j^2))
  First_Term_Var_Ni <- l_ik %*% S_k_N
  S_kk_Var_Ni <- t(p_jk * (v_j^2)) %*% p_jk
  Second_Term_Var_Ni <- rowSums((l_ik %*% S_kk_Var_Ni) * l_ik)
  Var_Ni <- First_Term_Var_Ni - Second_Term_Var_Ni
  
  # Compute single-cell scores using the final corrected formula
  cs <- (E_Mi / E_Ni) * (1 - (Cov_Mi_Ni / (E_Mi * E_Ni)) + (Var_Ni / (E_Ni^2)))
  
  return(list(
    cs = cs,
    ldsc_res_table = tau_display,
    peak_regions = peak_regions,
    cs_dat = list(J = nrow(p_jk), K = ncol(p_jk), I = nrow(l_ik), v_j = v_j,
                  E_Mi = E_Mi, E_Ni = E_Ni, Cov_Mi_Ni = Cov_Mi_Ni, Var_Ni = Var_Ni)
  ))
}






# get_cs <- function(topic_res, ldsc_res_dir, trait, nTopics) {
#   
#   # Extract matrices from topic_res
#   # p_jk <- topic_res$Pmat  # J x K matrix
#   p_jk <- as.matrix((topic_res$Pmat > 0.5) + 0)
#   l_ik <- topic_res$Lmat  # I x K matrix
#   
#   # Extract peak names and compute peak coordinates
#   peak_names <- rownames(topic_res$Fmat)
#   library(stringr)
#   peak_coords <- data.frame(
#     chr = str_extract(peak_names, "^chr[0-9XYM]+"),
#     start = as.numeric(str_extract(peak_names, "(?<=:)[0-9]+(?=-)")),
#     end = as.numeric(str_extract(peak_names, "(?<=-)[0-9]+"))
#   )
#   
#   # Read GWAS summary statistics
#   library(data.table)
#   gwas_file <- "/dartfs/rc/lab/S/Szhao/liyang/SCSviaTM/analysis/pbmc/scads_20241119/IBD_sumstats.txt.gz"
#   gwas_data <- fread(cmd = paste("zcat", gwas_file), header = TRUE)
#   gwas_data[, chr := paste0("chr", as.character(chr))]
#   gwas_data[, pos := as.numeric(pos)]
#   
#   # Create GRanges objects
#   library(GenomicRanges)
#   peaks_gr <- GRanges(
#     seqnames = peak_coords$chr,
#     ranges = IRanges(start = peak_coords$start, end = peak_coords$end)
#   )
#   variants_gr <- GRanges(
#     seqnames = gwas_data$chr,
#     ranges = IRanges(start = gwas_data$pos, end = gwas_data$pos)
#   )
#   
#   # Find overlaps and compute v_j
#   overlaps <- findOverlaps(peaks_gr, variants_gr)
#   variant_counts <- as.data.frame(table(queryHits(overlaps)))
#   v_j <- rep(0, length(peaks_gr))
#   v_j[as.numeric(variant_counts$Var1)] <- variant_counts$Freq
#   
#   # Ensure v_j matches p_jk
#   J <- nrow(p_jk)
#   if (length(v_j) != J) {
#     stop("Length of v_j does not match the number of peaks (J).")
#   }
#   
#   # Initialize an empty list to store enrichment values
#   tau <- list()
#   tau_display <- list()
#   
#   # Loop over each topic
#   for (k in 1:nTopics) {
#     # Construct the folder name for the current topic
#     topic_folder <- paste0("k", k, "_output")
#     
#     # Construct the full path to the .results file
#     results_file <- file.path(ldsc_res_dir, topic_folder, "results", sprintf("%s_enrichment.results", trait))
#     
#     # Check if the file exists
#     if (file.exists(results_file)) {
#       # Read the .results file
#       res <- read.table(
#         results_file,
#         header = TRUE,
#         sep = "\t",
#         stringsAsFactors = FALSE,
#         check.names = FALSE
#       )
#       
#       # Ensure that the data has at least one row
#       if (nrow(res) >= 1) {
#         
#         # Extract the tau value (per-snp h2 for annot) for the topic annotation
#         tau_value <- res$Prop._h2[1]
#         
#         # Store the enrichment value in the tau list
#         tau[[k]] <- tau_value
#         
#         # Store the LDSC results for display
#         if (is.null(tau_display)) {
#           tau_display <- res[1, ]
#         } else {
#           tau_display <- rbind(tau_display, res[1, ])
#         }
#         
#       } else {
#         warning(paste("File", results_file, "does not have at least one row. Skipping."))
#         tau[[k]] <- NA  # Assign NA if there are not enough rows
#       }
#     } else {
#       warning(paste("File", results_file, "does not exist. Skipping."))
#       tau[[k]] <- NA  # Assign NA if the file does not exist
#     }
#   }
#   
#   # Convert tau to a numeric vector
#   tau_Ck <- unlist(tau)
#   
#   # Print the tau vector
#   print(tau_display)
#   
#   # Number of topics (K) and cells (I)
#   K <- ncol(p_jk)
#   I <- nrow(l_ik)
#   
#   # get enrichment per topic
#   e_Ck <- c(tau_display$Enrichment)
#   
#   # Step 1: Compute variants in open peaks
#   v_j_p_jk <- p_jk * v_j  # Matrix of m by k
#   
#   # Step 2: Weight enrichment by variants in open peaks
#   v_j_p_jk_e_Ck <- e_Ck * t(v_j_p_jk) # Matrix of k by m
#   
#   # Step 3: Compute E(M_i)
#   E_Mi <- rowSums(l_ik %*% v_j_p_jk_e_Ck)  # Vector of length I
#   
#   # Step 4: Compute E(N_i)
#   E_Ni <- rowSums(l_ik %*% t(v_j_p_jk))  # Vector of length I
#   
#   # Step 5: 
#   h_j <- v_j * (p_jk %*% tau_Ck)
#   W_j <- c(h_j) * v_j
#   S_k_Cov1 <- colSums(p_jk * W_j)
#   First_Term_Cov_Mi_Ni <- l_ik %*% S_k_Cov1
#   S_kk_Cov2 <- t(p_jk * W_j) %*% p_jk
#   Second_Term_Cov_Mi_Ni <- rowSums((l_ik %*% S_kk_Cov2) * l_ik)
#   Cov_Mi_Ni <- First_Term_Cov_Mi_Ni - Second_Term_Cov_Mi_Ni
#   
#   # p_ij <- l_ik %*% t(p_jk)
#   # Cov_Mi_Ni <- h_j * v_j * p_ij * (1 - p_ij)
#   
#   # Step 6:
#   S_k_N <- colSums(p_jk * (v_j^2))
#   First_Term_Var_Ni <- l_ik %*% S_k_N
#   S_kk_Var_Ni <- t(p_jk * (v_j^2)) %*% p_jk
#   Second_Term_Var_Ni <- rowSums((l_ik %*% S_kk_Var_Ni) * l_ik)
#   Var_Ni <- First_Term_Var_Ni - Second_Term_Var_Ni
#   
#   # Step 7:
#   cs <- (E_Mi / E_Ni) * (1 - (Cov_Mi_Ni / (E_Mi * E_Ni)) + (Var_Ni / (E_Ni^2)))
# }
