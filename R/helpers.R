#' Parse Peak Range Names into Components
#'
#' Parses genomic peak range strings in various formats into their components
#' (chromosome, start, end). Supports underscore (chr1_100_200), colon-dash
#' (chr1:100-200), and dash (chr1-100-200) formats.
#'
#' @param peak_ranges A character vector of peak range strings.
#' @return A data.frame with columns: seqnames, start, end.
#' @examples
#' parse_peak_ranges(c("chr1_100_200", "chr2_300_400"))
#' parse_peak_ranges(c("chr1:100-200", "chr2:300-400"))
#' @export
parse_peak_ranges <- function(peak_ranges) {
  if (is.null(peak_ranges) || length(peak_ranges) == 0) {
    stop("peak_ranges must be a non-empty character vector")
  }

  if (all(grepl("^[^_]+_[0-9]+_[0-9]+$", peak_ranges))) {
    # underscore format: chr1_10234_10734
    parts <- strsplit(peak_ranges, "_", fixed = TRUE)
    seqnames <- vapply(parts, `[`, 1, FUN.VALUE = "")
    start    <- as.integer(vapply(parts, `[`, 2, FUN.VALUE = ""))
    end      <- as.integer(vapply(parts, `[`, 3, FUN.VALUE = ""))

  } else if (all(grepl("^[^:]+:[0-9]+-[0-9]+$", peak_ranges))) {
    # colon-dash format: chr1:10234-10734
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
    stop("peak_ranges must be one of: 'chr_start_end', 'chr:start-end', or 'chr-start-end'")
  }

  # Ensure chr prefix
  seqnames <- ifelse(grepl("^chr", seqnames), seqnames, paste0("chr", seqnames))

  data.frame(
    seqnames = seqnames,
    start    = start,
    end      = end,
    stringsAsFactors = FALSE
  )
}


#' Compute Topic P-values and Binarize Peak Assignments
#'
#' Takes z-scores from differential expression analysis, computes one-sided
#' p-values, applies FDR correction, and binarizes based on a cutoff.
#' This is the core logic extracted from \code{run_fastTopics()}.
#'
#' @param z_matrix A matrix of z-scores (peaks x topics) from DE analysis.
#' @param fdr_cutoff FDR threshold for significance (default: 0.05).
#' @param total_bins Total number of hypothetical bins for FDR correction
#'   (default: ceiling(3e9 / 500), i.e., genome-wide bins).
#' @return A list with:
#'   \item{pvals}{One-sided p-value matrix (peaks x topics)}
#'   \item{qvals}{FDR-adjusted p-value matrix}
#'   \item{p_jk}{Binary matrix (1 = significant peak-topic pair)}
#'   \item{n_sig_per_topic}{Named vector of significant peaks per topic}
#' @export
compute_topic_pvalues <- function(z_matrix,
                                  fdr_cutoff = 0.05,
                                  total_bins = ceiling(3e9 / 500)) {

  if (!is.matrix(z_matrix)) stop("z_matrix must be a matrix")
  if (fdr_cutoff <= 0 || fdr_cutoff >= 1) stop("fdr_cutoff must be between 0 and 1")

  # One-sided p-values (testing for enrichment, i.e., positive z)
  pvals <- 1 - pnorm(z_matrix)

  # FDR correction per topic, accounting for genome-wide testing

  qvals <- apply(pvals, 2, function(col) {
    p.adjust(col, method = "fdr", n = total_bins)
  })
  qvals[is.na(qvals)] <- 1

  # Binarize
  p_jk <- ifelse(qvals < fdr_cutoff, 1, 0)
  p_jk[p_jk < 0] <- 0

  n_sig <- colSums(p_jk)

  list(
    pvals           = pvals,
    qvals           = qvals,
    p_jk            = p_jk,
    n_sig_per_topic = n_sig
  )
}


#' Build S-LDSC Shell Commands
#'
#' Constructs the shell commands for running S-LDSC pipeline steps.
#' All paths are properly quoted with \code{shQuote()} for shell safety.
#'
#' @param step Character: one of "munge", "annot", "l2", "h2".
#' @param ... Named arguments specific to each step (see Details).
#' @return A character string containing the shell command.
#' @details
#' For \code{step = "munge"}: requires polyfun_path, sumstats_gz, n, munged_out.
#' For \code{step = "annot"}: requires ldsc_path, user_bed, bim_file, out_annot.
#' For \code{step = "l2"}: requires polyfun_path, bfile_chr, hm3_snps, out_annot, out_prefix.
#' For \code{step = "h2"}: requires polyfun_path, munged_out, ref_ld, frqfile, w_ld, final_out, use_chr (logical).
#' @export
build_sldsc_command <- function(step, ...) {
  args <- list(...)

  cmd <- switch(step,
    "munge" = paste(
      "python3", shQuote(file.path(args$polyfun_path, "munge_polyfun_sumstats.py")),
      "--sumstats", shQuote(args$sumstats_gz),
      "--n", args$n,
      "--out", shQuote(args$munged_out)
    ),
    "annot" = paste(
      "python3", shQuote(file.path(args$ldsc_path, "make_annot.py")),
      "--bed-file", shQuote(args$user_bed),
      "--bimfile", shQuote(args$bim_file),
      "--annot-file", shQuote(args$out_annot)
    ),
    "l2" = paste(
      "python3", shQuote(file.path(args$polyfun_path, "ldsc.py")),
      "--l2",
      "--bfile", shQuote(args$bfile_chr),
      "--print-snps", shQuote(args$hm3_snps),
      "--ld-wind-cm 1",
      "--annot", shQuote(args$out_annot),
      "--thin-annot",
      "--out", shQuote(args$out_prefix)
    ),
    "h2" = {
      if (args$use_chr) {
        paste(
          "python3", shQuote(file.path(args$polyfun_path, "ldsc.py")),
          "--h2", shQuote(args$munged_out),
          "--ref-ld-chr", paste0(shQuote(args$ref_ld), ",", shQuote(args$baseline_dir)),
          "--frqfile-chr", shQuote(args$frqfile),
          "--w-ld-chr", shQuote(args$w_ld),
          "--overlap-annot",
          "--print-coefficients",
          "--print-delete-vals",
          "--out", shQuote(args$final_out)
        )
      } else {
        paste(
          "python3", shQuote(file.path(args$polyfun_path, "ldsc.py")),
          "--h2", shQuote(args$munged_out),
          "--ref-ld", paste0(shQuote(args$ref_ld), ",", shQuote(args$baseline_dir)),
          "--frqfile", shQuote(args$frqfile),
          "--w-ld", shQuote(args$w_ld),
          "--overlap-annot",
          "--print-coefficients",
          "--print-delete-vals",
          "--out", shQuote(args$final_out)
        )
      }
    },
    stop("Unknown step: ", step)
  )

  cmd
}


#' Log a message with timestamp
#'
#' Prints a timestamped message to the console. Used throughout the scads
#' pipeline for progress tracking.
#'
#' @param ... Message components (passed to \code{cat}).
#' @param verbose Logical; if FALSE, message is suppressed. Default TRUE.
#' @export
scads_log <- function(..., verbose = TRUE) {
  if (verbose) {
    cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n")
  }
}
