#' Genome Coverage of Topic Annotations in 0.01 cM Blocks
#'
#' For each topic annotation, the percentage of the genome's 0.01 centimorgan
#' (cM) blocks that contain at least one annotated SNP. This is the size metric
#' Tashman et al. (2021) use to decide when S-LDSC's block-jackknife p-values
#' control type I error.
#'
#' @section Why blocks and not SNPs:
#' S-LDSC's significance tests rest on the block jackknife. For small
#' annotations the jackknife standard error is noisy and correlated with the
#' coefficient estimate, and the estimate itself is right-skewed, so very small
#' p-values are inflated even when P = 0.05 is fine. The inflation depends on how
#' many independent pieces of the genome an annotation touches, not only on how
#' many SNPs it contains: 0.5% of SNPs packed into a few regions is far less
#' reliable than the same number spread across the genome. The percentage of
#' 0.01 cM blocks overlapped captures both size and clustering, and of six
#' metrics compared it separated well- from badly-calibrated annotations best.
#'
#' In null simulations (3,750 phenotypes, 246 annotations) type I error was
#' controlled for annotations covering at least 1.7% of 0.01 cM blocks for the
#' one-sided test of the coefficient and the enrichment test (0.83% suffices for
#' the enrichment test alone), and at least 8.3% for the two-sided test of the
#' coefficient. All baseline-model annotations pass.
#'
#' @section Definition:
#' Blocks are `floor(CM / block_cm)` within each chromosome, using the genetic
#' map positions (`CM` column) of the reference SNPs in the baseline annotation
#' files. The genome's blocks are those containing at least one reference SNP;
#' an annotation overlaps a block if any SNP in the block has a non-zero
#' annotation value. The topic annotations must be thin annotations whose rows
#' match the baseline files, as written by [run_sldsc()].
#'
#' @param ldsc_res_dir Directory holding `k*_output/annotations/<trait>/` for each
#'   topic, as produced by [run_sldsc()].
#' @param trait Trait name used in the annotation filenames.
#' @param nTopics Number of topics (K).
#' @param baseline_prefix Prefix of the baseline annotation files (with `CHR`,
#'   `BP`, `SNP`, `CM` columns), matching the `--ref-ld-chr` of the S-LDSC runs.
#' @param chrs Chromosomes to use. Defaults to 1:22.
#' @param block_cm Block size in cM. Defaults to 0.01, as in Tashman et al.
#' @return A data frame with one row per topic: `topic`, `n_snps` (reference
#'   SNPs in the annotation), `prop_snps` (as a fraction of all reference SNPs;
#'   note S-LDSC's `Prop._SNPs` counts MAF 5-50% SNPs only), `n_blocks` (blocks
#'   overlapped) and `pct_blocks` (percentage of the genome's blocks overlapped).
#'   The number of blocks in the genome is attached as attribute
#'   `n_blocks_genome`.
#' @references Tashman KC, Cui R, O'Connor LJ, Neale BM, Finucane HK (2021).
#'   Significance testing for small annotations in stratified LD-Score
#'   regression. medRxiv 2021.03.13.21249938.
#' @seealso [get_cs()]
#' @export
annot_block_coverage <- function(ldsc_res_dir, trait, nTopics, baseline_prefix,
                                 chrs = 1:22, block_cm = 0.01) {

  if (!requireNamespace("data.table", quietly = TRUE))
    stop("data.table is required for annot_block_coverage()")
  fread <- data.table::fread

  kdir <- function(k) file.path(ldsc_res_dir, paste0("k", k, "_output"))
  n_snps <- numeric(nTopics); n_blocks <- numeric(nTopics)
  n_ref <- 0; n_blocks_genome <- 0
  for (chr in chrs) {
    b <- fread(paste0(baseline_prefix, chr, ".annot.gz"), select = c("CM"),
               showProgress = FALSE)
    cm <- b$CM; rm(b)
    if (all(cm == 0))
      stop("chr", chr, ": the CM column of the baseline annotation is all zero; ",
           "0.01 cM block coverage needs genetic map positions")
    blk <- floor(cm / block_cm)
    n_ref <- n_ref + length(blk)
    n_blocks_genome <- n_blocks_genome + length(unique(blk))
    for (k in seq_len(nTopics)) {
      f <- fread(file.path(kdir(k), "annotations", trait,
                           paste0(trait, ".", chr, ".annot.gz")), showProgress = FALSE)
      a <- as.numeric(f[[1]]); rm(f)
      if (length(a) != length(blk))
        stop("chr", chr, ", topic ", k, ": annotation has ", length(a),
             " rows but the baseline has ", length(blk),
             "; the topic annotations must be thin annotations matching baseline_prefix")
      on <- a != 0
      n_snps[k]   <- n_snps[k] + sum(on)
      n_blocks[k] <- n_blocks[k] + length(unique(blk[on]))
    }
  }

  out <- data.frame(topic = paste0("k", seq_len(nTopics)), n_snps = n_snps,
                    prop_snps = n_snps / n_ref, n_blocks = n_blocks,
                    pct_blocks = 100 * n_blocks / n_blocks_genome)
  attr(out, "n_blocks_genome") <- n_blocks_genome
  out
}
