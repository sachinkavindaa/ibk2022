# ============================================================
# STEP 2: Filter low-abundance genes and low-quality samples
# Separate thresholds per domain (virus is sparser than bacteria)
# ============================================================

library(ggplot2)

dir.create("01_QC_filtered", showWarnings = FALSE)

metadata <- read.table("groups_filtered.tsv",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)
colnames(metadata) <- c("SampleID", "Group")

# ---- Filtering thresholds (adjust after inspecting step 1 output) ----
# min_count  = minimum count value to consider a gene "present" in a sample
# min_samples = gene must be present in at least this many samples

filter_params <- list(
  bacteria  = list(min_count = 2,  min_samples = 3),
  virus     = list(min_count = 1,  min_samples = 2),   # virus is sparse, be lenient
  eukaryote = list(min_count = 2,  min_samples = 3)
)

domains <- c("bacteria", "virus", "eukaryote")

for (domain in domains) {
  
  cat("\n=============================\n")
  cat("Filtering:", domain, "\n")
  cat("=============================\n")
  
  counts <- read.table(paste0(domain, "_count_matrix.tsv"),
                       header = TRUE,
                       row.names = 1,
                       sep = "\t",
                       check.names = FALSE)
  
  colnames(counts) <- gsub("_vs_.*", "", colnames(counts))
  colnames(counts) <- gsub(" Read Count", "", colnames(counts))
  
  # Match to metadata
  common_samples <- intersect(colnames(counts), metadata$SampleID)
  counts <- counts[, common_samples, drop = FALSE]
  meta   <- metadata[match(common_samples, metadata$SampleID), ]
  
  counts[is.na(counts)] <- 0
  counts <- round(as.matrix(counts))
  
  cat("Before filtering:", nrow(counts), "genes,", ncol(counts), "samples\n")
  
  # ---- Sample-level QC: flag low library size samples ----
  lib_sizes <- colSums(counts)
  low_cutoff <- quantile(lib_sizes, 0.10)   # bottom 10% flagged
  low_samples <- names(lib_sizes[lib_sizes < low_cutoff])
  
  cat("Library size cutoff (10th percentile):", round(low_cutoff), "\n")
  if (length(low_samples) > 0) {
    cat("WARNING - Low library size samples (consider removing):\n")
    print(lib_sizes[low_samples])
  } else {
    cat("No extreme low-library-size samples detected.\n")
  }
  
  # ---- Gene filtering ----
  min_count   <- filter_params[[domain]]$min_count
  min_samples <- filter_params[[domain]]$min_samples
  
  keep <- rowSums(counts >= min_count) >= min_samples
  counts_filt <- counts[keep, , drop = FALSE]
  
  cat("After gene filtering:", nrow(counts_filt), "genes retained\n")
  cat("Genes removed:", sum(!keep), "\n")
  
  # Save filtered matrix
  write.table(counts_filt,
              paste0(domain, "_filtered_count_matrix.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
  
  # ---- Compare before/after distributions ----
  before_df <- data.frame(
    logsum = log2(rowSums(counts) + 1),
    Stage = "Before filtering"
  )
  after_df <- data.frame(
    logsum = log2(rowSums(counts_filt) + 1),
    Stage = "After filtering"
  )
  compare_df <- rbind(before_df, after_df)
  
  p <- ggplot(compare_df, aes(x = logsum, fill = Stage)) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.6, color = "white") +
    theme_bw() +
    xlab("log2(gene total counts + 1)") +
    ylab("Number of genes") +
    ggtitle(paste(domain, "- Gene count distribution before vs after filtering")) +
    scale_fill_manual(values = c("Before filtering" = "gray60",
                                 "After filtering"  = "steelblue"))
  
  ggsave(paste0("01_QC_filtered/", domain, "_filter_comparison.pdf"),
         p, width = 8, height = 5)
  
  # Summary table
  summary_df <- data.frame(
    Domain         = domain,
    Genes_before   = nrow(counts),
    Genes_after    = nrow(counts_filt),
    Genes_removed  = sum(!keep),
    Samples        = ncol(counts_filt),
    Min_libsize    = min(lib_sizes),
    Max_libsize    = max(lib_sizes),
    Median_libsize = median(lib_sizes)
  )
  
  write.table(summary_df,
              paste0("01_QC_filtered/", domain, "_filter_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Saved: ", domain, "_filtered_count_matrix.tsv\n")
}

cat("\nDONE. Check 01_QC_filtered/ and review *_filter_summary.tsv files.\n")
cat("Then proceed to step3_normalize.R\n")
