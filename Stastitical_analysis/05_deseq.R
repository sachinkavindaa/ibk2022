# ============================================================
# STEP 5: DESeq2 differential expression
# Input:  *_filtered_count_matrix.tsv  (RAW counts from step2)
#         groups_filtered.tsv
# IMPORTANT: Always feed RAW counts into DESeq2, never normalized!
# Output: 05_DESeq2/
# ============================================================

library(DESeq2)
library(ggplot2)
library(ggrepel)

dir.create("05_DESeq2", showWarnings = FALSE)

metadata <- read.table("groups_filtered.tsv",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)
colnames(metadata) <- c("SampleID", "Group")
metadata$Group <- factor(metadata$Group, levels = c("healthy", "infected"))

domains <- c("bacteria", "virus", "eukaryote")

for (domain in domains) {
  
  cat("\n=============================\n")
  cat("DESeq2:", domain, "\n")
  cat("=============================\n")
  
  # ---- Load RAW filtered counts from step2 ----
  counts <- read.table(paste0(domain, "_filtered_count_matrix.tsv"),
                       header = TRUE,
                       row.names = 1,
                       sep = "\t",
                       check.names = FALSE)
  
  counts <- as.matrix(counts)
  counts[is.na(counts)] <- 0
  counts <- round(counts)
  
  # Match samples to metadata
  common_samples <- intersect(colnames(counts), metadata$SampleID)
  counts <- counts[, common_samples, drop = FALSE]
  meta   <- metadata[match(common_samples, metadata$SampleID), ]
  rownames(meta) <- meta$SampleID
  
  cat("Samples:", ncol(counts), "| Genes:", nrow(counts), "\n")
  
  if (nrow(counts) < 10) {
    cat("Too few genes for DESeq2. Skipping", domain, "\n")
    next
  }
  
  # ==========================
  # A) Run DESeq2
  # ==========================
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = meta,
    design    = ~ Group
  )
  
  # Pre-filter: remove genes with very low total counts
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds  <- dds[keep, ]
  cat("Genes after DESeq2 pre-filter (>=10 counts in >=3 samples):",
      nrow(dds), "\n")
  
  dds <- DESeq(dds)
  
  # ==========================
  # B) Extract results
  # contrast = infected vs healthy (infected is numerator)
  # ==========================
  
  res <- results(dds,
                 contrast  = c("Group", "infected", "healthy"),
                 alpha     = 0.05)
  
  summary(res)
  
  # -------------------------------------------------------
  # NEW: LFC shrinkage using apeglm
  # This pulls noisy fold-change estimates (from low-count
  # genes) toward zero, making results more reliable.
  # It does NOT change which genes are significant.
  # The coef name must match your contrast exactly.
  # -------------------------------------------------------
  res <- lfcShrink(dds,
                   coef = "Group_infected_vs_healthy",
                   type = "apeglm")
  
  res_df <- as.data.frame(res)
  res_df$GeneID <- rownames(res_df)
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
  
  # Save full results
  write.table(res_df,
              paste0("05_DESeq2/", domain, "_all_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Significant genes: padj < 0.05 AND |log2FC| >= 1
  sig <- subset(res_df,
                !is.na(padj) &
                  padj < 0.05 &
                  abs(log2FoldChange) >= 1)
  
  sig_up   <- subset(sig, log2FoldChange > 0)
  sig_down <- subset(sig, log2FoldChange < 0)
  
  cat("Significant genes (padj<0.05, |LFC|>=1):", nrow(sig), "\n")
  cat("  Upregulated in infected:  ", nrow(sig_up), "\n")
  cat("  Downregulated in infected:", nrow(sig_down), "\n")
  
  write.table(sig,
              paste0("05_DESeq2/", domain, "_significant_genes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(sig_up,
              paste0("05_DESeq2/", domain, "_upregulated.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(sig_down,
              paste0("05_DESeq2/", domain, "_downregulated.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # ==========================
  # C) Volcano plot
  # ==========================
  
  res_df$Significance <- "Not significant"
  res_df$Significance[!is.na(res_df$padj) &
                        res_df$padj < 0.05 &
                        res_df$log2FoldChange >= 1]  <- "Up in infected"
  res_df$Significance[!is.na(res_df$padj) &
                        res_df$padj < 0.05 &
                        res_df$log2FoldChange <= -1] <- "Down in infected"
  
  # Label top 15 significant genes
  top_genes <- head(sig[order(sig$padj), ], 15)
  
  p_volcano <- ggplot(res_df,
                      aes(x = log2FoldChange,
                          y = -log10(pvalue),
                          color = Significance)) +
    geom_point(size = 1.2, alpha = 0.6) +
    geom_text_repel(data = top_genes,
                    aes(label = GeneID),
                    size = 2.5,
                    max.overlaps = 20,
                    color = "black") +
    scale_color_manual(values = c(
      "Not significant"  = "grey70",
      "Up in infected"   = "#F44336",
      "Down in infected" = "#2196F3"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme_bw() +
    xlab("log2 Fold Change (infected / healthy)") +
    ylab("-log10(p-value)") +
    ggtitle(paste(domain, "- Volcano plot"))
  
  ggsave(paste0("05_DESeq2/", domain, "_volcano.pdf"),
         p_volcano, width = 8, height = 6)
  
  # ==========================
  # D) VST PCA
  # ==========================
  
  # blind = TRUE is recommended for QC/exploratory PCA
  # (does not use group info, so gives an unbiased view)
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  
  pca_data <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
  var_exp  <- round(100 * attr(pca_data, "percentVar"), 1)
  
  p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                color = Group, label = name)) +
    geom_point(size = 3.5) +
    geom_text_repel(size = 2.5, max.overlaps = 20) +
    scale_color_manual(values = c("healthy" = "#2196F3",
                                  "infected" = "#F44336")) +
    theme_bw() +
    xlab(paste0("PC1 (", var_exp[1], "%)")) +
    ylab(paste0("PC2 (", var_exp[2], "%)")) +
    ggtitle(paste(domain, "- DESeq2 VST PCA"))
  
  ggsave(paste0("05_DESeq2/", domain, "_VST_PCA.pdf"),
         p_pca, width = 7, height = 6)
  
  # ==========================
  # E) MA plot
  # (now uses shrunken LFC — cleaner, less fan shape)
  # ==========================
  
  pdf(paste0("05_DESeq2/", domain, "_MA_plot.pdf"), width = 6, height = 5)
  plotMA(res, main = paste(domain, "MA plot (shrunken LFC)"), ylim = c(-5, 5))
  dev.off()
  
  cat("Done:", domain, "\n")
}

cat("\n\nDONE step5. Check 05_DESeq2/ folder.\n")
cat("Key outputs:\n")
cat("  *_all_results.tsv        - every gene with LFC, pvalue, padj\n")
cat("  *_significant_genes.tsv  - padj<0.05, |LFC|>=1\n")
cat("  *_upregulated.tsv        - up in infected\n")
cat("  *_downregulated.tsv      - down in infected\n")
cat("  *_volcano.pdf            - volcano plot\n")
cat("  *_VST_PCA.pdf            - PCA on VST counts\n")
cat("  *_MA_plot.pdf            - MA plot (shrunken LFC)\n")