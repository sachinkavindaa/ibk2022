# ============================================================
# STEP 3: Normalization
# - TMM via edgeR  → diversity/PCoA
# - DESeq2 VST     → PCA/heatmaps
# - Raw counts     → use later for DESeq2/edgeR testing
# ============================================================

library(edgeR)
library(DESeq2)
library(ggplot2)
library(ggrepel)

dir.create("02_QC_normalized", showWarnings = FALSE)

metadata <- read.table("groups_filtered.tsv",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)

colnames(metadata) <- c("SampleID", "Group")
metadata$Group <- factor(metadata$Group, levels = c("healthy", "infected"))

domains <- c("bacteria", "virus", "eukaryote")

for (domain in domains) {
  
  cat("\n=============================\n")
  cat("Normalizing:", domain, "\n")
  cat("=============================\n")
  
  infile <- paste0(domain, "_filtered_count_matrix.tsv")
  
  if (!file.exists(infile)) {
    cat("File not found:", infile, "- skipping\n")
    next
  }
  
  counts <- read.table(infile,
                       header = TRUE,
                       row.names = 1,
                       sep = "\t",
                       check.names = FALSE)
  
  counts[is.na(counts)] <- 0
  counts <- round(as.matrix(counts))
  
  common_samples <- intersect(colnames(counts), metadata$SampleID)
  
  if (length(common_samples) < 2) {
    cat("Too few matched samples for", domain, "- skipping\n")
    next
  }
  
  counts <- counts[, common_samples, drop = FALSE]
  meta <- metadata[match(common_samples, metadata$SampleID), , drop = FALSE]
  
  cat("Before zero-sample removal:",
      ncol(counts), "samples |", nrow(counts), "genes\n")
  
  # ------------------------------------------------------------
  # Remove samples with zero library size
  # Important for virus because some samples may have no viral reads
  # ------------------------------------------------------------
  
  lib_sizes <- colSums(counts)
  zero_samples <- names(lib_sizes[lib_sizes == 0])
  
  if (length(zero_samples) > 0) {
    cat("Removing zero-library samples:\n")
    print(zero_samples)
    
    keep_samples <- lib_sizes > 0
    counts <- counts[, keep_samples, drop = FALSE]
    meta <- meta[keep_samples, , drop = FALSE]
  }
  
  cat("After zero-sample removal:",
      ncol(counts), "samples |", nrow(counts), "genes\n")
  
  if (ncol(counts) < 2) {
    cat("Too few non-zero samples for", domain, "- skipping\n")
    next
  }
  
  if (length(unique(meta$Group)) < 2) {
    cat("Only one group remains for", domain, "- skipping\n")
    next
  }
  
  # Remove genes with zero counts after sample filtering
  keep_genes <- rowSums(counts) > 0
  counts <- counts[keep_genes, , drop = FALSE]
  
  cat("After removing zero genes:",
      ncol(counts), "samples |", nrow(counts), "genes\n")
  
  if (nrow(counts) < 2) {
    cat("Too few genes for", domain, "- skipping\n")
    next
  }
  
  # ------------------------------------------------------------
  # A) TMM normalization using edgeR
  # ------------------------------------------------------------
  
  dge <- DGEList(counts = counts, group = meta$Group)
  dge <- calcNormFactors(dge, method = "TMM")
  
  cat("TMM normalization factors:\n")
  print(dge$samples[, c("group", "lib.size", "norm.factors")])
  
  tmm_cpm <- cpm(dge,
                 normalized.lib.sizes = TRUE,
                 log = FALSE)
  
  log_tmm_cpm <- cpm(dge,
                     normalized.lib.sizes = TRUE,
                     log = TRUE,
                     prior.count = 1)
  
  write.table(tmm_cpm,
              paste0("02_QC_normalized/", domain, "_TMM_normalized.tsv"),
              sep = "\t",
              quote = FALSE,
              col.names = NA)
  
  write.table(log_tmm_cpm,
              paste0("02_QC_normalized/", domain, "_TMM_log2CPM.tsv"),
              sep = "\t",
              quote = FALSE,
              col.names = NA)
  
  # ------------------------------------------------------------
  # B) DESeq2 VST normalization
  # ------------------------------------------------------------
  
  vst_matrix <- NULL
  
  if (nrow(counts) >= 10) {
    
    dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData = meta,
      design = ~ Group
    )
    
    # Remove genes with zero total count
    dds <- dds[rowSums(counts(dds)) > 0, ]
    
    dds <- estimateSizeFactors(dds)
    
    cat("\nDESeq2 size factors:\n")
    print(sizeFactors(dds))
    
    dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData = meta,
      design = ~ Group
    )
    
    dds <- dds[rowSums(counts(dds)) > 0, ]
    
    dds <- estimateSizeFactors(dds)
    
    cat("\nDESeq2 size factors:\n")
    print(sizeFactors(dds))
    
    if (nrow(dds) < 1000) {
      cat("Small gene set detected. Using varianceStabilizingTransformation()\n")
      vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    } else {
      vsd <- vst(dds, blind = TRUE)
    }
    
    vst_matrix <- assay(vsd)
    vst_matrix <- assay(vsd)
    
    write.table(vst_matrix,
                paste0("02_QC_normalized/", domain, "_DESeq2_VST.tsv"),
                sep = "\t",
                quote = FALSE,
                col.names = NA)
    
  } else {
    cat("Too few genes for VST. Skipping VST for", domain, "\n")
  }
  
  # ------------------------------------------------------------
  # C) QC plots
  # ------------------------------------------------------------
  
  raw_melt <- stack(as.data.frame(log2(counts + 1)))
  raw_melt$Stage <- "Raw counts"
  colnames(raw_melt) <- c("value", "Sample", "Stage")
  
  tmm_melt <- stack(as.data.frame(log_tmm_cpm))
  tmm_melt$Stage <- "TMM normalized"
  colnames(tmm_melt) <- c("value", "Sample", "Stage")
  
  combined <- rbind(raw_melt, tmm_melt)
  
  p_box <- ggplot(combined,
                  aes(x = Sample, y = value, fill = Stage)) +
    geom_boxplot(outlier.size = 0.3,
                 position = position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     size = 6)) +
    ylab("log2 value") +
    xlab("") +
    ggtitle(paste(domain, "- Raw vs TMM normalized distribution"))
  
  ggsave(paste0("02_QC_normalized/", domain, "_raw_vs_TMM_boxplot.pdf"),
         p_box,
         width = 14,
         height = 6)
  
  p_dens <- ggplot(combined,
                   aes(x = value, color = Sample, linetype = Stage)) +
    geom_density() +
    theme_bw() +
    theme(legend.position = "none") +
    xlab("log2 value") +
    ggtitle(paste(domain, "- Count density raw vs TMM"))
  
  ggsave(paste0("02_QC_normalized/", domain, "_density_raw_vs_TMM.pdf"),
         p_dens,
         width = 8,
         height = 5)
  
  # ------------------------------------------------------------
  # D) PCA using VST matrix
  # ------------------------------------------------------------
  
  if (!is.null(vst_matrix)) {
    
    pca_res <- prcomp(t(vst_matrix), scale. = FALSE)
    
    pca_df <- as.data.frame(pca_res$x[, 1:2])
    pca_df$SampleID <- rownames(pca_df)
    pca_df <- merge(pca_df, meta, by = "SampleID")
    
    var_exp <- round(100 * summary(pca_res)$importance[2, 1:2], 1)
    
    p_pca <- ggplot(pca_df,
                    aes(x = PC1,
                        y = PC2,
                        color = Group,
                        label = SampleID)) +
      geom_point(size = 3) +
      geom_text_repel(size = 2.5, max.overlaps = 20) +
      theme_bw() +
      xlab(paste0("PC1 (", var_exp[1], "%)")) +
      ylab(paste0("PC2 (", var_exp[2], "%)")) +
      ggtitle(paste(domain, "- PCA on VST normalized counts"))
    
    ggsave(paste0("02_QC_normalized/", domain, "_VST_PCA.pdf"),
           p_pca,
           width = 8,
           height = 6)
  }
  
  # ------------------------------------------------------------
  # E) Save normalization summary
  # ------------------------------------------------------------
  
  summary_df <- data.frame(
    Domain = domain,
    Samples_used = ncol(counts),
    Genes_used = nrow(counts),
    Healthy_samples = sum(meta$Group == "healthy"),
    Infected_samples = sum(meta$Group == "infected"),
    Min_library_size = min(colSums(counts)),
    Median_library_size = median(colSums(counts)),
    Max_library_size = max(colSums(counts))
  )
  
  write.table(summary_df,
              paste0("02_QC_normalized/", domain, "_normalization_summary.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  cat("Saved normalization outputs for", domain, "\n")
}

cat("\nDONE.\n")
cat("Check 02_QC_normalized/ folder.\n")
cat("Use *_TMM_normalized.tsv for diversity/PCoA.\n")
cat("Use *_DESeq2_VST.tsv for PCA/heatmaps.\n")
cat("Use raw filtered count matrices for differential abundance testing.\n")