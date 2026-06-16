# ============================================================
# STEP 4: Alpha diversity, Bray-Curtis PCoA, PERMANOVA
# Input:  *_TMM_normalized.tsv  (from step3)
#         groups_filtered.tsv
# Output: 04_diversity_PCoA/
# ============================================================

library(vegan)
library(ggplot2)

dir.create("04_diversity_PCoA", showWarnings = FALSE)

metadata <- read.table("groups_filtered.tsv",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)
colnames(metadata) <- c("SampleID", "Group")
metadata$Group <- factor(metadata$Group, levels = c("healthy", "infected"))

domains <- c("bacteria", "virus", "eukaryote")

for (domain in domains) {
  
  cat("\n=============================\n")
  cat("Diversity analysis:", domain, "\n")
  cat("=============================\n")
  
  # ---- Load TMM normalized matrix ----
  tmm <- read.table(paste0("QC_normalized/", domain, "_TMM_normalized.tsv"),
                    header = TRUE,
                    row.names = 1,
                    sep = "\t",
                    check.names = FALSE)
  
  tmm <- as.matrix(tmm)
  
  # Match metadata to samples
  common_samples <- intersect(colnames(tmm), metadata$SampleID)
  tmm  <- tmm[, common_samples, drop = FALSE]
  meta <- metadata[match(common_samples, metadata$SampleID), ]
  
  cat("Samples:", ncol(tmm), "| Genes:", nrow(tmm), "\n")
  
  # Transpose: vegan expects samples as rows
  sample_table <- t(tmm)
  
  # ==========================
  # A) Alpha diversity
  # ==========================
  
  alpha <- data.frame(
    SampleID = rownames(sample_table),
    Observed  = specnumber(sample_table),       # number of detected genes
    Shannon   = diversity(sample_table, index = "shannon"),
    Simpson   = diversity(sample_table, index = "simpson"),
    Evenness  = diversity(sample_table, index = "shannon") / log(specnumber(sample_table))
  )
  
  alpha <- merge(alpha, meta, by = "SampleID")
  
  write.table(alpha,
              paste0("04_diversity_PCoA/", domain, "_alpha_diversity.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Shannon boxplot
  p_shannon <- ggplot(alpha, aes(x = Group, y = Shannon, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2.5) +
    theme_bw() +
    scale_fill_manual(values = c("healthy" = "#2196F3", "infected" = "#F44336")) +
    ggtitle(paste(domain, "- Shannon diversity")) +
    ylab("Shannon index") +
    xlab("")
  
  ggsave(paste0("04_diversity_PCoA/", domain, "_shannon_diversity.pdf"),
         p_shannon, width = 5, height = 4)
  
  # Simpson boxplot
  p_simpson <- ggplot(alpha, aes(x = Group, y = Simpson, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2.5) +
    theme_bw() +
    scale_fill_manual(values = c("healthy" = "#2196F3", "infected" = "#F44336")) +
    ggtitle(paste(domain, "- Simpson diversity")) +
    ylab("Simpson index") +
    xlab("")
  
  ggsave(paste0("04_diversity_PCoA/", domain, "_simpson_diversity.pdf"),
         p_simpson, width = 5, height = 4)
  
  # Wilcoxon test on Shannon
  healthy_shannon  <- alpha$Shannon[alpha$Group == "healthy"]
  infected_shannon <- alpha$Shannon[alpha$Group == "infected"]
  wtest <- wilcox.test(healthy_shannon, infected_shannon)
  
  cat("Wilcoxon test Shannon healthy vs infected:\n")
  cat("  W =", wtest$statistic, "| p-value =", round(wtest$p.value, 4), "\n")
  
  sink(paste0("04_diversity_PCoA/", domain, "_alpha_wilcoxon.txt"))
  cat("Wilcoxon test: Shannon diversity healthy vs infected\n\n")
  print(wtest)
  sink()
  
  # ==========================
  # B) Beta diversity — Bray-Curtis distance
  # ==========================
  
  bray <- vegdist(sample_table, method = "bray")
  
  # Save distance matrix
  write.table(as.matrix(bray),
              paste0("04_diversity_PCoA/", domain, "_bray_distance_matrix.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
  
  # ==========================
  # C) PCoA
  # ==========================
  
  pcoa <- cmdscale(bray, eig = TRUE, k = 2)
  
  # Variance explained by each axis
  eig_vals  <- pcoa$eig
  var_exp   <- round(100 * eig_vals / sum(eig_vals[eig_vals > 0]), 1)
  
  pcoa_df <- data.frame(
    SampleID = rownames(sample_table),
    PC1 = pcoa$points[, 1],
    PC2 = pcoa$points[, 2]
  )
  
  pcoa_df <- merge(pcoa_df, meta, by = "SampleID")
  
  write.table(pcoa_df,
              paste0("04_diversity_PCoA/", domain, "_bray_pcoa.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group, label = SampleID)) +
    geom_point(size = 3.5) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20) +
    theme_bw() +
    scale_color_manual(values = c("healthy" = "#2196F3", "infected" = "#F44336")) +
    xlab(paste0("PC1 (", var_exp[1], "%)")) +
    ylab(paste0("PC2 (", var_exp[2], "%)")) +
    ggtitle(paste(domain, "- Bray-Curtis PCoA"))
  
  ggsave(paste0("04_diversity_PCoA/", domain, "_bray_pcoa.pdf"),
         p_pcoa, width = 7, height = 6)
  
  # ==========================
  # D) PERMANOVA
  # ==========================
  
  set.seed(123)   # for reproducibility
  permanova <- adonis2(bray ~ Group, data = meta, permutations = 999)
  
  cat("PERMANOVA result:\n")
  print(permanova)
  
  sink(paste0("04_diversity_PCoA/", domain, "_PERMANOVA.txt"))
  cat("PERMANOVA: Bray-Curtis distance ~ Group\n")
  cat("Permutations: 999\n\n")
  print(permanova)
  sink()
  
  # ==========================
  # E) Sample clustering dendrogram
  # ==========================
  
  dist_clust <- dist(log2(sample_table + 1))
  hc <- hclust(dist_clust, method = "ward.D2")
  
  pdf(paste0("04_diversity_PCoA/", domain, "_sample_clustering.pdf"),
      width = 10, height = 6)
  plot(hc,
       main = paste(domain, "- Sample clustering (ward.D2, log2 TMM)"),
       xlab = "", sub = "")
  # Add colored bar for group
  group_colors <- ifelse(meta$Group[match(hc$labels, meta$SampleID)] == "healthy",
                         "#2196F3", "#F44336")
  colored_bars <- matrix(group_colors, nrow = 1)
  rownames(colored_bars) <- "Group"
  dev.off()
  
  cat("Done:", domain, "\n")
}

cat("\n\nDONE step4. Check 04_diversity_PCoA/ folder.\n")
cat("Key outputs:\n")
cat("  *_alpha_diversity.tsv     - Shannon, Simpson, Evenness per sample\n")
cat("  *_alpha_wilcoxon.txt      - statistical test healthy vs infected\n")
cat("  *_bray_pcoa.pdf           - PCoA plot\n")
cat("  *_PERMANOVA.txt           - is community difference significant?\n")
cat("  *_sample_clustering.pdf   - dendrogram\n")