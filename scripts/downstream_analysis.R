# downstream_analysis.R
# LungMAP AT2 Cell Cluster Analysis and Downstream RNA-seq Analyses
# Author: Patrick Ozark
# Date: 2026-03-13

# =============================================================================
# SETUP
# =============================================================================

# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(vsn)
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse gene annotations
library(enrichplot)
library(DOSE)
library(pathview)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(RColorBrewer)

# Set working directory
setwd("~/Projects/lungMAPmouse/")

# Load the rlog-normalized expression matrix (already computed in main analysis)
# If not already loaded, run: load("./lungMAPmouse.RData")
# The object rld_at2_matrix should be available in your environment

# =============================================================================
# CLUSTER SUBSETTING: Extract Cluster 6 Genes from rld_at2_matrix
# =============================================================================

# Read in cluster assignments from hierarchical clustering
clusters <- read.table("./data/AT2_gene_hclust_clusters_k12.csv", 
                       sep = ",", 
                       header = TRUE, 
                       stringsAsFactors = FALSE)

# Extract cluster 6 specifically
cluster_6 <- subset(clusters, cluster == 6)

# Get gene names as character vector
cluster6_genes <- as.character(cluster_6$gene)

# Verify genes exist in rld_at2_matrix
# (some genes may have been filtered during DESeq2 processing)
genes_in_matrix <- intersect(cluster6_genes, rownames(rld_at2_matrix))
genes_missing <- setdiff(cluster6_genes, rownames(rld_at2_matrix))

if (length(genes_missing) > 0) {
  message(paste0("Warning: ", length(genes_missing), " genes from cluster 6 not found in rld_at2_matrix"))
}

# Subset rld_at2_matrix to cluster 6 genes only
cluster6_matrix <- rld_at2_matrix[genes_in_matrix, , drop = FALSE]

# Save cluster6_matrix for future use
saveRDS(cluster6_matrix, "./data/AT2_cluster6_rlog_matrix.rds")

# =============================================================================
# DOWNSTREAM ANALYSIS 1: Heatmap Visualization of Cluster 6
# =============================================================================

# Create heatmap for cluster 6 genes
# Using z-score scaling for better visualization
cluster6_scaled <- t(scale(t(cluster6_matrix)))

# Define color palette
heat_colors <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

# ComplexHeatmap version (more customizable)
pdf("./plots/AT2_cluster6_heatmap.pdf", width = 8, height = 10)
ComplexHeatmap::Heatmap(
  cluster6_scaled,
  name = "Z-score",
  col = heat_colors,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_title = "Cluster 6 Genes Across Samples",
  row_title = paste0("Genes (n=", nrow(cluster6_scaled), ")"),
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D2",
  top_annotation = ComplexHeatmap::HeatmapAnnotation(
    Age = colData(rld_at2)$age,
    col = list(Age = c("P7" = "#4DAF4A", "P14" = "#377EB8", "P21" = "#E41A1C"))
  )
)
dev.off()

# =============================================================================
# DOWNSTREAM ANALYSIS 2: Gene Ontology (GO) Enrichment Analysis
# =============================================================================

# Convert gene symbols to Entrez IDs for cluster 6
cluster6_entrez <- bitr(
  cluster6_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

# Remove NAs
cluster6_entrez <- na.omit(cluster6_entrez)

# GO Biological Process enrichment
go_bp <- enrichGO(
  gene = cluster6_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# GO Molecular Function enrichment
go_mf <- enrichGO(
  gene = cluster6_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# GO Cellular Component enrichment
go_cc <- enrichGO(
  gene = cluster6_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Save GO results
write.csv(as.data.frame(go_bp), "./data/AT2_cluster6_GO_BP_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_mf), "./data/AT2_cluster6_GO_MF_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_cc), "./data/AT2_cluster6_GO_CC_results.csv", row.names = FALSE)

# Plot top GO terms
pdf("./plots/AT2_cluster6_GO_BP_dotplot.pdf", width = 10, height = 8)
dotplot(go_bp, showCategory = 20) + 
  ggtitle("Cluster 6: Top GO Biological Process Terms")
dev.off()

pdf("./plots/AT2_cluster6_GO_BP_barplot.pdf", width = 10, height = 8)
barplot(go_bp, showCategory = 20) + 
  ggtitle("Cluster 6: Top GO Biological Process Terms")
dev.off()

# =============================================================================
# DOWNSTREAM ANALYSIS 3: KEGG Pathway Enrichment
# =============================================================================

# KEGG pathway enrichment
kegg_enrich <- enrichKEGG(
  gene = cluster6_entrez$ENTREZID,
  organism = "mmu",  # Mouse
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Convert KEGG IDs to gene symbols for readability
kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

# Save KEGG results
write.csv(as.data.frame(kegg_enrich), "./data/AT2_cluster6_KEGG_results.csv", row.names = FALSE)

# Plot KEGG pathways
if (nrow(as.data.frame(kegg_enrich)) > 0) {
  pdf("./plots/AT2_cluster6_KEGG_dotplot.pdf", width = 10, height = 8)
  dotplot(kegg_enrich, showCategory = 15) + 
    ggtitle("Cluster 6: KEGG Pathway Enrichment")
  dev.off()
}

# Visualize specific pathway (example: if enriched)
# Replace "mmu04060" with actual enriched pathway ID
# pathview(gene.data = cluster6_entrez$ENTREZID,
#          pathway.id = "mmu04060",
#          species = "mmu",
#          out.suffix = "cluster6")

# =============================================================================
# DOWNSTREAM ANALYSIS 4: Reactome Pathway Enrichment
# =============================================================================

# Reactome pathway enrichment (alternative pathway database)
# library(ReactomePA)
# reactome_enrich <- enrichPathway(
#   gene = cluster6_entrez$ENTREZID,
#   organism = "mouse",
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.2,
#   readable = TRUE
# )

# =============================================================================
# DOWNSTREAM ANALYSIS 5: Gene Set Enrichment Analysis (GSEA)
# =============================================================================

# If you have a ranked list of genes (e.g., from differential expression)
# Create a ranked gene list from correlation results
# Using the Spearman correlation with Il34 from earlier analysis

# Read the correlation results (if saved previously)
il34_cor_results <- read.csv("./data/AT2_Il34_spearman_by_gene.csv", stringsAsFactors = FALSE)

# Create ranked list for GSEA
ranked_genes <- il34_cor_results %>%
  dplyr::select(gene, rho) %>%
  na.omit() %>%
  arrange(desc(rho))

# Convert to named vector for GSEA
ranked_vector <- setNames(ranked_genes$rho, ranked_genes$gene)

# Run GSEA with GO terms
gsea_result <- GSEA(
  geneList = ranked_vector,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "SYMBOL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Plot GSEA results
pdf("./plots/AT2_Il34_GSEA_dotplot.pdf", width = 10, height = 8)
dotplot(gsea_result, showCategory = 15) + 
  ggtitle("GSEA: Gene Sets Correlated with Il34 Expression")
dev.off()

# Ridge plot for enriched pathways
pdf("./plots/AT2_Il34_GSEA_ridgeplot.pdf", width = 10, height = 8)
ridgeplot(gsea_result, showCategory = 15) + 
  ggtitle("GSEA Ridge Plot: Il34-Correlated Pathways")
dev.off()

# =============================================================================
# DOWNSTREAM ANALYSIS 6: Gene Correlation Network Analysis
# =============================================================================

# Build gene-gene correlation network for cluster 6 genes
cor_matrix <- cor(t(cluster6_matrix), method = "spearman")

# Set diagonal to 0
diag(cor_matrix) <- 0

# Filter for strong correlations (|r| > 0.7)
threshold <- 0.7
cor_matrix_filtered <- cor_matrix
cor_matrix_filtered[abs(cor_matrix_filtered) < threshold] <- 0

# Create adjacency matrix for network analysis
adj_matrix <- abs(cor_matrix_filtered)

# Plot correlation heatmap
pdf("./plots/AT2_cluster6_gene_correlation_heatmap.pdf", width = 12, height = 12)
pheatmap(
  cor_matrix[1:min(50, nrow(cor_matrix)), 1:min(50, ncol(cor_matrix))],
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Cluster 6 Gene-Gene Correlation Matrix (Top 50 Genes)"
)
dev.off()

# =============================================================================
# DOWNSTREAM ANALYSIS 7: Principal Component Analysis (PCA)
# =============================================================================

# PCA on cluster 6 genes
pca_result <- prcomp(t(cluster6_matrix), center = TRUE, scale. = TRUE)

# Get variance explained
var_explained <- summary(pca_result)$importance[2, ] * 100

# Create PCA data frame
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = rownames(pca_result$x),
  Age = colData(rld_at2)$age[match(rownames(pca_result$x), colnames(rld_at2_matrix))]
)

# Plot PCA
pdf("./plots/AT2_cluster6_PCA.pdf", width = 8, height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Age)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 3) +
  labs(
    title = "PCA: Cluster 6 Gene Expression",
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")
  ) +
  scale_color_manual(values = c("P7" = "#4DAF4A", "P14" = "#377EB8", "P21" = "#E41A1C")) +
  theme_bw()
dev.off()

# Loadings plot (top contributing genes)
loadings_df <- data.frame(
  Gene = rownames(pca_result$rotation),
  PC1_loading = abs(pca_result$rotation[, 1]),
  PC2_loading = abs(pca_result$rotation[, 2])
) %>%
  arrange(desc(PC1_loading))

write.csv(loadings_df[1:50, ], "./data/AT2_cluster6_PCA_top_loadings.csv", row.names = FALSE)

# =============================================================================
# DOWNSTREAM ANALYSIS 8: Differential Expression Within Cluster
# =============================================================================

# If comparing age groups within cluster 6 genes
# This requires the original DESeq2 dataset

# Subset DESeq2 results for cluster 6 genes only
at2_results <- results(at2_dds, tidy = TRUE)
cluster6_deg <- at2_results %>%
  filter(row %in% cluster6_genes) %>%
  arrange(padj)

write.csv(cluster6_deg, "./data/AT2_cluster6_DESeq2_results.csv", row.names = FALSE)

# Volcano plot for cluster 6 genes
cluster6_deg <- cluster6_deg %>%
  mutate(
    significance = case_when(
      padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
      padj < 0.05 ~ "Adj. p < 0.05",
      abs(log2FoldChange) > 1 ~ "|log2FC| > 1",
      TRUE ~ "Not significant"
    )
  )

pdf("./plots/AT2_cluster6_volcano.pdf", width = 10, height = 8)
ggplot(cluster6_deg, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6) +
  scale_color_manual(
    values = c("Not significant" = "grey",
               "Adj. p < 0.05" = "blue",
               "|log2FC| > 1" = "orange",
               "Significant" = "red")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  labs(
    title = "Volcano Plot: Cluster 6 Genes",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_bw()
dev.off()

# =============================================================================
# DOWNSTREAM ANALYSIS 9: Time-course Analysis
# =============================================================================

# Identify genes with monotonic expression trends across timepoints
# Useful for developmental studies (P7 -> P14 -> P21)

# Calculate mean expression per timepoint for cluster 6
sample_info <- data.frame(
  Sample = colnames(cluster6_matrix),
  Age = colData(rld_at2)$age
)

# Group samples by age
p7_samples <- sample_info$Sample[sample_info$Age == "P7"]
p14_samples <- sample_info$Sample[sample_info$Age == "P14"]
p21_samples <- sample_info$Sample[sample_info$Age == "P21"]

# Calculate mean expression per age
cluster6_mean_expr <- data.frame(
  Gene = rownames(cluster6_matrix),
  P7_mean = rowMeans(cluster6_matrix[, p7_samples, drop = FALSE]),
  P14_mean = rowMeans(cluster6_matrix[, p14_samples, drop = FALSE]),
  P21_mean = rowMeans(cluster6_matrix[, p21_samples, drop = FALSE])
)

# Classify genes by trend
cluster6_mean_expr <- cluster6_mean_expr %>%
  mutate(
    trend = case_when(
      P21_mean > P14_mean & P14_mean > P7_mean ~ "Up-regulated",
      P7_mean > P14_mean & P14_mean > P21_mean ~ "Down-regulated",
      TRUE ~ "Non-monotonic"
    )
  )

# Count genes per trend category
table(cluster6_mean_expr$trend)

write.csv(cluster6_mean_expr, "./data/AT2_cluster6_expression_trends.csv", row.names = FALSE)

# Plot heatmap of up-regulated genes
up_genes <- cluster6_mean_expr$Gene[cluster6_mean_expr$trend == "Up-regulated"]
if (length(up_genes) > 1) {
  pdf("./plots/AT2_cluster6_upregulated_heatmap.pdf", width = 8, height = 10)
  pheatmap(
    t(scale(t(cluster6_matrix[up_genes, ]))),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = paste0("Up-regulated Genes (n=", length(up_genes), ")")
  )
  dev.off()
}

# =============================================================================
# DOWNSTREAM ANALYSIS 10: Export for External Tools
# =============================================================================

# Export cluster 6 gene list for external tools (e.g., STRING, Metascape)
write.table(
  cluster6_genes,
  "./data/AT2_cluster6_gene_list.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Export expression matrix for external analysis
write.csv(
  as.data.frame(cluster6_matrix),
  "./data/AT2_cluster6_expression_matrix.csv"
)

# Export as tab-delimited for more tools
write.table(
  as.data.frame(cluster6_matrix),
  "./data/AT2_cluster6_expression_matrix.txt",
  sep = "\t",
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)

# =============================================================================
# SESSION INFO
# =============================================================================

# Save session info for reproducibility
session_info <- sessionInfo()
writeLines(capture.output(session_info), "./data/downstream_analysis_session_info.txt")

message("Downstream analysis complete. Results saved to ./data/ and ./plots/")