# Load required libraries
library(readr)
library(dplyr)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(enrichplot)

# Step 1: Import quantification results from Salmon
quant <- read_tsv("quant_ENCLB172MPH/quant.sf")
quant$GeneID <- sub(".*\\|(ENSG[0-9]+)\\..*", "\\1", quant$Name)

# Summarize transcript-level expression to gene-level counts
my_gene_counts <- quant %>%
  group_by(GeneID) %>%
  summarise(my_counts = sum(NumReads)) %>%
  ungroup()

# Step 2: Import NCBI-generated count tables for treated and control groups
treated <- read_tsv("GSE219610_raw_counts_GRCh38.p13_NCBI.tsv.gz")
control <- read_tsv("GSE86657_raw_counts_GRCh38.p13_NCBI.tsv.gz")

# Extract one column from NCBI-treated data for comparison
col_to_use <- "GSM6782908"
ncbi_counts <- treated %>%
  select(GeneID, all_of(col_to_use)) %>%
  setNames(c("GeneID", "ncbi_counts")) %>%
  mutate(GeneID = as.character(GeneID))

# Convert Ensembl Gene IDs to Entrez IDs
id_map <- bitr(my_gene_counts$GeneID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
my_gene_entrez <- inner_join(my_gene_counts, id_map, by = c("GeneID" = "ENSEMBL"))

# Join my quantification and NCBI counts by Entrez ID
compare_df <- inner_join(
  dplyr::select(my_gene_entrez, ENTREZID, my_counts),
  dplyr::rename(ncbi_counts, ENTREZID = GeneID),
  by = "ENTREZID"
)

# Remove zero values before log transformation
compare_log <- compare_df %>% filter(my_counts > 0 & ncbi_counts > 0)

# Scatter plot comparison
p_counts <- ggplot(compare_log, aes(x = my_counts, y = ncbi_counts)) +
  geom_point(alpha = 0.3, size = 0.7, color = "steelblue") +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Gene-Level Count Comparison: My Salmon vs NCBI",
       x = "My Quantification (log10)",
       y = "NCBI Staff Quantification (log10)") +
  theme_minimal()
ggsave("salmon_vs_ncbi_counts.pdf", plot = p_counts, width = 7, height = 6)

# Spearman correlation
cor(compare_log$my_counts, compare_log$ncbi_counts, method = "spearman")

# Step 3: Differential expression analysis
merged_counts <- inner_join(treated, control, by = "GeneID")
count_matrix <- merged_counts[, -1]
rownames(count_matrix) <- merged_counts$GeneID
colData <- data.frame(row.names = colnames(count_matrix), condition = factor(c(rep("treated", 5), rep("control", 2))))
dds <- DESeqDataSetFromMatrix(countData = round(count_matrix), colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), "DESeq2_results_treated_vs_control.csv")

# Volcano plot for full dataset
res <- read.csv("DESeq2_results_treated_vs_control.csv", row.names = 1)
res$threshold <- "Not Sig"
res$threshold[res$padj < 0.05 & res$log2FoldChange > 1] <- "Up"
res$threshold[res$padj < 0.05 & res$log2FoldChange < -1] <- "Down"
res$log10padj <- -log10(pmax(res$padj, 1e-300))
res$ENTREZID <- rownames(res)
gene_info <- bitr(res$ENTREZID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
res_annotated <- left_join(res, gene_info, by = "ENTREZID")
top_up <- res_annotated %>% filter(threshold == "Up") %>% arrange(padj) %>% slice_head(n = 10)
top_down <- res_annotated %>% filter(threshold == "Down") %>% arrange(padj) %>% slice_head(n = 10)
top_genes <- bind_rows(top_up, top_down)
res_annotated$label <- ifelse(res_annotated$ENTREZID %in% top_genes$ENTREZID, res_annotated$SYMBOL, NA)
p <- ggplot(res_annotated, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(aes(color = threshold), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Down" = "blue", "Up" = "red", "Not Sig" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 100) +
  labs(title = "Volcano Plot (Top 10 Up & Down Genes Labeled)",
       x = "log2(Fold Change)",
       y = "-log10(adjusted p-value)",
       color = "Significance") +
  theme_minimal()
ggsave("volcano_top10_labeled.pdf", plot = p, width = 8, height = 6, dpi = 300)

# DESeq2 with 3 treated vs 2 control
treated_subset <- c("GSM6782904", "GSM6782905", "GSM6782906")
control_subset <- c("GSM2308412", "GSM2308413")
count_matrix_sub <- count_matrix[, c(treated_subset, control_subset)]
colData_sub <- data.frame(row.names = colnames(count_matrix_sub), condition = factor(c(rep("treated", 3), rep("control", 2))))
dds_sub <- DESeqDataSetFromMatrix(countData = round(count_matrix_sub), colData = colData_sub, design = ~ condition)
dds_sub <- DESeq(dds_sub)
res_sub <- results(dds_sub)
write.csv(as.data.frame(res_sub), "DESeq2_results_3treated_vs_2control.csv")

# Compare DEG counts
sig_5 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
sig_3 <- subset(res_sub, padj < 0.05 & abs(log2FoldChange) > 1)
cat("Number of significant DEGs (5 replicates):", nrow(sig_5), "\n")
cat("Number of significant DEGs (3 replicates):", nrow(sig_3), "\n")

# log2FC correlation plot
res$gene <- rownames(res)
res_sub$gene <- rownames(res_sub)
res_merge <- inner_join(
  res %>% select(gene, log2FoldChange) %>% rename(logFC_5 = log2FoldChange),
  res_sub %>% select(gene, log2FoldChange) %>% rename(logFC_3 = log2FoldChange),
  by = "gene"
)
p <- ggplot(res_merge, aes(x = logFC_3, y = logFC_5)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Comparison of log2 Fold Changes (3 vs 5 replicates)",
       x = "log2FC (3 treated vs 2 control)",
       y = "log2FC (5 treated vs 5 control)") +
  theme_minimal()
ggsave("log2FC_comparison_3vs5.pdf", plot = p, width = 7, height = 6)

# Venn diagram of DEGs
venn.plot <- draw.pairwise.venn(
  area1 = nrow(sig_5),
  area2 = nrow(sig_3),
  cross.area = length(intersect(rownames(sig_5), rownames(sig_3))),
  category = c("5 replicates", "3 replicates"),
  fill = c("skyblue", "pink"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2
)
ggsave("Venn.pdf", plot = venn.plot, width = 10, height = 7)

# GO enrichment (BP, CC, MF)
for (ont in c("BP", "MF", "CC")) {
  ego <- enrichGO(
    gene = gene_up,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  if (!is.null(ego) && nrow(ego@result) > 0) {
    pdf(paste0("GO_", ont, "_dotplot.pdf"), width = 12, height = 8)
    print(dotplot(ego, showCategory = 15, title = paste("GO", ont, "Enrichment (Upregulated)")))
    dev.off()
  }
}

# KEGG enrichment
ekegg_up <- enrichKEGG(
  gene = gene_up,
  organism = 'hsa',
  pvalueCutoff = 0.05
)
write.csv(as.data.frame(ekegg_up), "KEGG_upregulated.csv")
pdf("KEGG_barplot.pdf", width = 9, height = 6)
barplot(ekegg_up, showCategory = 15, title = "KEGG Pathway Enrichment (Upregulated)", font.size = 12)
dev.off()
