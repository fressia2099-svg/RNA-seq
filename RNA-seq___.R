install.packages(c(
  "tidyverse",   # 数据处理
  "data.table",  # 大数据处理
  "ggplot2",     # 可视化
  "readr",       # 数据导入
  "dplyr"        # 数据操作
))

install.packages("tidyverse")
install.packages("data.table")
install.packages("readxl")
install.packages("ggrepel")

# Bioconductor
install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("ComplexHeatmap")
BiocManager::install("tximport")
BiocManager::install("edgeR")


setwd("F:/R_Studio Files")

getwd()

# 重新读取，跳过第一行，用第二行作为列名
df <- read.csv("RNAseq GENES.csv", skip = 1)

# 验证
colnames(df)
head(df)

library(dplyr)

# Filter DEGs
deg <- df %>%
  filter(FDR < 0.05, abs(logFC) > 1)

# Summary
cat("Total DEGs:", nrow(deg), "\n")
cat("Upregulated:", sum(deg$logFC > 0), "\n")
cat("Downregulated:", sum(deg$logFC < 0), "\n")

# Export
write.csv(deg, "DEGs_filtered.csv", row.names = FALSE)


library(ggplot2)
library(dplyr)
library(ggrepel)

# Top 20 by FDR
top_fdr <- df %>%
  filter(significance != "Not Significant") %>%
  arrange(FDR) %>%
  head(20)

# Top 20 by logFC (10 most up + 10 most down)
top_fc <- df %>%
  filter(significance != "Not Significant") %>%
  arrange(desc(abs(logFC))) %>%
  head(20)

# Combine and remove duplicates
top_genes <- bind_rows(top_fdr, top_fc) %>%
  distinct(EnsemblID, .keep_all = TRUE)

# Volcano plot with labels
ggplot(df, aes(x = logFC, y = -log10(FDR), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated"     = "#E41A1C",
    "Downregulated"   = "#377EB8",
    "Not Significant" = "grey60"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_label_repel(
    data = top_genes,
    aes(label = gene_name),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.4,
    show.legend = FALSE
  ) +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(FDR)",
    color = "Regulation"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave("volcano_plot_labeled.png", width = 10, height = 7, dpi = 300)








BiocManager::install("ComplexHeatmap")


library(dplyr)
library(ComplexHeatmap)
library(circlize)

# Select top 50 DEGs by FDR
top50 <- deg %>%
  arrange(FDR) %>%
  head(50)

# Build expression matrix
mat <- as.matrix(top50[, c("GroupA_AvgCPM", "GroupB_AvgCPM")])
rownames(mat) <- top50$gene_name
colnames(mat) <- c("Control", "Ampicillin")

# Z-score scaling by row
mat_scaled <- t(scale(t(mat)))

# Color scale
col_fun <- colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C"))

# Save heatmap
png("heatmap_top50.png", width = 6, height = 10, units = "in", res = 300)

Heatmap(
  mat_scaled,
  name = "Z-score",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_title = "Top 50 DEGs (by FDR)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "Z-score")
)

dev.off()



#  GO/KEGG富集分析


BiocManager::install("clusterProfiler")
BiocManager::install("org.Sc.sgd.db")  # 酵母注释数据库
BiocManager::install("org.CalbSC5314.eg.db")

head(deg$EnsemblID)

# Export gene name list for GO analysis
write(deg$gene_name, "gene_list_for_GO.txt", sep = "\n")





#GO KEGG 需要在网站上分析,网站分析后获得tsv表格,下载后继续:


go_results <- read.delim("go_term_finder_results.tsv", header = TRUE)
head(go_results)
colnames(go_results)


nrow(go_results)



#Biological Process画图:

library(ggplot2)

# Select top 20 by FDR
top20_go <- go_results %>%
  arrange(FDR) %>%
  head(20)

# Dotplot
ggplot(top20_go, aes(
  x = Fold.Enrichment,
  y = reorder(GO.Term, Fold.Enrichment),
  size = Query.Count,
  color = FDR
)) +
  geom_point() +
  scale_color_gradient(low = "#E41A1C", high = "#377EB8") +
  scale_size(range = c(3, 10)) +
  labs(
    title = "GO Enrichment - Biological Process (Top 20)",
    x = "Fold Enrichment",
    y = "GO Term",
    size = "Gene Count",
    color = "FDR"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8)
  )

ggsave("GO_dotplot.png", width = 10, height = 8, dpi = 800)







#Molecular Function 和 Cellular Component 画图

getwd()

mf_results <- read.delim("go_term_finder_Molecular Function.tsv", header = TRUE)
cc_results <- read.delim("go_term_finder_Cellular Component.tsv", header = TRUE)

nrow(mf_results)
nrow(cc_results)



library(ggplot2)
library(dplyr)

# Add ontology label to each
bp_top <- go_results %>% arrange(FDR) %>% head(10) %>% mutate(Ontology = "Biological Process")
mf_top <- mf_results %>% arrange(FDR) %>% head(10) %>% mutate(Ontology = "Molecular Function")
cc_top <- cc_results %>% arrange(FDR) %>% head(10) %>% mutate(Ontology = "Cellular Component")

# Combine
all_go <- bind_rows(bp_top, mf_top, cc_top)

# Dotplot
ggplot(all_go, aes(
  x = Fold.Enrichment,
  y = reorder(GO.Term, Fold.Enrichment),
  size = Query.Count,
  color = FDR
)) +
  geom_point() +
  scale_color_gradient(low = "#E41A1C", high = "#377EB8") +
  scale_size(range = c(3, 8)) +
  facet_wrap(~ Ontology, scales = "free_y", ncol = 1) +
  labs(
    title = "GO Enrichment Analysis",
    x = "Fold Enrichment",
    y = "GO Term",
    size = "Gene Count",
    color = "FDR"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 10)
  )

ggsave("GO_dotplot_all.png", width = 12, height = 14, dpi = 800)

















#KEGG
library(clusterProfiler)

summary(kegg_result)

# Check what your IDs look like vs KEGG format
# Your format:   C1_00070W_A
# KEGG format:   CAALFM_C100070WA

# Try converting by reformatting
deg$KEGG_ID <- gsub("_", "", deg$EnsemblID)
deg$KEGG_ID <- paste0("CAALFM_", deg$KEGG_ID)

head(deg$KEGG_ID)



kegg_result <- enrichKEGG(
  gene = deg$KEGG_ID,
  organism = "cal",
  pvalueCutoff = 0.05
)

summary(kegg_result)
nrow(kegg_result)


# Dotplot
dotplot(kegg_result, showCategory = 15, title = "KEGG Pathway Enrichment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("KEGG_dotplot.png", width = 10, height = 8, dpi = 800)
