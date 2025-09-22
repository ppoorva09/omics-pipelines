
###############################################
# Differential Expression Analysis Pipeline
# Using DESeq2 and edgeR on RNA-seq Count Data
###############################################

# Load required libraries (Alot of Libraries are present in BiocManager library)
library(DESeq2)    
library(ggplot2)    
library(edgeR)      
library(AnnotationDbi)   
library(org.Mm.eg.db)    
library(ggrepel)         
library(pheatmap)        


# Using FastQC and FeatureCount I was able to get the counts data 
# and later performed the QC and removed the NA values 
# then performed the differential analysis 


# Step 1: Load and filter count data
##########################

Counts <- read.delim("counts_Geneids.csv", header = TRUE, row.names = 1, sep = ",")
head(Counts)

# Remove genes with very low counts across all samples
Counts <- Counts[which(rowSums(Counts) > 15), ]

# Keep only first 11 columns (in case extra metadata is present)
Counts <- Counts[, 1:11]


# Step 2: sample information
##########################

condition <- factor(c("WT", "WT", "WT", "WT", "WT", "WT", 
                      "TX", "TX", "TX", "TX", "TX"))
coldata <- data.frame(row.names = colnames(Counts), condition)
coldata


# Step 3: Run DESeq2 analysis
##########################

dds <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = coldata,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)

res <- results(dds)

# Filter significant results (padj < 0.05 & |log2FC| >= 1)
res_filtered <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) >= 1), ]
write.csv(as.data.frame(res_filtered), "DESeq2_results.csv")

summary(res)


# Step 4: Run edgeR analysis (for comparison)
##########################

group <- coldata$condition
y <- DGEList(counts = Counts, group = group)
y <- calcNormFactors(y)
y <- estimateDisp(y, design = model.matrix(~group))
fit <- glmQLFit(y, design = model.matrix(~group))
qlf <- glmQLFTest(fit)
topTags(qlf, n = 20)


# Step 5: Annotate genes
##########################

res$ENSEMBL <- rownames(res)

res$ENTREZID <- mapIds(org.Mm.eg.db,
                       key = res$ENSEMBL,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")

res$SYMBOL <- mapIds(org.Mm.eg.db,
                     key = res$ENSEMBL,
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")


# Step 6: Volcano Plot
##########################

res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 0.6 & res$pvalue < 0.05] <- "UP"
res$diffexpressed[res$log2FoldChange < -0.6 & res$pvalue < 0.05] <- "DOWN"

p_graph <- ggplot(data = as.data.frame(res), 
                  aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "darkred", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dashed") +
  geom_text_repel(aes(label = ifelse(diffexpressed != "NO", SYMBOL, "")), size = 3) +
  labs(title = "Volcano Plot of Differential Expression", 
       x = "Log2 Fold Change", 
       y = "-Log10(p-value)") +
  xlim(-10, 10) +  
  ylim(0, 10) +     
  theme_minimal(base_size = 14)

ggsave("volcano_plot.png", plot = p_graph, width = 10, height = 8, dpi = 300, bg = "white")

# Step 7: PCA Plot (sample clustering)
##########################

# Transform counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Run PCA
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# PCA plot using ggplot
p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  labs(title = "PCA of RNA-seq Samples",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 14)

ggsave("PCA_plot.png", plot = p_pca, width = 8, height = 6, dpi = 300, bg = "white")


# Step 8: Heatmap of Top DE Genes
##########################

# Order results by padj (most significant first)
res_ordered <- res[order(res$padj), ]

# Select top 5 upregulated and top 5 downregulated genes
top_up <- head(res_ordered[res_ordered$log2FoldChange > 0, ], 5)
top_down <- head(res_ordered[res_ordered$log2FoldChange < 0, ], 5)
top_genes <- c(rownames(top_up), rownames(top_down))

# Extract normalized counts for these genes
norm_counts <- assay(vsd)[top_genes, ]

# Add gene symbols if available
rownames(norm_counts) <- ifelse(is.na(res$SYMBOL[top_genes]), 
                                top_genes, 
                                res$SYMBOL[top_genes])

# Plot heatmap
pheatmap(norm_counts,
         scale = "row", 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         annotation_col = coldata,
         main = "Top 5 Up & Down Regulated Genes",
         fontsize_row = 10,
         fontsize_col = 12)

ggsave("heatmap_top_genes.png", width = 8, height = 6, dpi = 300, bg = "white")





