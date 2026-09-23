#Libraries
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)

#===================================================
#Volcano plots
#===================================================
#ZT3
res.3.df$sig <- "Not sig"
res.3.df$sig[res.3.df$padj < 0.05 & res.3.df$log2FoldChange > 1  & res.3.df$baseMean > 20] <- "Up"
res.3.df$sig[res.3.df$padj < 0.05 & res.3.df$log2FoldChange < -1 & res.3.df$baseMean > 20] <- "Down"

v3 = ggplot(res.3.df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not sig" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  )
labs(
  x = "log2 Fold Change",
  y = "-log10 adjusted p-value",
  color = "Significance"
)

# most significant gene overall
res.3.annotated[which.min(res.3.annotated$padj), ]

# top 10 by significance
head(res.3.annotated[order(res.3.annotated$padj), ], 10)




#ZT7
res.7.df$sig <- "Not sig"
res.7.df$sig[res.7.df$padj < 0.05 & res.7.df$log2FoldChange > 1  & res.7.df$baseMean > 20] <- "Up"
res.7.df$sig[res.7.df$padj < 0.05 & res.7.df$log2FoldChange < -1 & res.7.df$baseMean > 20] <- "Down"

v7 = ggplot(res.7.df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not sig" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  )
labs(
  x = "log2 Fold Change",
  y = "-log10 adjusted p-value",
  color = "Significance"
)

# most significant gene overall
res.7.annotated[which.min(res.7.annotated$padj), ]

# top 10 by significance
head(res.7.annotated[order(res.7.annotated$padj), ], 10)





#ZT17
res.17.df$sig <- "Not sig"
res.17.df$sig[res.17.df$padj < 0.05 & res.17.df$log2FoldChange > 1  & res.17.df$baseMean > 20] <- "Up"
res.17.df$sig[res.17.df$padj < 0.05 & res.17.df$log2FoldChange < -1 & res.17.df$baseMean > 20] <- "Down"

v17 = ggplot(res.17.df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not sig" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  )
labs(
  x = "log2 Fold Change",
  y = "-log10 adjusted p-value",
  color = "Significance"
)

# most significant gene overall
res.17.annotated[which.min(res.17.annotated$padj), ]

# top 10 by significance
head(res.17.annotated[order(res.17.annotated$padj), ], 10)








#ZT21
res.21.df = as.data.frame(res.21)
res.21.df$gene = rownames(res.21.df)

res.21.df$sig <- "Not sig"
res.21.df$sig[res.21.df$padj < 0.05 & res.21.df$log2FoldChange > 1  & res.21.df$baseMean > 20] <- "Up"
res.21.df$sig[res.21.df$padj < 0.05 & res.21.df$log2FoldChange < -1 & res.21.df$baseMean > 20] <- "Down"

v21 = ggplot(res.21.df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not sig" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  )
labs(
  x = "log2 Fold Change",
  y = "-log10 adjusted p-value",
  color = "Significance"
)

# most significant gene overall
res.21.annotated[which.min(res.21.annotated$padj), ]

# top 10 by significance
head(res.21.annotated[order(res.21.annotated$padj), ], 10)
head(res.21.annotated[order(res.21.annotated$log2FoldChange), ], 10)












#========================
#PCA plots/loadings
#========================
#ZT3
vsd3 = vst(dds.3)

pca3 = plotPCA(vsd3, intgroup = "treatment", returnData = TRUE)

pca3$sample = rownames(pca3)

pcplot3 = ggplot(pca3, aes(PC1, PC2, color = treatment, label = sample)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(100 * attr(pca3, "percentVar")[1]), "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pca3, "percentVar")[2]), "% variance")) +
  theme_bw() +
  scale_color_manual(values = c("D" = "#6E6A6A", "L" = "#F5B800"))




#ZT7
vsd7 = vst(dds.7)

pca7 = plotPCA(vsd7, intgroup = "treatment", returnData = TRUE)

pca7$sample = rownames(pca7)

pcplot7 = ggplot(pca7, aes(PC1, PC2, color = treatment, label = sample)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(100 * attr(pca7, "percentVar")[1]), "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pca7, "percentVar")[2]), "% variance")) +
  theme_bw() +
  scale_color_manual(values = c("D" = "#6E6A6A", "L" = "#F5B800"))





#ZT17
vsd17 = vst(dds.17)

pca17 = plotPCA(vsd17, intgroup = "treatment", returnData = TRUE)

pca17$sample = rownames(pca17)

pcplot17 = ggplot(pca17, aes(PC1, PC2, color = treatment, label = sample)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(100 * attr(pca17, "percentVar")[1]), "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pca17, "percentVar")[2]), "% variance")) +
  theme_bw() +
  scale_color_manual(values = c("D" = "#6E6A6A", "L" = "#F5B800"))






#ZT21
vsd21 = vst(dds.21)

pca21 = plotPCA(vsd21, intgroup = "treatment", returnData = TRUE)

pca21$sample = rownames(pca21)

pcplot21 = ggplot(pca21, aes(PC1, PC2, color = treatment, label = sample)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(100 * attr(pca21, "percentVar")[1]), "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pca21, "percentVar")[2]), "% variance")) +
  theme_bw() +
  scale_color_manual(values = c("D" = "#6E6A6A", "L" = "#F5B800"))



#LOADINGS
get.pc.loadings = function(vsd, pc = "PC1", ntop = 500, top.n = 50) {
  mat = assay(vsd)
  select = order(rowVars(mat), decreasing = TRUE)[seq_len(min(ntop, nrow(mat)))]
  loadings = prcomp(t(mat[select, ]))$rotation[, pc]
  
  loadings.df = data.frame(gene = names(loadings), loading = loadings, abs.loading = abs(loadings))
  loadings.df = loadings.df[order(loadings.df$abs.loading, decreasing = TRUE)[1:top.n], ]
  
  merge(loadings.df, genes[, c("Gene.Name", "Description")],
        by.x = "gene", by.y = "Gene.Name", all.x = TRUE, sort = FALSE)
}

#Run on samples
loadings.3 = get.pc.loadings(vsd3, pc = "PC1", top.n = 50)
loadings.7 = get.pc.loadings(vsd7, pc = "PC1", top.n = 50)
loadings.17 = get.pc.loadings(vsd17, pc = "PC1", top.n = 50)
loadings.21 = get.pc.loadings(vsd21, pc = "PC1", top.n = 50)

#Combine into one for comparison
combined.loadings = rbind(
  cbind(loadings.3, timepoint = "zt3"),
  cbind(loadings.7, timepoint = "zt7"),
  cbind(loadings.17, timepoint = "zt17"),
  cbind(loadings.21, timepoint = "zt21")
)

#which genes are driving more than one timepoint
combined.loadings[combined.loadings$timepoint %in% c("zt21", "zt17") &
                    combined.loadings$gene %in% intersect(loadings.21$gene, loadings.17$gene), ]


