#=============
# Gene x timepoint matrix
#=============
library(DESeq2)
library(dplyr)
library(tibble)

#extract logFC and p value at each timepoint (L vs D) — used only to define light-regulated genes
get.res <- function(dds, tp) {
  res <- results(dds, contrast = c("treatment", "L", "D"))
  as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    transmute(gene, !!paste0("LFC_", tp) := log2FoldChange,
              !!paste0("padj_", tp) := padj)
}

contrast.list <- Map(get.res, dds.list, names(dds.list))
lfc.table <- Reduce(function(x, y) full_join(x, y, by = "gene"), contrast.list)


#Pull significant (light-regulated) genes
padj.cols <- grep("^padj_", colnames(lfc.table), value = TRUE)
lfc.cols  <- grep("^LFC_",  colnames(lfc.table), value = TRUE)

sig.genes <- lfc.table %>%
  filter(if_any(all_of(padj.cols), ~ !is.na(.x) & .x < 0.05) ) %>%
  filter(if_any(all_of(lfc.cols), ~ abs(.x) > 1)) %>%
  pull(gene)

length(sig.genes) #5789

#=============
# Light-only expression trajectory (for clustering shape)
#=============

# Pull normalized counts from LIGHT samples only, average replicates per timepoint
get.light.expr <- function(dds, tp) {
  norm.counts <- counts(dds, normalized = TRUE)
  l.samples <- colnames(dds)[colData(dds)$treatment == "L"]
  l.counts  <- norm.counts[, l.samples, drop = FALSE]
  
  log.expr <- log2(l.counts + 1)
  data.frame(gene = rownames(log.expr),
             mean_log_expr = rowMeans(log.expr)) %>%
    dplyr::rename(!!paste0("LFC_", tp) := mean_log_expr)
}

expr.list <- Map(get.light.expr, dds.list, names(dds.list))
expr.table <- Reduce(function(x, y) full_join(x, y, by = "gene"), expr.list)

expr.cols <- grep("^LFC_", colnames(expr.table), value = TRUE)

#Build matrix: genes x timepoints, values = mean log2(normalized expr) in light
response.mat <- expr.table %>%
  filter(gene %in% sig.genes) %>%
  column_to_rownames("gene") %>%
  dplyr::select(all_of(expr.cols)) %>%
  as.matrix()

colnames(response.mat) <- c("ZT3", "ZT7", "ZT17", "ZT21")  # clean names

# Drop genes with any NA (failed test / low counts at a timepoint)
response.mat <- response.mat[complete.cases(response.mat), ]

# Z-score each gene's row across the 4 timepoints (captures shape, not magnitude)
z.mat <- t(scale(t(response.mat)))
z.mat <- z.mat[complete.cases(z.mat), ]  # drop genes with zero variance

dim(z.mat) #5787 x 4







# ============================================================
# Mfuzz (fuzzy c-means) — short time courses
# ============================================================

library(Mfuzz)
library(Biobase)

# Build MFuzz expression set from non z-score matrix, does own standardisation
eset <- ExpressionSet(assayData = response.mat)
eset <- filter.std(eset, min.std = 0, visu = FALSE)
eset.std <- standardise(eset)

# Estimate fuzzifier m
m.est <- mestimate(eset.std)
m.est

# Choose number of clusters: check within-cluster distance for a
# range of k, look for an elbow (Dmin is Mfuzz's built-in helper)
par(cex.lab = 1.4, cex.axis = 1.3, cex.main = 1.4, cex = 1.3)
Dmin(eset.std, m = m.est, crange = 2:12, repeats = 10, visu = TRUE)

# k picked by elbow and stability (5)
set.seed(42)
mfuzz.res <- mfuzz(eset.std, c = 5, m = m.est)

# Plot cluster trajectories
timepoints <- colnames(z.mat)
mfuzz.plot2(eset.std, cl = mfuzz.res, mfrow = c(2, 4),
            time.labels = timepoints, x11 = FALSE,
            cex.main = 1.7,   # cluster title size
            cex.lab = 1.5,    # axis title size (x/y labels)
            cex.axis = 1.5)   # tick label size

# Extract hard cluster assignment (highest membership) and
# membership score (how well each gene fits its cluster, 0-1)
gene.clusters.fuzz <- data.frame(
  gene = names(mfuzz.res$cluster),
  cluster = mfuzz.res$cluster,
  membership = apply(mfuzz.res$membership, 1, max)
)

# identify "core" members of each cluster (higher confidence)
gene.clusters.core <- gene.clusters.fuzz %>% filter(membership > 0.5)

table(gene.clusters.fuzz$cluster)
table(gene.clusters.core$cluster)


## ============================================================
## OPTION B: Hierarchical clustering, correlation distance
## ============================================================

# Distance = 1 - Pearson correlation between gene response profiles
# (captures shape/timing of response, ignores overall magnitude)
cor.mat  <- cor(t(z.mat), method = "pearson")
dist.cor <- as.dist(1 - cor.mat)

hc <- hclust(dist.cor, method = "average")  # average linkage

# Look at dendrogram to choose number of clusters
plot(hc, labels = FALSE, main = "Hierarchical clustering (1 - Pearson correlation)")
rect.hclust(hc, k = 5, border = "red")
rect.hclust(hc, k = 9, border = "blue")

#loop to test which number of clusters provides the best silhouette width
library(cluster)
sil.widths <- sapply(2:12, function(k) {
  cl <- cutree(hc, k = k)
  mean(silhouette(cl, dist.cor)[, "sil_width"])
})

plot(2:12, sil.widths, type = "b", xlab = "k", ylab = "avg silhouette width")

clusters.hc <- cutree(hc, k = 9) 

gene.clusters.hc <- data.frame(
  gene = rownames(z.mat),
  cluster = clusters.hc
)

table(gene.clusters.hc$cluster)



#============================================================
#Mean trajectory per cluster
#============================================================
timepoints <- colnames(z.mat)
plot.cluster.trajectories <- function(z.mat, cluster.assignments, timepoints) {
  df <- as.data.frame(z.mat) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "timepoint", values_to = "z_LFC") %>%
    mutate(timepoint = factor(timepoint, levels = timepoints)) %>%
    left_join(cluster.assignments, by = "gene")   # <- fixed: dot, not underscore
  
  summary_df <- df %>%
    group_by(cluster, timepoint) %>%
    summarise(mean_z = mean(z_LFC), se = sd(z_LFC) / sqrt(n()), .groups = "drop")
  
  ggplot(summary_df, aes(x = timepoint, y = mean_z, group = cluster, color = factor(cluster))) +
    geom_line(linewidth = 1) +
    geom_errorbar(aes(ymin = mean_z - se, ymax = mean_z + se), width = 0.1) +
    facet_wrap(~cluster) +
    theme_minimal() +
    labs(y = "z-scored light-vs-dark LFC", x = "Zeitgeber time", color = "Cluster")
}

# For Mfuzz result:
plot.cluster.trajectories(z.mat, gene.clusters.fuzz[, c("gene", "cluster")], timepoints) +
  ggtitle("Mfuzz cluster trajectories")
#Cluster 5 upreg at zt17

# For hclust result:
plot.cluster.trajectories(z.mat, gene.clusters.hc[, c("gene", "cluster")], timepoints) +
  ggtitle("Hierarchical clustering trajectories")
#Clusters 4, 5, 6, 8

plot.cluster.trajectories(z.mat, gene.clusters.hc[, c("gene", "cluster")], timepoints) +
  ggtitle("Hierarchical clustering trajectories") +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),   # facet labels, e.g. "Cluster 4"
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
#=========================================================
#=========================================================
#Comparisons with other lists - Mfuzz
cluster5.genes.fuzz = gene.clusters.fuzz %>%
  filter(cluster == 5) %>%
  pull(gene)

#Using interaction test
intersect(full.comp.2a, cluster5.genes.fuzz) #4


#incorporating DESeq test
Reduce(intersect, list(full.comp.2a, cluster5.genes.fuzz, comp2b)) #4

#DESeq test only
Reduce(intersect, list(full.comp.2b, cluster5.genes.fuzz)) #18

intersect(sig_rhythmic, cluster5.genes.fuzz) #412
intersect(comp1, cluster5.genes.fuzz) #65

#Limit to core genes and retest all
cluster5.genes.core = gene.clusters.core %>%
  filter(cluster == 5) %>%
  pull(gene)

intersect(full.comp.2a, cluster5.genes.core) #2


Reduce(intersect, list(full.comp.2a, cluster5.genes.core, comp2b)) #2


intersect(sig_rhythmic, cluster5.genes.core) #261
intersect(comp1, cluster5.genes.core) #36

#==============================================================
#Comparisons with other lists - HC

#CLUSTER 4
cluster4.genes.hc = gene.clusters.hc %>%
  filter(cluster == 4) %>%
  pull(gene)

#rhythmic
intersect(sig_rhythmic, cluster4.genes.hc) #123
intersect(comp1, cluster4.genes.hc) #7

#comp1
intersect(comp1, cluster4.genes.hc) #7

#Using interaction test
intersect(full.comp.2a, cluster4.genes.hc) #1

#incorporating DESeq test
Reduce(intersect, list(full.comp.2a, cluster4.genes.hc, comp2b)) #1

#DESeq only
Reduce(intersect, list(full.comp.2b, cluster4.genes.hc)) #2



#CLUSTER 5
cluster5.genes.hc = gene.clusters.hc %>%
  filter(cluster == 5) %>%
  pull(gene)

#rhythmic
intersect(sig_rhythmic, cluster5.genes.hc) #26

#comp1
intersect(comp1, cluster5.genes.hc) #3

#Using interaction test
intersect(full.comp.2a, cluster5.genes.hc) #0

#incorporating DESeq test
Reduce(intersect, list(full.comp.2a, cluster5.genes.hc, comp2b)) #0

#DESeq only
Reduce(intersect, list(full.comp.2b, cluster5.genes.hc)) #0




#CLUSTER 6
cluster6.genes.hc = gene.clusters.hc %>%
  filter(cluster == 6) %>%
  pull(gene)

#rhythmic
intersect(sig_rhythmic, cluster6.genes.hc) #302

#comp1
intersect(comp1, cluster6.genes.hc) #41

#Using interaction test
intersect(full.comp.2a, cluster6.genes.hc) #3

#incorporating DESeq test
Reduce(intersect, list(full.comp.2a, cluster6.genes.hc, comp2b)) #3

#DESeq only
Reduce(intersect, list(full.comp.2b, cluster6.genes.hc)) #13






#CLUSTER 8
cluster8.genes.hc = gene.clusters.hc %>%
  filter(cluster == 8) %>%
  pull(gene)

#rhythmic
intersect(sig_rhythmic, cluster8.genes.hc) #56

#comp1
intersect(comp1, cluster8.genes.hc) #18

#Using interaction test
intersect(full.comp.2a, cluster8.genes.hc) #0

#incorporating DESeq test
Reduce(intersect, list(full.comp.2a, cluster8.genes.hc, comp2b)) #0

#DESeq only
Reduce(intersect, list(full.comp.2b, cluster8.genes.hc)) #3










