#========================
#Pearson correlation
#========================
library(limma)
library(edgeR)

#Check column names same in counts and coldata
stopifnot(all(colnames(counts) == rownames(coldata)))

# ---- Normalize ----
dge.GI <- DGEList(counts = counts)
keep.GI <- filterByExpr(dge.GI, group = interaction(coldata$treatment, coldata$timepoint)) #remove noisy genes and use smallest group to decide sample count threshold (conservative)
dge.GI <- dge.GI[keep.GI, , keep.lib.sizes = FALSE]
dge.GI <- calcNormFactors(dge.GI)

design.GI <- model.matrix(~ timepoint + treatment, data = coldata)
v.GI <- voom(dge.GI, design.GI)   # v$E = normalized log2-CPM, genes x samples

# ---- Pull out GI's expression ----
GI_expr <- v.GI$E['Bradi2g05226', ]

# ---- Residualize against timepoint + treatment ----
fit.GI <- lmFit(v.GI$E, design.GI)
resid_expr.GI <- residuals(fit.GI, v.GI$E)
GI_resid <- resid_expr.GI['Bradi2g05226', ]

# ---- Correlate every gene's residuals against GI's residuals ----
cor_results.GI <- t(apply(resid_expr.GI, 1, function(g) {
  test <- cor.test(g, GI_resid, method = "spearman")
  c(rho = unname(test$estimate), p = test$p.value)
}))

cor_results.GI <- as.data.frame(cor_results.GI)
cor_results.GI$padj <- p.adjust(cor_results.GI$p, method = "BH")
cor_results.GI <- cor_results.GI[order(cor_results.GI$padj), ]


cor_results.GI <- cor_results.GI %>%
  tibble::rownames_to_column("Gene.Name") %>%
  left_join(genes, by = "Gene.Name") %>%
  arrange(padj)


head(cor_results.GI, 20)

#Look at expression across timepoints
library(ggplot2)
library(patchwork)

top_genes <- c("Bradi1g76580", "Bradi2g47540")

# confirm both exist before proceeding
top_genes %in% rownames(resid_expr.GI)

# ---- Build plotting dataframe using the correct objects ----
plot_df <- data.frame(
  sample = colnames(v.GI$E),
  GI = GI_resid,
  treatment = coldata$treatment,
  timepoint = coldata$timepoint
)

for (g in top_genes) {
  plot_df[[g]] <- resid_expr.GI[g, ]
}

head(plot_df)   # sanity check columns are all there

# ---- Scatter plots ----
make_scatter <- function(gene) {
  ggplot(plot_df, aes(x = GI, y = .data[[gene]], color = treatment, shape = timepoint)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
    labs(title = gene, x = "GI (residual expr)", y = paste(gene, "(residual expr)")) +
    theme_minimal()
}

plots <- lapply(top_genes, make_scatter)
wrap_plots(plots, ncol = 2)



#====================
#GENIE3 - AI
#====================
library(GENIE3)

exprMat.GI <- v.GI$E

# Or run against all genes as potential regulators (not isolating TFs)
#weightMat.GI <- GENIE3(exprMat.GI)

# GI as regulator: how well GI predicts every other gene
#GI_as_regulator <- weightMat.GI['Bradi2g05226', ]
#GI_as_regulator <- sort(GI_as_regulator, decreasing = TRUE)

GI_as_regulator <- readRDS("//storage.hcs-p01.otago.ac.nz/biochemistry/Lab_Groups/brownfieldlab/Documents/Riley/GI_as_regulator.rds")
head(GI_as_regulator, 20)
summary(GI_as_regulator)
hist(GI_as_regulator, breaks = 100)

# GI as target: how well every other gene predicts GI
#GI_as_target <- weightMat.GI[, 'Bradi2g05226']
#GI_as_target <- sort(GI_as_target, decreasing = TRUE)

GI_as_target <- readRDS("//storage.hcs-p01.otago.ac.nz/biochemistry/Lab_Groups/brownfieldlab/Documents/Riley/GI_as_target.rds")
head(GI_as_target, 20)
summary(GI_as_target)
hist(GI_as_target, breaks = 100)







#======================
#Comparison
#======================
cor_hits <- rownames(cor_results.GI)[cor_results.GI$padj < 0.05]
cor_hits <- setdiff(cor_hits, "Bradi2g05226")  # drop GI self-hit

top_n <- 50   # cast a slightly wide net on the GENIE3 side
top_reg <- names(sort(GI_as_regulator, decreasing = TRUE))[1:top_n]
top_tgt <- names(sort(GI_as_target, decreasing = TRUE))[1:top_n]

overlap_as_regulator <- intersect(cor_hits, top_reg)
overlap_as_target <- intersect(cor_hits, top_tgt)

overlap_as_regulator #0
overlap_as_target #0


gene_ranks <- function(genes) {
  idx <- match(genes, rownames(cor_results.GI))
  
  data.frame(
    gene = genes,
    reg_rank = rank(-GI_as_regulator)[genes],
    tgt_rank = rank(-GI_as_target)[genes],
    rho = cor_results.GI$rho[idx],
    padj = cor_results.GI$padj[idx],
    row.names = NULL
  )
}

gene_ranks(cor_hits)
gene_ranks(top.candidates)
gene_ranks(tiptop.candidates)
