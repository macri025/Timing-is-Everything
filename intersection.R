#=============================
#Interaction model
#=============================
library(DESeq2)
library(limma)
library(ggVennDiagram)

#test timepoint, treatment effect, and intersection
# LRT p-values used to test for any interaction, not just one (eg 3 vs 21)
dds.int <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = coldata,
                                  design = ~ timepoint + treatment + timepoint:treatment)
dds.int <- DESeq(dds.int, test = "LRT", reduced = ~ timepoint + treatment)
res.int <- results(dds.int)  # padj comes from LRT, ignore its log2FoldChange column

# Pull each interaction coefficient's fold change separately (Wald, but only using LFC not p-value)
lfc_zt7  <- results(dds.int, name = "timepointzt_7.treatmentL")$log2FoldChange
lfc_zt17 <- results(dds.int, name = "timepointzt_17.treatmentL")$log2FoldChange
lfc_zt21 <- results(dds.int, name = "timepointzt_21.treatmentL")$log2FoldChange

lfc_mat <- cbind(zt7 = lfc_zt7, zt17 = lfc_zt17, zt21 = lfc_zt21)
rownames(lfc_mat) <- rownames(res.int)

# Max absolute interaction effect per gene, and which timepoint drove it
max_abs_lfc <- apply(lfc_mat, 1, function(x) max(abs(x), na.rm = TRUE))
driver_tp   <- apply(lfc_mat, 1, function(x) c("zt7","zt17","zt21")[which.max(abs(x))])

res.int$max_interaction_lfc <- max_abs_lfc[rownames(res.int)]
res.int$driver_timepoint    <- driver_tp[rownames(res.int)]

res.int <- res.int[order(res.int$padj),]

sig.int <- rownames(res.int)[which(res.int$padj < 0.05 &
                                     res.int$max_interaction_lfc > 1 &
                                     res.int$baseMean > 20)]

#=========================
#Genes with constant (circadian) pattern in controls
#========================

#make a vector of clock times sampled and isolate normalised counts for D (limma=linear)
vst_mat <- assay(vst(dds.int, blind = FALSE))
d.cols <- grep("D", colnames(vst_mat), value = TRUE)
expr.d <- vst_mat[, d.cols]
time <- as.numeric(gsub("zt_", "", colData(dds.int)[d.cols, "timepoint"]))

#24 = 24 hour circadian period
design <- model.matrix(~ sin(2*pi*time/24) + cos(2*pi*time/24))

fit <- lmFit(expr.d, design)
fit <- eBayes(fit)

# F-test on both harmonic terms jointly = test for rhythmicity
res_rhythmic <- topTable(fit, coef = 2:3, number = Inf, sort.by = "none")
sig_rhythmic <- rownames(res_rhythmic)[res_rhythmic$adj.P.Val < 0.05]



#ADD AMPLITUDE THRESHOLD
# Extract sin/cos coefficients per gene from the fitted model
coefs <- coef(fit)  # matrix: genes x [intercept, sin, cos]
beta_sin <- coefs[, 2]
beta_cos <- coefs[, 3]

# Amplitude and phase
amplitude   <- sqrt(beta_sin^2 + beta_cos^2)          # log2 units (half peak-to-trough)
phase_rad   <- atan2(beta_cos, beta_sin)
peak_time   <- (phase_rad / (2*pi)) * 24 %% 24        # approx phase in ZT hours

# Fold-change interpretation: peak-to-trough fold change
p2t_fc <- 2^(2 * amplitude)

# Attach to results table
res_rhythmic$amplitude   <- amplitude[rownames(res_rhythmic)]
res_rhythmic$peak_to_trough_FC <- p2t_fc[rownames(res_rhythmic)]
res_rhythmic$peak_time_ZT <- peak_time[rownames(res_rhythmic)]

# Filter: significant AND a real oscillation size
sig_rhythmic <- rownames(res_rhythmic)[
  res_rhythmic$adj.P.Val < 0.05 &
    res_rhythmic$peak_to_trough_FC > 1.5
]


#=============================
#Find interaction (circadian and light responsive) for candidate genes
#=============================
intersection = intersect(sig.int, sig_rhythmic)
length(intersection) #614

ggVennDiagram(list('Light \nresponsive' = sig.int, 'Rhythmic' = sig_rhythmic))


#GO enrichment
#library(clusterProfiler)
go.enrich.int <- enricher(
  gene         = intersection,
  universe     = filtered.loci,
  TERM2GENE    = goterm2gene,
  TERM2NAME    = goterm2name,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

go.int.df=as.data.frame(go.enrich.int)
dotplot(go.enrich.int, showCategory = 20)

#Find gene names from Bradi IDs
intersection.df = unlist(intersection)
intersection.df <- data.frame(gene = intersection)
intersection.annot <- merge(intersection.df,
                             genes[, c("Gene.Name", "Description")],
                             by.x = "gene",
                             by.y = "Gene.Name",
                             all.x = TRUE)
head(intersection.annot)


























#=============================
#Comparison 1: upregulated at zt17 (upreg.17) and circadian reg in dark samples (sig_rhythmic)
#Comparison 2a: upregulated at zt17 and more induced by light at 17 compared to 7 (interaction test)
#Comparison 2b: upregulated at zt17 and more induced by light at 17 compared to 7 (DE test)
#Comparison 3: upregulated at zt17 and not DE at 7
#===============================



#====================
#Up at zt17 AND circadian influenced (comp1)
#====================
upreg.17 = sig.17[sig.17$log2FoldChange>0,]
upreg.17 = rownames(upreg.17)


  ggVennDiagram(list('upreg 17' = upreg.17, "circadian" = sig_rhythmic)) +
  ggtitle('Comp 1')
comp1 = intersect (upreg.17,sig_rhythmic)

#==========================
#Up at zt17 but not DE at zt 7 (comp3)
#==========================
nonsig.7 = res.7.df[!(rownames(res.7.df) %in% rownames(sig.7)), ]
nonsig.7 = nonsig.7[complete.cases(nonsig.7), ] #remove NAs
nonsig.7 = rownames(nonsig.7)


comp3 = intersect (upreg.17, nonsig.7)
ggVennDiagram(list(zt7.circ.reg = comp1, "zt17 not zt7" = comp3)) +
  ggtitle('Comp 3')


#==========================
#different influence of light at zt7 vs zt17 (interaction test) - comp2a
#==========================
library(DESeq2)
dds = DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ timepoint + treatment + timepoint:treatment 
)

dds = dds[rowSums(counts(dds)) >= 10,]
dds.7.17 = dds[, dds$timepoint %in% c('zt_7', 'zt_17')]

dds.7.17$timepoint <- droplevels(as.factor(dds.7.17$timepoint))
dds.7.17$treatment <- droplevels(as.factor(dds.7.17$treatment))

design(dds.7.17) <- ~ timepoint + treatment + timepoint:treatment
dds.7.17 <- DESeq(dds.7.17, test = "LRT", reduced = ~ timepoint + treatment)

#find sig results
res.int.717 <- as.data.frame(results(dds.7.17))

sig.int.717 <- res.int.717 |>
  subset(!is.na(padj) & padj < 0.05)

sig.int.717.genes <- rownames(sig.int.717)

ggVennDiagram(list('upreg 17' = upreg.17, "DE at 17 vs 7" = sig.int.717.genes))+
  ggtitle('Comp 2a')

comp2a = intersect(upreg.17, sig.int.717.genes)


#different influence of light at zt7 vs zt17 (DE test) - comp2b

de.7.17.samples = coldata$treatment == "L" & coldata$timepoint %in% c("zt_17", "zt_7") #subset samples of interest

counts7.17  = counts[, de.7.17.samples]
coldata7.17 = coldata[de.7.17.samples, ]


coldata7.17$timepoint = droplevels(factor(coldata7.17$timepoint)) #remove old factors

# Set zt_7 as the reference level so results are 17L vs 7L
coldata7.17$timepoint <- relevel(coldata7.17$timepoint, ref = "zt_7")

# Build DESeq dataset
dds.7.17 <- DESeqDataSetFromMatrix(
  countData = counts7.17,
  colData   = coldata7.17,
  design    = ~ timepoint
)

#Filter low count genes
dds.7.17 <- dds.7.17[rowSums(counts(dds.7.17)) >= 10, ]

# Run DESeq2
dds.7.17 <- DESeq(dds.7.17)

# Extract results: 17L vs 7L
res.7.17 = results(dds.7.17, contrast = c("timepoint", "zt_17", "zt_7"))
res.7.17 = res.7.17[order(res.7.17$padj), ]

sig.7.17 = res.7.17[!is.na(res.7.17$padj) & res.7.17$padj < 0.05, ]
up.17.vs.7 = sig.7.17[sig.7.17$log2FoldChange > 0, ]
up.17.vs.7.genes = rownames(up.17.vs.7)

ggVennDiagram(list('upreg 17' = upreg.17, "DE at 17 vs 7" = up.17.vs.7.genes))+
  ggtitle('Comp 2b')

comp2b = intersect(upreg.17, up.17.vs.7.genes)

ggVennDiagram(list('interaction' = comp2a, "DESeq" = comp2b))+
  ggtitle('Interaction vs DESeq')

#======================
#overall comparisons
#=======================

ggVennDiagram(list('up at 17 and circ' = comp1, "DE at 17 vs 7" = comp2a, "up at 17 not at 7" = comp3))+
  ggtitle('Overall comp with 2a')

ggVennDiagram(list('up at 17 and circ' = comp1, "DE at 17 vs 7" = comp2b, "up at 17 not at 7" = comp3))+
  ggtitle('Overall comp with 2b')

full.comp.2a =  Reduce(intersect, list(comp1, comp2a, comp3))
full.comp.2b =  Reduce(intersect, list(comp1, comp2b, comp3))

all.comp <- list(comp1, comp2a, comp2b, comp3)
complist <- table(unlist(lapply(all.comp, unique)))
sum(complist >= 2) #316

