#Libraries
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)

#===================================================
#ZT3
#===================================================
#DE expression between treatment and control
zt3.dat = coldata$timepoint == "zt_3"

zt3.counts  = counts[, zt3.dat]
zt3.conds = coldata[zt3.dat, ]

dds.3 = DESeqDataSetFromMatrix(
  countData = zt3.counts,
  colData   = zt3.conds,
  design    = ~ treatment
)

dds.3 = DESeq(dds.3)

#save significant genes and remove NAs
res.3 = results(dds.3, contrast = c("treatment", "L", "D"))
sig.3 = res.3[which(res.3$padj < 0.05 & abs(res.3$log2FoldChange) > 1 
                    & res.3$baseMean > 20), ] #filter by biological, statistical sig and min counts
sig.3 = sig.3[!is.na(sig.3$padj), ]
nrow(sig.3) #1375

#make into lists for later comparison (match gene names)
sig.3.genes = rownames(sig.3)

#add gene names
res.3.df = as.data.frame(res.3)
res.3.df$gene = rownames(res.3.df)

res.3.annotated <- merge(res.3.df,
                         genes[, c("Gene.Name", "Description")],
                         by.x = "gene",
                         by.y = "Gene.Name",
                         all.x = TRUE)



#===================================================
#ZT7
#===================================================
#DE expression between treatment and control
zt7.dat = coldata$timepoint == "zt_7"

zt7.counts  = counts[, zt7.dat]
zt7.conds = coldata[zt7.dat, ]

dds.7 = DESeqDataSetFromMatrix(
  countData = zt7.counts,
  colData   = zt7.conds,
  design    = ~ treatment
)

dds.7 = DESeq(dds.7)

#save significant genes and remove NAs
res.7 = results(dds.7, contrast = c("treatment", "L", "D"))
sig.7 = res.7[which(res.7$padj < 0.05 & abs(res.7$log2FoldChange) > 1 & res.7$baseMean > 20), ] #filter by biological, statistical sig and min counts
sig.7 = sig.7[!is.na(sig.7$padj), ]
nrow(sig.7) #3153

#make into vector for later comparison (just gene names)
sig.7.genes = rownames(sig.7)

#add gene names
res.7.df = as.data.frame(res.7)
res.7.df$gene = rownames(res.7.df)

res.7.annotated <- merge(res.7.df,
                         genes[, c("Gene.Name", "Description")],
                         by.x = "gene",
                         by.y = "Gene.Name",
                         all.x = TRUE)

#===================================================
#ZT17
#===================================================
#DE expression between treatment and control
zt17.dat = coldata$timepoint == "zt_17"

zt17.counts  = counts[, zt17.dat]
zt17.conds = coldata[zt17.dat, ]

dds.17 = DESeqDataSetFromMatrix(
  countData = zt17.counts,
  colData   = zt17.conds,
  design    = ~ treatment
)

dds.17 = DESeq(dds.17)

#save significant genes and remove NAs
res.17 = results(dds.17, contrast = c("treatment", "L", "D"))
sig.17 = res.17[which(res.17$padj < 0.05 & abs(res.17$log2FoldChange) > 1 & res.17$baseMean > 20), ] #filter by biological, statistical sig and min counts
sig.17 = sig.17[!is.na(sig.17$padj), ]
nrow(sig.17) #2131

#make into vector for later comparison (just gene names)
sig.17.genes = rownames(sig.17)

#add gene names
res.17.df = as.data.frame(res.17)
res.17.df$gene = rownames(res.17.df)

res.17.annotated <- merge(res.17.df,
                         genes[, c("Gene.Name", "Description")],
                         by.x = "gene",
                         by.y = "Gene.Name",
                         all.x = TRUE)

#===================================================
#ZT21
#===================================================
#DE expression between treatment and control
zt21.dat = coldata$timepoint == "zt_21"

zt21.counts  = counts[, zt21.dat]
zt21.conds = coldata[zt21.dat, ]

dds.21 = DESeqDataSetFromMatrix(
  countData = zt21.counts,
  colData   = zt21.conds,
  design    = ~ treatment
)

dds.21 = DESeq(dds.21)

#save significant genes and remove NAs
res.21 = results(dds.21, contrast = c("treatment", "L", "D"))
sig.21 = res.21[which(res.21$padj < 0.05 & abs(res.21$log2FoldChange) > 1 
                    & res.21$baseMean > 20), ] #filter by biological, statistical sig and min counts
sig.21 = sig.21[!is.na(sig.21$padj), ]
nrow(sig.21) #2104

#make into vector for later comparison (just gene names)
sig.21.genes = rownames(sig.21)

#add gene names
res.21.df = as.data.frame(res.21)
res.21.df$gene = rownames(res.21.df)

res.21.annotated <- merge(res.21.df,
                          genes[, c("Gene.Name", "Description")],
                          by.x = "gene",
                          by.y = "Gene.Name",
                          all.x = TRUE)

