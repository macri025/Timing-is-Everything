#============================================================
#Venn diagram
#============================================================
#install.packages("ggVennDiagram")
library(ggVennDiagram)

ggVennDiagram(list(zt3 = sig.3.genes, zt7 = sig.7.genes, zt17 = sig.17.genes, zt21 = sig.21.genes)) +
  scale_fill_gradient(low = "white", high = "steelblue")

ggVennDiagram(list(zt3 = sig.3.genes,  zt17 = sig.17.genes)) +
  scale_fill_gradient(low = "white", high = "steelblue")

ggVennDiagram(list(zt7 = sig.7.genes,  zt17 = sig.17.genes)) +
  scale_fill_gradient(low = "white", high = "steelblue")

#Lists of DE genes only during day or only during night (Bradi)
DE.day = setdiff(intersect(sig.3.genes, sig.7.genes), union(sig.17.genes, sig.21.genes))
DE.night = setdiff(intersect(sig.17.genes, sig.21.genes), union(sig.3.genes, sig.7.genes))

length(DE.day) #should match Venn
length(DE.night)

#Bradi -> gene names
DE.day.names = genes$Description[match(DE.day, genes$Gene.Name)]
DE.night.names = genes$Description[match(DE.night, genes$Gene.Name)]

#==============================================================
#KEGG and GO analysis - zt3+7  vs zt17+21
#==============================================================
#BiocManager::install("clusterProfiler")
#BiocManager::install("readr")
library(clusterProfiler)
library(dplyr)
library(readr)

annot = read_tsv("annotation_info.txt") #has KO IDs
annot.sub <- annot[annot$locusName %in% genes$Gene.Name, ] #only genes in genes

#collapse duplicate transcripts

sum(duplicated(annot$locusName)) #24408
length(unique(annot$locusName)) #32439

annot.sub <- annot.sub %>%
  group_by(locusName) %>%
  summarise(
    across(c(KO, GO, Pfam, Panther, KOG, ec,
             best_arabi_gene, best_arabi_defline,
             best_rice_gene, best_rice_defline),
           ~ paste(unique(na.omit(.x)), collapse = "; ")),
    .groups = "drop"
  ) %>%
  mutate(across(everything(), ~ na_if(.x, "")))

sum(duplicated(annot.sub$locusName)) #0

sum(is.na(annot.sub$GO)) #12441/29764
sum(is.na(annot.sub$KO)) #21835/29764


#get KO mapping from KEGG
ko2pathway <- download_KEGG("ko")
term2gene <- ko2pathway$KEGGPATHID2EXTID
term2name <- ko2pathway$KEGGPATHID2NAME

# Filter annotation to genes significant at each timepoint
annot.am <- annot.sub %>% filter(locusName %in% DE.day)
annot.pm <- annot.sub %>% filter(locusName %in% DE.night)

# Check all IDs match
setdiff(DE.day, annot.sub$locusName)
setdiff(DE.night, annot.sub$locusName)

# Pull KO terms, drop blanks, split any multi-KO entries
DE.day.ko <- annot.am %>% filter(!is.na(KO) & KO != "") %>% pull(KO)
DE.day.ko <- unlist(strsplit(DE.day.ko, "[ ,]+"))

DE.night.ko <- annot.pm %>% filter(!is.na(KO) & KO != "") %>% pull(KO)
DE.night.ko <- unlist(strsplit(DE.night.ko, "[ ,]+"))

#run compare cluster
am.pm.lists <- list('Daytime timepoints' = DE.day.ko,
                   'Nighttime timepoints'= DE.night.ko)

ck1 <- compareCluster(geneCluster = am.pm.lists,
                     fun = "enrichKEGG",
                     organism = "ko",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH")

as.data.frame(ck1)
dotplot(ck1, showCategory = 15)+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8))
ck1[grepl("Photosynthesis", ck.result$Description, ignore.case = TRUE), ]
#==============================================================
#KEGG and GO analysis - zt7 vs zt 17
#==============================================================
#run DESeq between light at zt7 and light at zt17
library(DESeq2)

#Identify which samples belong to 7L and 17L
L.only.7.17 <- rownames(coldata)[
  grepl("^17L", rownames(coldata)) | grepl("^7L", rownames(coldata))
]

# 2. Subset counts and coldata to just those samples
counts.7.17  <- counts[, L.only.7.17]
coldata.7.17 <- coldata[L.only.7.17, , drop = FALSE]

# 3. Create a "line" grouping factor if you don't already have one
coldata.7.17$line <- ifelse(grepl("^17", rownames(coldata.7.17)), "17", "7")
coldata.7.17$line <- factor(coldata.7.17$line, levels = c("7", "17"))  # "7" = reference level

# 4. Build DESeqDataSet
dds.7.17 <- DESeqDataSetFromMatrix(
  countData = counts.7.17,
  colData   = coldata.7.17,
  design    = ~ line
)

dds.7.17 <- DESeq(dds.7.17)

# 7. Extract results: 17 vs 7
res.7.17 <- results(dds.7.17, contrast = c("line", "17", "7"))
sig.7.17 = res.7.17[which(res.7.17$padj < 0.05 & abs(res.7.17$log2FoldChange) > 1 
                    & res.7.17$baseMean > 20), ] #filter by biological, statistical sig and min counts
sig.7.17 = sig.7.17[!is.na(sig.7.17$padj), ]
sig.7.17 = rownames(sig.7.17)
length(sig.7.17) #1556

#KO analysis
# Filter annotation to genes significant at each timepoint
annot.7.17 <- annot %>% filter(locusName %in% sig.7.17)

# Sanity check - did all IDs match? (catches typos, case mismatches, .1 suffixes etc.)
setdiff(sig.7.17, annot$locusName)

# Pull KO terms, drop blanks, split any multi-KO entries
sig.7.17.ko <- annot.7.17 %>% filter(!is.na(KO) & KO != "") %>% pull(KO)
sig.7.17.ko <- unlist(strsplit(sig.7.17.ko, "[ ,]+"))

# Run enrichment
enrich.7.17 <- enricher(
  gene          = sig.7.17.ko,      # vector of KO terms
  TERM2GENE     = term2gene,
  TERM2NAME     = term2name,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# View results
head(as.data.frame(enrich.7.17))

# Dot plot of enriched pathways
dotplot(enrich.7.17, showCategory = 20)+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8))











#==================================================
#Individual KEGG analysis - zt3
#==================================================
#run compare cluster
# Filter annotation to genes significant at each timepoint
annot.3 <- annot %>% filter(locusName %in% sig.3.genes)

# Pull KO terms, drop blanks, split any multi-KO entries
sig.3.ko <- annot.3 %>% filter(!is.na(KO) & KO != "") %>% pull(KO)
sig.3.ko <- unlist(strsplit(sig.3.ko, "[ ,]+"))

enrich.3 <- enricher(
  gene          = sig.3.ko,      # vector of KO terms
  TERM2GENE     = term2gene,
  TERM2NAME     = term2name,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2)

dotplot(enrich.3, showCategory = 20)




#==================================================
#Individual KEGG analysis - zt7
#==================================================
#run compare cluster
# Filter annotation to genes significant at each timepoint
annot.7 <- annot %>% filter(locusName %in% sig.7.genes)

# Pull KO terms, drop blanks, split any multi-KO entries
sig.7.ko <- annot.7 %>% filter(!is.na(KO) & KO != "") %>% pull(KO)
sig.7.ko <- unlist(strsplit(sig.7.ko, "[ ,]+"))

enrich.7 <- enricher(
  gene          = sig.7.ko,      # vector of KO terms
  TERM2GENE     = term2gene,
  TERM2NAME     = term2name,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2)

dotplot(enrich.7, showCategory = 20)


#==================================================
#Individual KEGG analysis - zt17
#==================================================
#run compare cluster
# Filter annotation to genes significant at each timepoint
annot.17 <- annot %>% filter(locusName %in% sig.17.genes)

# Pull KO terms, drop blanks, split any multi-KO entries
sig.17.ko <- annot.17 %>% filter(!is.na(KO) & KO != "") %>% pull(KO)
sig.17.ko <- unlist(strsplit(sig.17.ko, "[ ,]+"))

enrich.17 <- enricher(
  gene          = sig.17.ko,      # vector of KO terms
  TERM2GENE     = term2gene,
  TERM2NAME     = term2name,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2)

dotplot(enrich.17, showCategory = 20)


#==================================================
#Individual KEGG analysis - zt21
#==================================================
#run compare cluster
# Filter annotation to genes significant at each timepoint
annot.21 <- annot %>% filter(locusName %in% sig.21.genes)

# Pull KO terms, drop blanks, split any multi-KO entries
sig.21.ko <- annot.21 %>% filter(!is.na(KO) & KO != "") %>% pull(KO)
sig.21.ko <- unlist(strsplit(sig.21.ko, "[ ,]+"))

enrich.21 <- enricher(
  gene          = sig.21.ko,      # vector of KO terms
  TERM2GENE     = term2gene,
  TERM2NAME     = term2name,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2)

dotplot(enrich.21, showCategory = 20)


ck3 <- compareCluster(geneCluster = list(zt3 = sig.3.ko,
                                         zt7 = sig.7.ko,
                                         zt17 = sig.17.ko,
                                         zt21 = sig.21.ko),
                      fun = "enrichKEGG",
                      organism = "ko",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH")

as.data.frame(ck3)
dotplot(ck3, showCategory = 15)+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8))
















#===========================
#KEGG 71 candidates
#===========================
# Check all IDs match annotation table?
setdiff(narrowed_strict, annot$locusName)

# Filter annotation to gene list
annot.candidates<- annot %>% filter(locusName %in% narrowed_strict)

# Pull KO terms, drop blanks, split any multi-KO entries
candidate.ko <- annot.candidates %>% filter(!is.na(KO) & KO != "") %>% pull(KO)
candidate.ko <- unlist(strsplit(candidate.ko, "[ ,]+"))

# Run KEGG enrichment
enrich.candidates <- enricher(
  gene          = candidate.ko,
  TERM2GENE     = term2gene,
  TERM2NAME     = term2name,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

head(as.data.frame(enrich.candidates))

dotplot(enrich.candidates, showCategory = 20) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))

