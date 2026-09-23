#==========================
#GO enrichment analysis
#=============================
library(tidyr)
library(GO.db)
library(clusterProfiler)
library(dplyr)

#Build GO term to gene mapping from annotation file
goterm2gene <- annot.sub %>%
  filter(!is.na(GO) & GO != "") %>%
  dplyr::select(GO, locusName) %>%
  separate_rows(GO, sep = "[ ,;]+") %>%   # split multi-GO entries into one row each
  rename(term = GO, gene = locusName)

#Get GO term names for presentation
goterm2name <- AnnotationDbi::select(GO.db, keys = unique(goterm2gene$term),
                                      columns = "TERM", keytype = "GOID") %>%
  rename(term = GOID, name = TERM)

#Enrichment per timepoint
go.enrich.3 <- enricher(
  gene         = sig.3.genes,
  TERM2GENE    = goterm2gene,
  TERM2NAME    = goterm2name,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

as.data.frame(go.enrich.3)
dotplot(go.enrich.3, showCategory = 20)



go.enrich.7 <- enricher(
  gene         = sig.7.genes,
  TERM2GENE    = goterm2gene,
  TERM2NAME    = goterm2name,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

as.data.frame(go.enrich.7)
dotplot(go.enrich.7, showCategory = 20)


go.enrich.17 <- enricher(
  gene         = sig.17.genes,
  TERM2GENE    = goterm2gene,
  TERM2NAME    = goterm2name,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

as.data.frame(go.enrich.17)
dotplot(go.enrich.17, showCategory = 20)


go.enrich.21 <- enricher(
  gene         = sig.21.genes,
  TERM2GENE    = goterm2gene,
  TERM2NAME    = goterm2name,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

as.data.frame(go.enrich.21)
dotplot(go.enrich.21, showCategory = 20)

#Side by side of all 4 timepoints
ck.go <- compareCluster(
  geneCluster = list(
    zt3  = sig.3.genes,
    zt7  = sig.7.genes,
    zt17 = sig.17.genes,
    zt21 = sig.21.genes
  ),
  fun = "enricher",
  TERM2GENE = goterm2gene,
  TERM2NAME = goterm2name,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

dotplot(ck.go, showCategory = 20) + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8))

#Day timepoints compared to night timepoints
ck.go.daynight <- compareCluster(
  geneCluster = list(
    day   = DE.day,
    night = DE.night
  ),
  fun          = "enricher",
  TERM2GENE    = goterm2gene,
  TERM2NAME    = goterm2name,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

dotplot(ck.go.daynight, showCategory = 20) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8))

#7/17 DE comparison
go.enrich.7.17 <- enricher(
  gene         = sig.7.17,
  TERM2GENE    = goterm2gene,
  TERM2NAME    = goterm2name,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

as.data.frame(go.enrich.21)
dotplot(go.enrich.21, showCategory = 20)

#Number of genes with no GO terms
sum(DE.day %in% goterm2gene$gene)
#127
length(DE.day)
#203
sum(DE.night %in% goterm2gene$gene)
#144
length(DE.night)
#220

#=======================
#PCA loading GO
#======================
pca.genes.3  <- loadings.3$gene
pca.genes.7  <- loadings.7$gene
pca.genes.17 <- loadings.17$gene
pca.genes.21 <- loadings.21$gene

# --- define background (universe) ---
universe <- unique(annot$locusName)

# --- run enrichment for each timepoint ---
pca.3 <- enricher(
  gene          = pca.genes.3,
  universe      = universe,
  TERM2GENE     = goterm2gene,
  TERM2NAME     = goterm2name,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)


pca.7 <- enricher(
  gene          = pca.genes.7,
  universe      = universe,
  TERM2GENE     = goterm2gene,
  TERM2NAME     = goterm2name,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

pca.17 <- enricher(
  gene          = pca.genes.17,
  universe      = universe,
  TERM2GENE     = goterm2gene,
  TERM2NAME     = goterm2name,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

pca.21 <- enricher(
  gene          = pca.genes.21,
  universe      = universe,
  TERM2GENE     = goterm2gene,
  TERM2NAME     = goterm2name,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# inspect results
head(as.data.frame(pca.7), 20)
head(as.data.frame(pca.17), 20)

#visualize
library(enrichplot)
dotplot(pca.3,  showCategory = 15, title = "GO enrichment — ZT3")
dotplot(pca.7,  showCategory = 15, title = "GO enrichment — ZT7")
dotplot(pca.17, showCategory = 15, title = "GO enrichment — ZT17")
dotplot(pca.21, showCategory = 15, title = "GO enrichment — ZT21")



#=======================
#71 final candidates GO
#=======================
go.candidates <- enricher(
  gene         = narrowed_strict,
  TERM2GENE    = goterm2gene,
  TERM2NAME    = goterm2name,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

as.data.frame(go.candidates)
dotplot(go.candidates, showCategory = 20)
  
