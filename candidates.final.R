#Upset plot
library(UpSetR)

# Build a named list of all your gene sets
gene_sets <- list(
  Cluster5_fuzz = cluster5.genes.fuzz,
  Cluster5_core = cluster5.genes.core,
  Cluster4_HC = cluster4.genes.hc,
  Cluster5_HC = cluster5.genes.hc,
  Cluster6_HC = cluster6.genes.hc,
  Cluster8_HC = cluster8.genes.hc,
  rhythmic.upreg17 = comp1,
  rhythmic = sig_rhythmic
)

upset(fromList(gene_sets), 
      nsets = length(gene_sets),
      order.by = "freq",
      nintersects = 20)




#Summary table
library(dplyr)

summary_tbl <- tibble(
  set = names(gene_sets)[1:6],
  set_size = sapply(gene_sets[1:6], length),
  n_comp1 = sapply(gene_sets[1:6], function(x) length(intersect(x, comp1))),
  n_sig_rhythmic = sapply(gene_sets[1:6], function(x) length(intersect(x, sig_rhythmic))),
  n_both = sapply(gene_sets[1:6], function(x) length(Reduce(intersect, list(x, comp1, sig_rhythmic))))
) %>%
  mutate(
    pct_comp1 = round(100 * n_comp1 / set_size, 1),
    pct_sig_rhythmic = round(100 * n_sig_rhythmic / set_size, 1)
  )

summary_tbl


recurring_genes <- lapply(gene_sets[1:6], function(x) {
  Reduce(intersect, list(x, comp1, sig_rhythmic))
})

sapply(recurring_genes, length)   # quick counts





# Just pool all genes that show up in sig_rhythmic AND any of these six sets
# Genes in sig_rhythmic AND each individual set
rhythmic_overlap <- lapply(gene_sets[1:6], function(x) intersect(x, sig_rhythmic))

sapply(rhythmic_overlap, length)   # should match your n_sig_rhythmic column exactly

# Pool into one narrowed candidate list (unique genes across all six sets)
narrowed_candidates <- unique(unlist(rhythmic_overlap))
length(narrowed_candidates)


#Have to be in comp 1 as well
strict_candidates <- lapply(gene_sets[1:6], function(x) {
  Reduce(intersect, list(x, comp1, sig_rhythmic))
})

sapply(strict_candidates, length)
narrowed_strict <- unique(unlist(strict_candidates))
length(narrowed_strict)



























library(dplyr)
library(tibble)
library(openxlsx)

# --------------------------------------------------------
# 1. Start with narrowed_strict gene IDs
# --------------------------------------------------------
master_df <- tibble(Gene.Name = narrowed_strict)

# --------------------------------------------------------
# 2. Gene annotation
# --------------------------------------------------------
master_df <- master_df %>%
  left_join(genes, by = "Gene.Name")

# --------------------------------------------------------
# 3. Cluster5_fuzz membership score
#    (from your mfuzz object's $membership matrix)
# --------------------------------------------------------
# Adjust "c.fuzz" to whatever your mfuzz() output object is called
membership_matrix <- mfuzz.res$membership

fuzz_membership_df <- membership_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Gene.Name") %>%
  select(Gene.Name, cluster5_fuzz_membership = `5`)  # column "5" = cluster 5 score

master_df <- master_df %>%
  left_join(fuzz_membership_df, by = "Gene.Name")

# --------------------------------------------------------
# 4. Cluster5_core TRUE/FALSE
# --------------------------------------------------------
master_df <- master_df %>%
  mutate(cluster5_core = Gene.Name %in% cluster5.genes.core)

# --------------------------------------------------------
# 5. Which HC cluster each gene is in
# --------------------------------------------------------
# gene.clusters.hc should have columns: gene, cluster
hc_lookup <- gene.clusters.hc %>%
  select(Gene.Name = gene, hc_cluster = cluster)

master_df <- master_df %>%
  left_join(hc_lookup, by = "Gene.Name")

# --------------------------------------------------------
# 6. comp2a / comp2b / comp3 membership (TRUE/FALSE)
# --------------------------------------------------------
master_df <- master_df %>%
  mutate(
    'Interaction' = Gene.Name %in% comp2a,
    'DE'       = Gene.Name %in% comp2b,
    'Up at 17 not 7'       = Gene.Name %in% comp3
  )

# --------------------------------------------------------
# 7. logFC at ZT17 from res.17.df
# --------------------------------------------------------
# Adjust column names below to match res.17.df — check with colnames(res.17.df)
# Assuming it has a gene ID column and a log2FoldChange column
logfc_lookup <- res.17.df %>%
  rownames_to_column(var = "Gene.Name") %>%   # skip this line if gene ID is already a column, not rownames
  select(Gene.Name, logFC_ZT17 = log2FoldChange)  # adjust "log2FoldChange" to actual column name

master_df <- master_df %>%
  left_join(logfc_lookup, by = "Gene.Name")

# --------------------------------------------------------
# 8. Write to Excel
# --------------------------------------------------------
write.xlsx(master_df, "Book 3.xlsx")



















tiptop.candidates = c('Bradi1g34470','Bradi1g61450','Bradi1g63690','Bradi2g04810',
                   'Bradi2g05226','Bradi2g11640','Bradi2g52680','Bradi3g20740',
                   'Bradi3g27700','Bradi5g19610')
top.candidates = c('Bradi1g34470','Bradi1g61450','Bradi1g63690','Bradi2g04810',
                   'Bradi2g05226','Bradi2g11640','Bradi2g52680','Bradi3g20740',
                   'Bradi3g27700','Bradi5g19610','Bradi1g18720','Bradi1g52374',
                   'Bradi1g60670','Bradi2g04900','Bradi4g30240')


#=====================
#Plotting candidates (good expression profiles)
#=====================

#OS07G0661900 protein - Bradi1g18720 (Pentatricopeptide repeat-containing protein)
res.3.df[res.3.df$gene=='Bradi1g18720',] #logFC = 0.899
res.7.df[res.7.df$gene=='Bradi1g18720',] #logFC = 2.502
res.17.df[res.17.df$gene=='Bradi1g18720',] #logFC = 1.45
res.21.df[res.21.df$gene=='Bradi1g18720',] #logFC = 1.267

plot.gene("Bradi1g18720", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )

#MYB-LIKE HTH TRANSCRIPTIONAL REGULATOR FAMILY PROTEIN - Bradi1g34470
res.3.df[res.3.df$gene=='Bradi1g34470',] #logFC = -0.249
res.7.df[res.7.df$gene=='Bradi1g34470',] #logFC = 3.209
res.17.df[res.17.df$gene=='Bradi1g34470',] #logFC = 3.949
res.21.df[res.21.df$gene=='Bradi1g34470',] #logFC = 2.125

plot.gene("Bradi1g34470", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )


#OS10G0100500 PROTEIN - Bradi1g61450 (U-box domain containing protein)
res.3.df[res.3.df$gene=='Bradi1g61450',] #logFC =  0.487
res.7.df[res.7.df$gene=='Bradi1g61450',] #logFC = 1.816
res.17.df[res.17.df$gene=='Bradi1g61450',] #logFC = 1.345
res.21.df[res.21.df$gene=='Bradi1g61450',] #logFC = 0.487

plot.gene("Bradi1g61450", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )


#MYB FAMILY TRANSCRIPTION FACTOR APL - Bradi1g63690
res.3.df[res.3.df$gene=='Bradi1g63690',] #logFC = 0.887
res.7.df[res.7.df$gene=='Bradi1g63690',] #logFC = 1.474
res.17.df[res.17.df$gene=='Bradi1g63690',] #logFC = 1.944
res.21.df[res.21.df$gene=='Bradi1g63690',] #logFC = 0.924

plot.gene("Bradi1g63690", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )


#MYB FAMILY TRANSCRIPTION FACTOR EFM - Bradi2g04810
res.3.df[res.3.df$gene=='Bradi2g04810',] #logFC =  2.417
res.7.df[res.7.df$gene=='Bradi2g04810',] #logFC = 3.655
res.17.df[res.17.df$gene=='Bradi2g04810',] #logFC =  2.611
res.21.df[res.21.df$gene=='Bradi2g04810',] #logFC = 2.311

plot.gene("Bradi2g04810", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )


#GIGANTEA (GI) - Bradi2g05226
res.3.df[res.3.df$gene=='Bradi2g05226',] #logFC = 3.148
res.7.df[res.7.df$gene=='Bradi2g05226',] #logFC = 1.519
res.17.df[res.17.df$gene=='Bradi2g05226',] #logFC = 1.827
res.21.df[res.21.df$gene=='Bradi2g05226',] #logFC = 2.370

plot.gene("Bradi2g05226", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )



#SERINE_THREONINE-PROTEIN KINASE AFC2 - Bradi2g11640
res.3.df[res.3.df$gene=='Bradi2g11640',] #logFC = 0.209
res.7.df[res.7.df$gene=='Bradi2g11640',] #logFC = 0.826
res.17.df[res.17.df$gene=='Bradi2g11640',] #logFC = 1.449
res.21.df[res.21.df$gene=='Bradi2g11640',] #logFC = 1.066

plot.gene("Bradi2g11640", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")+
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )


#WAT1-RELATED PROTEIN - Bradi2g52680
res.3.df[res.3.df$gene=='Bradi2g52680',] #logFC = -1.986
res.7.df[res.7.df$gene=='Bradi2g52680',] #logFC = -0.155
res.17.df[res.17.df$gene=='Bradi2g52680',] #logFC = 6.591
res.21.df[res.21.df$gene=='Bradi2g52680',] #logFC = 3.932

plot.gene("Bradi2g52680", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )


#Bradi3g20740 (Ambiguous, maybe pentatricopeptide repeat containing protein)
res.3.df[res.3.df$gene=='Bradi3g20740',] #logFC = -0.130
res.7.df[res.7.df$gene=='Bradi3g20740',] #logFC = 0.598
res.17.df[res.17.df$gene=='Bradi3g20740',] #logFC = 1.102
res.21.df[res.21.df$gene=='Bradi3g20740',] #logFC = 0.962

plot.gene("Bradi3g20740", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )


#Predicted E3 ubiquitin ligase // Predicted E3 ubiquitin ligase - Bradi3g27700
res.3.df[res.3.df$gene=='Bradi3g27700',] #logFC =  0.970
res.7.df[res.7.df$gene=='Bradi3g27700',] #logFC = 0.735
res.17.df[res.17.df$gene=='Bradi3g27700',] #logFC = 1.317
res.21.df[res.21.df$gene=='Bradi3g27700',] #logFC = 0.512

plot.gene("Bradi3g27700", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )


#OS04G0586000 PROTEIN - Bradi5g19610 (FCS-Like zinc finger)
res.3.df[res.3.df$gene=='Bradi5g19610',] #logFC = 1.656
res.7.df[res.7.df$gene=='Bradi5g19610',] #logFC =  2.342
res.17.df[res.17.df$gene=='Bradi5g19610',] #logFC = 1.665
res.21.df[res.21.df$gene=='Bradi5g19610',] #logFC = 2.234

plot.gene("Bradi5g19610", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )
#=====================
#Plotting candidates (meh expression profiles)
#=====================
#OS07G0661900 PROTEIN - Bradi1g18720 (Pentatricopeptide repeat-containing protein)
res.3.df[res.3.df$gene=='Bradi1g18720',] #logFC = 0.899
res.7.df[res.7.df$gene=='Bradi1g18720',] #logFC = 2.502
res.17.df[res.17.df$gene=='Bradi1g18720',] #logFC = 1.452
res.21.df[res.21.df$gene=='Bradi1g18720',] #logFC = 1.267

plot.gene("Bradi1g18720", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")



#T6K12.7 PROTEIN - Bradi1g52374 (Acetyltransferase, U-box domain containing protein (acts with E3 ligase)
res.3.df[res.3.df$gene=='Bradi1g52374',] #logFC = 0.115
res.7.df[res.7.df$gene=='Bradi1g52374',] #logFC = 0.966
res.17.df[res.17.df$gene=='Bradi1g52374',] #logFC = 1.072
res.21.df[res.21.df$gene=='Bradi1g52374',] #logFC = 0.765

plot.gene("Bradi1g52374", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")


#TRANSCRIPTION FACTOR SAC51-RELATED - Bradi1g60670
res.3.df[res.3.df$gene=='Bradi1g60670',] #logFC = 0.774
res.7.df[res.7.df$gene=='Bradi1g60670',] #logFC =  2.160
res.17.df[res.17.df$gene=='Bradi1g60670',] #logFC =  1.262
res.21.df[res.21.df$gene=='Bradi1g60670',] #logFC =  1.262

plot.gene("Bradi1g60670", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")


#SERINE_THREONINE PROTEIN PHOSPHATASE 2A REGULATORY SUBUNIT - Bradi2g04900
res.3.df[res.3.df$gene=='Bradi2g04900',] #logFC =  2.182
res.7.df[res.7.df$gene=='Bradi2g04900',] #logFC = 3.058
res.17.df[res.17.df$gene=='Bradi2g04900',] #logFC = 2.294
res.21.df[res.21.df$gene=='Bradi2g04900',] #logFC = 1.791

plot.gene("Bradi2g04900", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")


#ZINC-FINGER HOMEODOMAIN PROTEIN 9 - Bradi4g30240
res.3.df[res.3.df$gene=='Bradi4g30240',] #logFC = 0.542
res.7.df[res.7.df$gene=='Bradi4g30240',] #logFC = 1.516
res.17.df[res.17.df$gene=='Bradi4g30240',] #logFC = 1.258
res.21.df[res.21.df$gene=='Bradi4g30240',] #logFC = 1.252

plot.gene("Bradi4g30240", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")
