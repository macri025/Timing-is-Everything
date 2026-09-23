library(ggplot2)

#===========================================================
#make objects to use for plotting
#===========================================================
counts.norm3 = counts(dds.3, normalized = TRUE)
counts.norm7 = counts(dds.7, normalized = TRUE)
counts.norm17 = counts(dds.17, normalized = TRUE)
counts.norm21 = counts(dds.21, normalized = TRUE)

norm.counts.list = list(counts.norm3, counts.norm7, counts.norm17, counts.norm21)
res.list = list(res.3.df, res.7.df, res.17.df, res.21.df)
dds.list <- list(zt3 = dds.3, zt7 = dds.7, zt17 = dds.17, zt21 = dds.21)
timepoint.labels <- c("ZT3", "ZT7", "ZT17", "ZT21")


#======================================
#Define plot function
#======================================
plot.gene = function(gene.id, norm.counts.list, dds.list, timepoint.labels, res.list,
                     annotated.table = NULL, id.col = "gene", desc.col = "Description") {
  all.df <- do.call(rbind, lapply(seq_along(norm.counts.list), function(i) {
    nc <- norm.counts.list[[i]]
    if (!(gene.id %in% rownames(nc))) return(NULL)
    data.frame(
      sample = colnames(nc),
      expression = nc[gene.id, ],
      treatment = colData(dds.list[[i]])$treatment,
      timepoint = timepoint.labels[i]
    )
  }))
  
  all.df$timepoint <- factor(all.df$timepoint, levels = timepoint.labels)
  
  # pull padj for this gene at each timepoint, convert to stars
  padj.vals <- sapply(res.list, function(res) {
    val <- res$padj[res$gene == gene.id]
    if (length(val) == 0) return(NA)
    val
  })
  
  lfc.vals <- sapply(res.list, function(res) {
    val <- res$log2FoldChange[res$gene == gene.id]
    if (length(val) == 0) return(NA)
    val
  })
  
  basemean.vals <- sapply(res.list, function(res) {
    val <- res$baseMean[res$gene == gene.id]
    if (length(val) == 0) return(NA)
    val
  })
  
  sig.stars <- mapply(function(p, lfc, bm) {
    if (is.na(p) || is.na(lfc) || is.na(bm)) return("")
    if (abs(lfc) <= 1) return("")
    if (bm <= 10) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    return("")
  }, padj.vals, lfc.vals, basemean.vals)
  
  # build a small data frame for star annotation, positioned above the tallest bar per timepoint
  star.df <- data.frame(
    timepoint = factor(timepoint.labels, levels = timepoint.labels),
    label = sig.stars,
    y = sapply(seq_along(timepoint.labels), function(i) {
      max(all.df$expression[all.df$timepoint == timepoint.labels[i]], na.rm = TRUE) * 1.1
    })
  )
  
  if (!is.null(annotated.table)) {
    gene.label <- annotated.table[[desc.col]][annotated.table[[id.col]] == gene.id][1]
    if (is.na(gene.label) || length(gene.label) == 0) gene.label <- gene.id
  } else {
    gene.label <- gene.id
  }
  ggplot(all.df, aes(x = timepoint, y = expression, fill = treatment)) +
    stat_summary(fun = mean, geom = "bar", position = "dodge", alpha = 0.7) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 position = position_dodge(width = 0.9), width = 0.2) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
               size = 1.5, alpha = 0.6, color = "black") +
    geom_text(data = star.df, aes(x = timepoint, y = y, label = label),
              inherit.aes = FALSE, size = 6, vjust = 0) +
    labs(title = gene.label, x = "Timepoint", y = "Normalized counts") +
    scale_fill_manual(values = c("D" = "#615E5E", "L" = "#FFD60A")) +
    theme_minimal()
}


#=================================================================
#Individual gene plotting
#=================================================================

#PPD1/PRR37 - Bradi1g16490
res.3.df[res.3.df$gene=='Bradi1g16490',] #logFC = 2.60
res.7.df[res.7.df$gene=='Bradi1g16490',] #logFC = 1.74
res.17.df[res.17.df$gene=='Bradi1g16490',] #logFC = 1.52
res.21.df[res.21.df$gene=='Bradi1g16490',] #logFC = 2.88

plot.gene("Bradi1g52374", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")

plot.gene("Bradi1g16490", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))

Filter(isTRUE, sapply(gene.lists, function(x) "Bradi1g16490" %in% x))
Filter(isTRUE, sapply(gene.lists.zt3, function(x) "Bradi1g16490" %in% x))


#ELF3 - Bradi2g14290
res.3.df[res.3.df$gene=='Bradi2g14290',] #not sig
res.7.df[res.7.df$gene=='Bradi2g14290',] #not sig
res.17.df[res.17.df$gene=='Bradi2g14290',] #not sig
res.21.df[res.21.df$gene=='Bradi2g14290',] #not sig

plot.gene("Bradi2g14290", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")+
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))

Filter(isTRUE, sapply(gene.lists, function(x) "Bradi2g14290" %in% x))
Filter(isTRUE, sapply(gene.lists.zt3, function(x) "Bradi2g14290" %in% x))

#HOMEODOMAIN-LIKE SUPERFAMILY PROTEIN-RELATED - LUX
res.3.df[res.3.df$gene=='Bradi2g62067',] #not sig
res.7.df[res.7.df$gene=='Bradi2g62067',] #LogFC = -3.55
res.17.df[res.17.df$gene=='Bradi2g62067',] #LogFC = -1.11
res.21.df[res.21.df$gene=='Bradi2g62067',] #not sig
plot.gene("Bradi2g62067", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))

Filter(isTRUE, sapply(gene.lists, function(x) "Bradi2g62067" %in% x))
Filter(isTRUE, sapply(gene.lists.zt3, function(x) "Bradi2g62067" %in% x))

#PROTEIN ELF4-LIKE 4
res.3.df[res.3.df$gene=='Bradi4g13227',] #not sig
res.7.df[res.7.df$gene=='Bradi4g13227',] #not sig
res.17.df[res.17.df$gene=='Bradi4g13227',] #not sig
res.21.df[res.21.df$gene=='Bradi4g13227',] #not sig

plot.gene("Bradi4g13227", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")+
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))


Filter(isTRUE, sapply(gene.lists, function(x) "Bradi4g13227" %in% x))
Filter(isTRUE, sapply(gene.lists.zt3, function(x) "Bradi4g13227" %in% x))


#phytochrome B (PHYB)
res.3.df[res.3.df$gene=='Bradi1g64360',] #not sig
res.7.df[res.7.df$gene=='Bradi1g64360',] #-1.66
res.17.df[res.17.df$gene=='Bradi1g64360',] #-1.00
res.21.df[res.21.df$gene=='Bradi1g64360',] #LogFC = -2.11
plot.gene("Bradi1g64360", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))


Filter(isTRUE, sapply(gene.lists, function(x) "Bradi1g64360" %in% x))
Filter(isTRUE, sapply(gene.lists.zt3, function(x) "Bradi1g64360" %in% x))

#PHYTOCHROME C // PHYTOCHROME 1
res.3.df[res.3.df$gene=='Bradi1g08400',] #not sig
res.7.df[res.7.df$gene=='Bradi1g08400',] #-1.66
res.17.df[res.17.df$gene=='Bradi1g08400',] #-1.00
res.21.df[res.21.df$gene=='Bradi1g08400',] #LogFC = -2.11
plot.gene("Bradi1g08400", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")+
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))

Filter(isTRUE, sapply(gene.lists, function(x) "Bradi1g08400" %in% x))
Filter(isTRUE, sapply(gene.lists.zt3, function(x) "Bradi1g08400" %in% x))



#GIGANTEA (GI) - Bradi2g05226
res.3.df[res.3.df$gene=='Bradi2g05226',] #LogFC = 3.15
res.7.df[res.7.df$gene=='Bradi2g05226',] #LogFC = 1.52
res.17.df[res.17.df$gene=='Bradi2g05226',] #LogFC = 1.83
res.21.df[res.21.df$gene=='Bradi2g05226',] #LogFC = 2.02

plot.gene("Bradi2g05226", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
            theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 20),
                  axis.text.x  = element_text(size = 18),
                  axis.text.y  = element_text(size = 18))

Filter(isTRUE, sapply(gene.lists, function(x) "Bradi2g05226" %in% x))
Filter(isTRUE, sapply(gene.lists.zt3, function(x) "Bradi2g05226" %in% x))


#protein FLOWERING LOCUS T (FT) - Bradi4g39730
res.3.df[res.3.df$gene=='Bradi4g39730',] #not sig
res.7.df[res.7.df$gene=='Bradi4g39730',] #not sig
res.17.df[res.17.df$gene=='Bradi4g39730',] #not sig
res.21.df[res.21.df$gene=='Bradi4g39730',] #not sig - no expression

plot.gene("Bradi4g39730", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))

ply(gene.lists, function(x) "Bradi4g39730" %in% x)
sapply(gene.lists.zt3, function(x) "Bradi4g39730" %in% x)

#FCS-LIKE ZINC FINGER 15 - Bradi3g47530
res.3.df[res.3.df$gene=='Bradi3g47530',] #LogFC = 1.02
res.7.df[res.7.df$gene=='Bradi3g47530',] #LogFC = 1.13
res.17.df[res.17.df$gene=='Bradi3g47530',] #not sig
res.21.df[res.21.df$gene=='Bradi3g47530',] #not sig

plot.gene("Bradi3g47530", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))


sapply(gene.lists, function(x) "Bradi3g47530" %in% x)
sapply(gene.lists.zt3, function(x) "Bradi3g47530" %in% x)

#FCS-LIKE ZINC FINGER 8 - Bradi1g12960
res.3.df[res.3.df$gene=='Bradi1g12960',] #not sig
res.7.df[res.7.df$gene=='Bradi1g12960',] #LogFC = 2.00
res.17.df[res.17.df$gene=='Bradi1g12960',] #not sig
res.21.df[res.21.df$gene=='Bradi1g12960',] #not sig

plot.gene("Bradi1g12960", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")+
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))


sapply(gene.lists, function(x) "Bradi1g12960" %in% x)
sapply(gene.lists.zt3, function(x) "Bradi1g12960" %in% x)

#CONSTANS-LIKE PROTEIN DAYS TO HEADING ON CHROMOSOME 2 - Bradi3g56490
res.3.df[res.3.df$gene=='Bradi3g56490',] #not sig
res.7.df[res.7.df$gene=='Bradi3g56490',] #not sig
res.17.df[res.17.df$gene=='Bradi3g56490',] #LogFC = -1.17
res.21.df[res.21.df$gene=='Bradi3g56490',] #LogFC = -1.24

plot.gene("Bradi3g56490", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        plot.margin  = margin(t = 20, r = 10, b = 10, l = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))


sapply(gene.lists, function(x) "Bradi3g56490" %in% x)
sapply(gene.lists.zt3, function(x) "Bradi3g56490" %in% x)

#ZINC FINGER PROTEIN CONSTANS-LIKE 5 - Bradi3g05800
res.3.df[res.3.df$gene=='Bradi3g05800',] #not sig
res.7.df[res.7.df$gene=='Bradi3g05800',] #not sig
res.17.df[res.17.df$gene=='Bradi3g05800',] #not sig
res.21.df[res.21.df$gene=='Bradi3g05800',] #LogFC = -1.28

plot.gene("Bradi3g05800", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")  +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        plot.margin  = margin(t = 20, r = 10, b = 10, l = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))

sapply(gene.lists, function(x) "Bradi3g05800" %in% x)
sapply(gene.lists.zt3, function(x) "Bradi3g05800" %in% x)

#	MYB-related transcription factor LHY (LHY) - morning complex Bradi3g48880
res.3.df[res.3.df$gene=='Bradi3g16515',] #LogFC = -1,24
res.7.df[res.7.df$gene=='Bradi3g16515',] #LogFC = -1.64
res.17.df[res.17.df$gene=='Bradi3g16515',] #LogFC = -1.60
res.21.df[res.21.df$gene=='Bradi3g16515',] #LogFC = -2.12
plot.gene("Bradi3g16515", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")+
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))
#suggests that 3 and 17 are super different, shows that clock working in the dark

#	pseudo-response regulator 1 (TOC1, APRR1) - clock gene, night time
res.3.df[res.3.df$gene=='Bradi3g48880',] #not sig
res.7.df[res.7.df$gene=='Bradi3g48880',] #LogFC= -1.06
res.17.df[res.17.df$gene=='Bradi3g48880',] #not sig
res.21.df[res.21.df$gene=='Bradi3g48880',] #not sig
plot.gene("Bradi3g48880", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))


#VRN1
res.3.df[res.3.df$gene=='Bradi1g08340',] #not sig
res.7.df[res.7.df$gene=='Bradi1g08340',] #LogFC = 1.04
res.17.df[res.17.df$gene=='Bradi1g08340',] #LogFC = 1.39
res.21.df[res.21.df$gene=='Bradi1g08340',] #LogFC = 2.24
plot.gene("Bradi1g08340", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))


#VRN2
res.3.df[res.3.df$gene=='Bradi3g10010 ',] 
res.7.df[res.7.df$gene=='Bradi3g10010 ',] 
res.17.df[res.17.df$gene=='Bradi3g10010 ',]
res.21.df[res.21.df$gene=='Bradi3g10010 ',]
plot.gene("Bradi3g10010 ", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))



#VRN3
res.3.df[res.3.df$gene=='Bradi1g48830 ',] 
res.7.df[res.7.df$gene=='Bradi1g48830 ',] 
res.17.df[res.17.df$gene=='Bradi1g48830 ',] 
res.21.df[res.21.df$gene=='Bradi1g48830 ',]
plot.gene("Bradi1g48830 ", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))

#FTL9: Bradi2g19670 - checked with BLAST homologous to Taylor 2025
res.3.df[res.3.df$gene=='Bradi2g19670',] #not sig
res.7.df[res.7.df$gene=='Bradi2g19670',] #LogFC = 1.04
res.17.df[res.17.df$gene=='Bradi2g19670',] #LogFC = 1.39
res.21.df[res.21.df$gene=='Bradi2g19670',] #LogFC = 2.24
plot.gene("Bradi2g19670", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))




#Prr73
res.3.df[res.3.df$gene=='Bradi1g65910',] #not sig
res.7.df[res.7.df$gene=='Bradi1g65910',] #not sig
res.17.df[res.17.df$gene=='Bradi1g65910',] #not sig
res.21.df[res.21.df$gene=='Bradi1g65910',] #not sig
plot.gene("Bradi1g65910", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")  +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))

#Prr95
res.3.df[res.3.df$gene=='Bradi4g36077',] #not sig
res.7.df[res.7.df$gene=='Bradi4g36077',] #not sig
res.17.df[res.17.df$gene=='Bradi4g36077',] #not sig
res.21.df[res.21.df$gene=='Bradi4g36077',] #not sig
plot.gene("Bradi4g36077", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18))


#CO10
plot.gene("Bradi1g11310", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description")

#=======================
#LpCO13 potential homologs
#========================
plot.gene("Bradi1g62420", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        plot.margin  = margin(t = 20, r = 10, b = 10, l = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))


plot.gene("Bradi1g187407", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.17.annotated, id.col = "gene", desc.col = "Description") #Not in gene list


plot.gene("Bradi1g52360", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        plot.margin  = margin(t = 20, r = 10, b = 10, l = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))

plot.gene("Bradi1g18407", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        plot.margin  = margin(t = 20, r = 10, b = 10, l = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))

plot.gene("Bradi3g41500", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        plot.margin  = margin(t = 20, r = 10, b = 10, l = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))

plot.gene("Bradi3g05800", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        plot.margin  = margin(t = 20, r = 10, b = 10, l = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))

plot.gene("Bradi3g48447", norm.counts.list, dds.list, timepoint.labels, res.list,
          res.3.annotated, id.col = "gene", desc.col = "Description") +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        plot.margin  = margin(t = 20, r = 10, b = 10, l = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))
