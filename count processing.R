# ----------------------------------------------------------
# Provenance / AI assistance statement:
# ----------------------------------------------------------
# This script was developed by the author with assistance from
# Claude for code streamlining, naming clarity,
# and minor robustness improvements (e.g., type standardisation,
# safer subsetting, and improved annotation spacing). All analytical
# decisions (candidate gene list, contrasts, thresholds, and
# interpretation) were made by the author, and results were validated
# by the author.

setwd("//storage.hcs-p01.otago.ac.nz/biochemistry/Lab_Groups/brownfieldlab/Documents/Riley")

fcData = read.delim("jgi.counts.txt", comment.char="#")

fcData$Geneid <- gsub("\\.v3\\.2$", "", fcData$Geneid) #remove .v3.2 from gene IDs
colnames(fcData) = gsub(".*BAM\\.|\\.bam", "", colnames(fcData))
counts.unfiltered = fcData[,7:35]
gene_loci = fcData[,1:6]
rownames(counts.unfiltered) = fcData$Geneid

#change 3D5 and 3D6 to 3D1 and 3D3; 3L5 and 3L6 to 3L1 and 3L2 to keep replicate as only having 4 factors
colnames(counts.unfiltered)[colnames(counts.unfiltered) == "3D5"] = "3D1"
colnames(counts.unfiltered)[colnames(counts.unfiltered) == "3D6"] = "3D3"
colnames(counts.unfiltered)[colnames(counts.unfiltered) == "3L5"] = "3L1"
colnames(counts.unfiltered)[colnames(counts.unfiltered) == "3L6"] = "3L2"

#How many rows without count data?
sum(rowSums(counts.unfiltered) == 0)
#2675 rows 

#remove rows with no count data
counts=counts.unfiltered[rowSums(counts.unfiltered)>0,]

#remove 21D4 and 21L4 because multiQC showed they had really low reads
counts <- counts[, !colnames(counts) %in% c("21D4", "21L4")]

#change 3D5 and 3D6 to 3D1 and 3D3; 3L5 and 3L6 to 3L1 and 3L2 to keep replicate as only having 4 factors
colnames(counts.unfiltered)[colnames(counts.unfiltered) == "3D5"] = "3D1"
colnames(counts.unfiltered)[colnames(counts.unfiltered) == "3D6"] = "3D3"
colnames(counts.unfiltered)[colnames(counts.unfiltered) == "3L5"] = "3L1"
colnames(counts.unfiltered)[colnames(counts.unfiltered) == "3L6"] = "3L2"


#swap 7L3 and 7D2 (swapped in protocol)
match(c("7L3", "7D2"), names(counts))
colnames(counts)[c(26, 21)] <- colnames(counts)[c(21, 26)]

#swap 7D4 and 7L4 (swapped in protocol)
match(c("7D4", "7L4"), names(counts))
colnames(counts)[c(23, 27)] <- colnames(counts)[c(27, 23)]

#makes a list of the sample names
samples = colnames(counts)

## TO PUT TREATMENT AND TIMEPOINT AS METADATA FACTORS
# Split into timepoint + treatment + replicate
parts <- regmatches(samples, regexec("^([0-9]+)([A-Za-z]+)([0-9]+)$", samples))
parts <- do.call(rbind, parts)

# Build metadata frame
coldata = data.frame(
  sample = samples,
  treatment = factor(parts[,3]),
  timepoint = factor(parts[,2]),
  replicate = factor(parts[,4]),
  row.names = samples
)

# Reorder timepoints properly: NB+1 → NB+5 → NB+17, day 1 2 3
coldata$timepoint = factor(coldata$timepoint,
                           levels = c("3", "7", "17","21"),
                           labels = c("zt_3", "zt_7", "zt_17","zt_21"))

#to confirm it worked
levels(coldata$timepoint)

# Ensure counts is a numeric matrix
counts = as.matrix(counts)
storage.mode(counts) = "numeric"

# Library sizes (total reads per sample)
lib_sizes = colSums(counts)


#match Bradi ID to gene names to add to DESeq results table later
filtered.loci = rownames(counts)
writeLines(filtered.loci, "all.genes.txt") #send to Phytozome biomart to get gene names
genes=read.csv("mart.export.txt")
genes$Description <- sub("^[^-]*-\\s*", "", genes$Description) #remove PANTHER/PFAM ID

head(genes)

