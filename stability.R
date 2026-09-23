# ============================================================
# Seed stability check for mfuzz
# ============================================================
library(Mfuzz)

n.runs <- 20
seeds <- 1:n.runs
c.val <- 5  # or loop over a few candidate c values, see below

# Store the hard cluster assignment (highest membership) for each run
membership.list <- vector("list", n.runs)

for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  res <- mfuzz(eset.std, c = c.val, m = m.est)
  membership.list[[i]] <- res$cluster  # named vector: feature -> cluster label
}

# Build a feature x run matrix of cluster assignments
feat.ids <- names(membership.list[[1]])
assign.mat <- sapply(membership.list, function(x) x[feat.ids])
rownames(assign.mat) <- feat.ids
colnames(assign.mat) <- paste0("seed_", seeds)



# Since cluster numbers are arbitrary per run, align them first
# by matching centroids across runs (highest correlation)
align.clusters <- function(ref.centers, target.centers) {
  cor.mat <- cor(t(ref.centers), t(target.centers))
  match.idx <- apply(cor.mat, 1, which.max)  # for each ref cluster, best match in target
  match.idx
}

centers.list <- vector("list", n.runs)
for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  res <- mfuzz(eset.std, c = c.val, m = m.est)
  centers.list[[i]] <- res$centers
}

ref.centers <- centers.list[[1]]
aligned.mat <- matrix(NA, nrow = length(feat.ids), ncol = n.runs,
                      dimnames = list(feat.ids, paste0("seed_", seeds)))
aligned.mat[,1] <- membership.list[[1]][feat.ids]

for (i in 2:n.runs) {
  map <- align.clusters(ref.centers, centers.list[[i]])
  # map: names = ref cluster idx, values = matching target cluster idx
  relabel <- setNames(as.integer(names(map)), map)  # invert: target label -> ref label
  aligned.mat[,i] <- relabel[as.character(membership.list[[i]][feat.ids])]
}

# ------------------------------------------------------------
# Stability metric: for each feature, % of runs matching the mode (most common) cluster
# ------------------------------------------------------------
stability <- apply(aligned.mat, 1, function(row) {
  tab <- table(row)
  max(tab) / length(row)
})

summary(stability)
hist(stability, breaks = 20, main = "Per-feature cluster stability across seeds",
     xlab = "Proportion of runs agreeing with modal cluster")

# How many features are "stable" (e.g. >=90% agreement) vs "unstable"
table(stability >= 0.9)
mean(stability >= 0.9)  # fraction of stable features

# Look at which features/clusters are least stable
unstable.feats <- names(stability)[stability < 0.9]
length(unstable.feats)

# ------------------------------------------------------------
# Repeat across a few candidate c values to compare overall stability
# ------------------------------------------------------------
check.stability.for.c <- function(c.val, n.runs = 20) {
  membership.list <- vector("list", n.runs)
  centers.list <- vector("list", n.runs)
  for (i in 1:n.runs) {
    set.seed(i)
    res <- mfuzz(eset.std, c = c.val, m = m.est)
    membership.list[[i]] <- res$cluster
    centers.list[[i]] <- res$centers
  }
  feat.ids <- names(membership.list[[1]])
  ref.centers <- centers.list[[1]]
  aligned.mat <- matrix(NA, nrow = length(feat.ids), ncol = n.runs,
                        dimnames = list(feat.ids, NULL))
  aligned.mat[,1] <- membership.list[[1]][feat.ids]
  for (i in 2:n.runs) {
    map <- align.clusters(ref.centers, centers.list[[i]])
    relabel <- setNames(as.integer(names(map)), map)
    aligned.mat[,i] <- relabel[as.character(membership.list[[i]][feat.ids])]
  }
  stability <- apply(aligned.mat, 1, function(row) max(table(row)) / length(row))
  mean(stability >= 0.9)
}

c.values <- 4:10
stability.by.c <- sapply(c.values, check.stability.for.c)
names(stability.by.c) <- c.values
barplot(stability.by.c, xlab = "c (number of clusters)",
        ylab = "Fraction of features stable (>=90% agreement)",
        main = "Cluster stability vs c")