# Load required libraries
library(WGCNA)
library(clusterProfiler)
library(org.At.tair.db)
# Preprocess data
data_for_wgcna <- (data_for_pca)  
data_for_wgcna <- data_for_wgcna[colSums(is.na(data_for_wgcna)) == 0]  # Remove genes with missing values

# Choose a soft-thresholding power
powers = c(1:10)
sft = pickSoftThreshold(data_for_wgcna, powerVector = powers, verbose = 5)

# Assuming soft-thresholding power is 6
power = 6
data_for_wgcna_matrix <- as.matrix(data_for_wgcna)
adjacency = adjacency(data_for_wgcna_matrix, power = power)
TOM = TOMsimilarity(adjacency)  # Transform the adjacency into Topological Overlap Matrix (TOM)
dissTOM = 1-TOM

# Perform gene clustering
geneTree = hclust(as.dist(dissTOM), method="average")

# Define modules
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

# Gene IDs
# gene_ids <- c("gene1", "gene2", "gene3", ...)
# Get gene identifiers from the column names
gene_ids <- b$ID


# Convert module assignments into a list of gene vectors
geneLists = split(gene_ids, dynamicMods)

# Run GO term enrichment for each module
enrichment_results <- lapply(geneLists, function(genes) {
  enrichGO(gene = genes,
           OrgDb = org.At.tair.db,
           keyType = "TAIR",
           ont = "BP",
           pAdjustMethod = "BH",
           qvalueCutoff = 0.05)
})

names(enrichmentResults) = paste0("Module", seq_along(enrichmentResults))

# Print the results
print(enrichmentResults)
