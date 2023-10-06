##PROTEOMIC META-ANALYSIS 
#7/11/2023
#AUTHOR: GBOLAGA OLANREWAJU (Wyatt Lab, Ohio University)

#Disclaimer::: Not a pipeline, code fragments to achieve my aims
#future prospect to develop into pipeline
library(ggfortify)
library(ggbiplot)
library(tidyr)
library(ggplot2)
library(vegan)
library(readr)
library(cluster)
library(dendextend)
library(pheatmap)

BRIC20<-read.csv("BRIC20.csv")
abrs<-read.csv("abrs R.csv")
mazar<-read.csv("Mazar.csv")
BRICR<-read.csv("BRIC LED root.csv")
BRICS<-read.csv("BRIC LED shoot.csv")


###################################
#Change column 1 of all dataset to "ID"
##################################################

# Get all object names in the environment
all_obj_names <- ls()

# Loop over each object name
for (obj_name in all_obj_names) {
  
  # Check if the object is a data frame
  if (is.data.frame(get(obj_name))) {
    
    # Get the data frame
    df <- get(obj_name)
    
    # Change the name of the first column
    names(df)[1] <- "ID"
    
    # Assign the modified data frame back to the variable
    assign(obj_name, df)
  }
}


###################################
#Merge all dataset in environment by gene ID
##################################################

# Get all object names in the environment
all_obj_names <- ls()

# Create an empty list to store all the data frames
dfs <- list()

# Loop over each object name
for (obj_name in all_obj_names) {
  
  # Check if the object is a data frame
  if (is.data.frame(get(obj_name))) {
    
    # Get the data frame
    df <- get(obj_name)
    
    # Prefix the column names with the data frame name, except for the ID column
    names(df)[-1] <- paste0(obj_name, "_", names(df)[-1])
    
    # Add the modified data frame to the list
    dfs[[obj_name]] <- df
  }
}

# Use Reduce to iteratively merge all data frames in the list by "ID"
merged_df <- Reduce(function(x, y) merge(x, y, by="ID", all=TRUE), dfs)


##################################
#remove redundant column and retain only LFC
###############################################
merged_df <- merged_df[,-c(5,6,8,10,11,12,14)]

#remove all ROWS with no gene ID
a<-merged_df[merged_df$ID != "", ]

#convert column 2 to numeric from chr
a[, 2] <- as.numeric(a[, 2])

####################################
#Change column names to what I want
############################################
# Get the current column names
current_names <- colnames(a)

# Specify the new names
new_names <- c("ABRS leaf", "ABRS root", "BRIC 20", "BRIC LED root", 
               "BRIC LED shoot", "EMCS")

# Replace the old names with the new ones
current_names[2:7] <- new_names

# Assign the new column names back to the data frame
colnames(a) <- current_names


####################################
#Make PCA
####################################


# Get the data excluding the gene ID column
data_for_pca <- a[,-1]

# Ensure the data are numeric
data_for_pca <- data.frame(lapply(data_for_pca, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)))


####################################
#Transposed PCA
####################################


data_for_pca_t <- t(data_for_pca)



pca_result <- prcomp(data_for_pca_t, scale. = TRUE)

# Convert PCA result to a data frame
pca_df <- as.data.frame(pca_result$x)

# Add the group column to the PCA data frame
pca_df$group <- rownames(pca_df)

# Calculate percentage variance explained by each PC
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100



ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group), size = 4) +  # Map color to group here
  geom_text(aes(label = group), size = 5, vjust = 2, hjust = 0.4) +
  xlab(paste0("PC1: ", round(var_explained[1], 2), "% variance")) + 
  ylab(paste0("PC2: ", round(var_explained[2], 2), "% variance")) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 14)) 

# Adjust plot limits to avoid cutting off labels
p <- p + expand_limits(x = c(min(pca_df$PC1, na.rm = TRUE) - 1, max(pca_df$PC1, na.rm = TRUE) + 1),
                       y = c(min(pca_df$PC2, na.rm = TRUE) - 1, max(pca_df$PC2, na.rm = TRUE) + 1))




####################################
#Using PCoA instead of PCA
####################################


b<-a

exp_data<-b[,2:7]

# Replace NA values with 0
exp_data[is.na(exp_data)] <- 0

# Transpose data
transposed_data <- t(exp_data)

# Calculate dissimilarity matrix
diss <- vegdist(transposed_data, method = "horn")

# Perform PCoA
pcoa <- cmdscale(diss, eig = TRUE, k = 4) #k is the number of dimensions, it's usually 2 or 3
# Bind the samples back to the ordination

# Calculate variance explained
var_explained <- pcoa$eig/sum(pcoa$eig)


df <- as.data.frame(pcoa$points)
df$Samples <- rownames(df)

# Plot using ggplot2
ggplot(df, aes(x = V1, y = V2, color = Samples)) +
  geom_point(size = 4) +
  geom_text(aes(label = Samples), vjust = 1.5, hjust = 0.8, size = 5, check_overlap = TRUE) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold")  # Adjust size and face here
  ) +
  labs(
    x = paste("PC1: ", round(var_explained[1] * 100, 2), "% variance"),
    y = paste("PC2: ", round(var_explained[2] * 100, 2), "% variance"),
    color = "Experimental Conditions"
  )

####################################
#Make Scree plot to see the variables with highest contribution
#################################### 
# Get variance explained by each PCoA coordinate
var_explained <- pcoa$eig / sum(pcoa$eig)

# Create the scree plot
par(mar = c(5, 5, 4, 2) + 0.1, oma = c(0, 0, 2, 0) + 0.1)

# Create the plot
plot(var_explained, type = "b", xlab = "PCoA Coordinate",
     ylab = "Proportion of Variance Explained",
     main = "Scree Plot for PCoA",
     col = "red")


####################################
#Make Eucledean hierachical clustering
#################################### 
# Create a heatmap with dendrograms
pheatmap(data_for_pca,
         scale = "row", # scale data by row (consider removing this if data is already scaled)
         clustering_distance_rows = dist(data_for_pca, method = "euclidean"), # Euclidean distance for row clustering
         clustering_distance_cols = dist(t(data_for_pca), method = "euclidean"), # Euclidean distance for column clustering
         clustering_method = "complete", # linkage method for clustering
         show_rownames = F,
         show_colnames = T)






# ###############
#WCGNA
#########################

f<-a

#There are duplicated Gene ID 
duplicated_rows <- a[duplicated(a$ID) | duplicated(a$ID, fromLast = TRUE),]
# Keep only the first instance of each duplicated ID
f <- a[!duplicated(a$ID), ]
rownames(f) <- f[,1]  
f <- f[,-1]


# Load required package
library(WGCNA)

# Allow multi-threading
allowWGCNAThreads()

# Choose a power for which the scale-free topology fit index reaches a satisfactory value (say 0.9)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(f, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Choose a power based on the plot
power = 6

# Calculate # Calculate adjacency matrix
adjacency = adjacency(f, power = power)

# Turn adjacency into Topological Overlap Matrix (TOM)
TOM = TOMsimilarity(adjacency, TOMType = "unsigned", verbose = 5)

# Cluster the genes using hierarchical clustering
geneTree = hclust(as.dist(1-TOM), method = "average")
moduleColors = cutreeDynamic(dendro = geneTree, distM = 1-TOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)

# Plot the resulting clustering with colors representing modules
plotDendroAndColors(geneTree, moduleColors, "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

