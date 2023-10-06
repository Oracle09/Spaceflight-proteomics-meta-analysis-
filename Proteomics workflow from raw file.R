BiocManager::install("mzR")
BiocManager::install("S4Vectors")
BiocManager::install("BiocParallel")
BiocManager::install("RforProteomics")
BiocManager::install("Peptides")
BiocManager::install("Biostrings")
BiocManager::install("MSnbase")
BiocManager::install(c("MSnID", "limma"))
library(mzR)
library(S4Vectors)
library(BiocParallel)
library(RforProteomics)
library(Peptides)
library(Biostrings)
library(MSnbase)
library(MSnID)
library(limma)
library(tidyverse)
library(edgeR)
library(stringr)

#Input file will be from MaxQuant
protein_data <- read_delim("proteinGroups.txt")
# A brief look at the data
head(protein_data)

# Cleaning the data
# Keeping only columns required for the analysis - the columns with LFQ intensities
# Please replace "LFQ intensity" with the appropriate column name
# Adjust the pattern in str_detect to match   actual column names
data_lfq <- protein_data %>%
  select(starts_with("Protein"), str_detect(names(protein_data), "LFQ intensity"))

# Rename columns if necessary
# In this example, assume column names follow the format LFQ intensity_<Group>_<ReplicateNumber>
names(data_lfq) <- str_replace(names(data_lfq), "LFQ intensity_", "")

# A brief look at the cleaned data
head(data_lfq)

# Transposing the data so that samples are in rows and proteins are in columns
data_lfq_t <- t(data_lfq)

# Create a group vector for the samples
# The group order should correspond to the order of samples in the data_lfq_t matrix
groups <- factor(c(rep("earth", 3), rep("spaceflight", 3)))

# Create a DGEList object for edgeR analysis
y <- DGEList(counts = data_lfq_t, group = groups)

# Filter to remove rows with too few observations
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalize data
y <- calcNormFactors(y)

# Design matrix
design <- model.matrix(~groups)

# Estimate dispersion
y <- estimateDisp(y, design)

# Differential expression analysis
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

# Get the results
results <- topTags(qlf, n = Inf)