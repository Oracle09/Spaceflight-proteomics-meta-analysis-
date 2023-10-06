###7/6/2023
#AUTHOR: Gbolaga
# This is after MaxQuant have given the count numbers

library(limma)
library(edgeR)
library(DESeq2)

df<-read.csv("GLDS-38_soluble_proteome_WyattLab_soluble.csv")
# 1. Eliminate rows where columns 11-16 values are all 0
f <- df[!apply(df[,11:16] == 0, 1, all), ]
# 2. Calculate average for each unique value in column 2 for columns 11-16
library(dplyr)

f <- f %>%
  group_by(f[,2]) %>%
  summarise(across(11:16, mean, na.rm = TRUE))


#Change column names
colnames(f)[2:4]<-c("GC1", "GC2","GC3")
colnames(f)[5:7]<-c("F1", "F2","F3")

#conduct differential analysis to calculate p-value and log fold change
#Establish metadata
sample_info <- data.frame(Sample = c("GC1", "GC2","GC3","F1","F2","F3"),
                          condition = c("Ground","Ground","Ground","Flight", "Flight", "Flight"))
write.table(sample_info, file = "sample_info.txt", quote = FALSE, sep = "\t", row.names = FALSE)

sample_info$condition<-as.factor(sample_info$condition)
sample_info$condition<-relevel(sample_info$condition, ref="Ground")
### 
###set geneID as a row header
# Convert 'f' from a tibble to a data frame
f <- as.data.frame(f)
# Set the first column of 'f' as row names
row.names(f) <- f[,1]
# Remove the first column of 'f'
f <- f[,-1]

##########################################
#*Differential analysis using Limma-voom because proteomics data are not integers, they are continous
#*############################################

#Convert your data to a DGEList object for edgeR:
dge <- DGEList(counts = f)
dge <- estimateDisp(dge, sample_info)
vdata <- voom(dge, normalize.method="quantile")
fit <- lmFit(vdata, design=model.matrix(~condition, data=metadata))
fit <- eBayes(fit)
results <- topTable(fit, coef="conditionGround_control", n=Inf)

#make row header the first column

g<-rownames_to_column(f, "Row")
h<-rownames_to_column(results, "Row")

j<-merge(h,g, by ="Row")

write.csv(j, "BRIC20 mem protein.csv")



##########################################
#*Differential analysis using DESeq2, just out of curiosity
#*############################################
#convert protein intensity to integer
# Assuming dge is your DGEList object
dge$counts <- round(dge$counts) 

# Now create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = dge$counts,
                              colData = sample_info,
                              design = ~ condition)
# Run the DESeq differential expression pipeline
dds <- DESeq(dds)

# Get the results
res <- results(dds)

#adjusting for signifcance
res <- results(dds, alpha = 0.05)

#Volcano plot
res$col <- "black"  # default color
res$col[res$padj <= 0.05 & res$log2FoldChange > 0] <- "blue"
res$col[res$padj <= 0.05 & res$log2FoldChange < 0] <- "red"

# Now use this new column in the plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, col=col, main="Volcano plot", xlab="log2(FoldChange)", ylab="-log10(pvalue)"))

write.csv(res, file = "BRIC20 soluble.csv")
