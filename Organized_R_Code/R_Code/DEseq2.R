setwd("/Users/jeffersonchen/programming/YuceBio/YuceBio/Organized_R_Code/R_Code")
library("devtools")
library("magrittr") # needs to be run every time you start R and want to use %>%
library("dplyr")    # alternative
library("clusterProfiler")
library( "DESeq2" )
library("ggplot2")
library("dplyr")
library("htmltools")
library("limma")
library("edgeR")
# BiocManager::install("DESeq2")
# install_github("jmzeng1314/GEOmirror")

# 加载featureCount 结果
data <- read.delim("./../Data/GSE135222.hg38.genename.tpm.txt", sep = '\t')
rownames(data) <- data$X1
countData <- as.matrix(data[, -1])
rownames(countData) <- rownames(data)
# data <- round(data, digits = 3)
# 读取样本的基础信息
phenotype_info <- read.delim("./../Data/GSE135222_series_matrix.txt", sep = '\t')

sample_id <- read.delim("./../Data/Sample_Names.txt", sep = '\t', header = F)
colnames(sample_id) <- c("Samples", "Condition", "Dataset")
rownames(sample_id) <- sample_id$Samples

sample_id$Condition <- ifelse(sample_id$Condition == "Nonresponder", "NR", "R") %>% factor( levels = c("NR", "R")) %>% relevel(ref = "NR")


group_list <- sample_id$Condition
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(countData)
v <- voom(countData, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
fit
output <- topTable(fit, coef=2,n=Inf)
sum(output$adj.P.Val<0.05)
