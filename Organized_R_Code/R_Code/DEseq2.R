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
library(grid)
# BiocManager::install("DESeq2")
# install_github("jmzeng1314/GEOmirror")

# 加载featureCount 结果
# data <- read.delim("./../Data/GSE135222.hg38.genename.tpm.txt", sep = '\t')
# rownames(data) <- data$X1
# countData <- as.matrix(data[, -1])
# rownames(countData) <- rownames(data)
# # data <- round(data, digits = 3)
# # 读取样本的基础信息
# phenotype_info <- read.delim("./../Data/GSE135222_series_matrix.txt", sep = '\t')
# 
# sample_id <- read.delim("./../Data/Sample_Names.txt", sep = '\t', header = F)
# colnames(sample_id) <- c("Samples", "Condition", "Dataset")
# rownames(sample_id) <- sample_id$Samples
# 
# sample_id$Condition <- ifelse(sample_id$Condition == "Nonresponder", "NR", "R") %>% factor( levels = c("NR", "R")) %>% relevel(ref = "NR")
# 
# 
# 



# 加载featureCount 结果
# 20
#  [1] "ENSG00000224689.9"  "ENSG00000248485.2"  "ENSG00000016602.9"  "ENSG00000125740.14"
# [5] "ENSG00000189221.10" "ENSG00000121552.4"  "ENSG00000162772.17" "ENSG00000215158.9" 
# [9] "ENSG00000183347.15" "ENSG00000131203.13" "ENSG00000224163.4"  "ENSG00000149968.12"
# [13] "ENSG00000121904.18" "ENSG00000106701.13" "ENSG00000115414.21" "ENSG00000112837.17"
# [17] "ENSG00000115290.10" "ENSG00000060718.22" "ENSG00000100285.10" "ENSG00000100097.12"
with_p_exon <- read.delim("/Users/jeffersonchen/programming/YuceBio/Feature_Count/With_P_Exon.txt", sep = '\t')

# 20
# [1] "ENSG00000224689.9"  "ENSG00000248485.2"  "ENSG00000121552.4"  "ENSG00000189221.10"
# [5] "ENSG00000016602.9"  "ENSG00000162772.17" "ENSG00000183347.15" "ENSG00000186197.15"
# [9] "ENSG00000215158.9"  "ENSG00000116141.17" "ENSG00000224163.4"  "ENSG00000149968.12"
# [13] "ENSG00000121904.18" "ENSG00000115414.21" "ENSG00000106701.13" "ENSG00000112837.17"
# [17] "ENSG00000060718.22" "ENSG00000115290.10" "ENSG00000100285.10" "ENSG00000100097.12"
with_p_gene <- read.delim("/Users/jeffersonchen/programming/YuceBio/Feature_Count/With_P_Gene.txt", sep = '\t')

data <- with_p_gene

# 重新命名数据列名
sample_names <- c("SRR1382274", "SRR1382275", "SRR1382276", "SRR1382277", "SRR1382278", "SRR1382279", "SRR1382280", "SRR1382281", "SRR1382282")
names(data)[7:15] <- sample_names
data <- data[, c(1, 7:15)]
data <- aggregate(data, .~ Geneid, FUN = "sum")
rownames(data) <- data$Geneid
data <- data[, -1]
# 取出所有count值
countData <- as.matrix(data)

# write.table(countData, "/Users/jeffersonchen/programming/YuceBio/YuceBio/Day6/Count_Data.txt", sep = '\t', row.names = TRUE)

# 略微调整数据
rownames(countData) <- rownames(data)

# 读取样本的基础信息
database <- read.delim("../Data/Original_Version_Info.txt", sep = " ")
colnames(database) <- c("SAMPLE_NAME", "SAMPLE_TYPE")

type <- c(factor(database$SAMPLE_TYPE))
type <- relevel(type, ref = "ND")
# 整理数据
database <- data.frame(name=sample_names, condition=type)
write.table(database, file = "../Data/Original_Version_Info.txt")
# ***********************
# implement DESEQ2


dds <- DESeqDataSetFromMatrix(countData, colData=database, design = ~ condition)

# 取行合大于一的数据进行筛选
dds <- dds[ rowSums(counts(dds)) > 1, ]
# dds$condition <- factor(dds$condition, levels = c("D","ND"))
# dds$condition <- relevel(dds$condition, ref = "ND")
countData <- countData[rowSums(countData) > 1,]
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds)

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)


resdata <- resdata[!is.na(resdata$padj), ]
top_genes <- resdata[(resdata$padj <0.05) & (abs(resdata$log2FoldChange) > 1), ]
fpkm_data <- read.delim("../Results/With_P_Gene/With_P_Genefpkm.txt", sep = '\t')
colnames(fpkm_data)[2:10] <- sample_names
fpkm_data <- aggregate(fpkm_data, .~ Geneid, FUN = "sum")
rownames(fpkm_data) <- fpkm_data$Geneid

fpkm_data <-fpkm_data[intersect(fpkm_data$Geneid, top_genes$Row.names), ]
top_genes<- data.frame(fpkm_data[, -1], top_genes[, c(1:7)])

top_genes <- top_genes[order(top_genes$log2FoldChange, decreasing = TRUE),]
rld <- rlog(dds, blind=FALSE)
assay(rld)
mat <- assay(rld)[top_genes$Row.names, sample_names]
ncol(top_genes)
colnames(mat) <- sample_names
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale))
colnames(mat.scaled) <- colnames(mat)

num_keep <- nrow(mat.scaled)
row_keep <- c(1:nrow(mat.scaled))

l2_val <- as.matrix(top_genes[row_keep,]$log2FoldChange)
colnames(l2_val) <- "logFC"

mean <- as.matrix(top_genes[row_keep,]$baseMean)
colnames(mean) <- "AveExpr"


# BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library(dendextend)

row_dend <- as.dendrogram(hclust(dist(mat.scaled)))
row_dend <- color_branches(row_dend, k = 2)


col_logFC <- colorRamp2(c(min(l2_val), 0, max(l2_val)), c("blue", "white", "red"))

col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), height = unit(2, "cm")))

rl = top_genes$Row.names[row_keep]
h1 <- Heatmap(mat.scaled[row_keep,], 
              column_labels = sample_names, 
              name = "Z-score", 
              cluster_columns = TRUE, 
              row_names_gp = gpar(fontsize = 3), 
              row_labels = rl, 
              cluster_rows = row_dend, 
              column_km = 2, 
              column_names_centered = T,
              column_title = c("ND", "D"),
              column_title_gp = gpar(fill = c("purple", "orange"), 
                                     label = textGrob(c("ND", "D"))))



# row_labels = top_genes$Row.names[row_keep]
# cell_fun = function(j, i, x, y, w, h, col){grid.text(round(l2_val[i, j], 2), x, y)})
# h2 <- Heatmap(l2_val, cluster_rows = F, name = "logFC", top_annotation = ha, col = col_logFC)

 # row_labels = top_genes$Row.names[row_keep]
# cell_fun = function(j, i, x, y, w, h, col){grid.text(round(l2_val[i, j], 2), x, y)})
# h3 <- Heatmap(mean, cluster_rows = F, name = "AveExpr", col = col_AveExpr)

# h <- h1 + h2 +h3
pdf("HeatMap.pdf", width = 10, height = 12)
h1
dev.off()

write.table(res, "./../Results/With_P_Gene/res_output.txt", sep = '\t')

newdata <- read.delim("./../Results/With_P_Gene/res_output.txt", sep = '\t')

p1 <- ggplot(newdata, aes(log2FoldChange, -log(padj,10))) + # -log10 conversion  s
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1


# 进行分类，根据p值和log值将数据分为三等
new<- newdata %>% 
  mutate(
    Expression = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )


new$Genes <- rownames(new)
p2 <- ggplot(new, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2



# 列出前10最具差异性的基因
library(knitr)
top <- 10
p3 <- ggplot(new, aes(log2FoldChange, -log(padj, 10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))

top_genes <- bind_rows(
  new %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top),
  new %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top)
)
library(ggrepel)
p3 <-  p3 +
  geom_label_repel(data = top_genes,
                   mapping = aes(log2FoldChange, -log(padj,10), label = Genes),
                   size = 2)

pdf("Volcano.pdf")
p3
dev.off()




