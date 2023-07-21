setwd("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day4/R_Code/DeSeq2")
library("devtools")
library("magrittr") # needs to be run every time you start R and want to use %>%
library("dplyr")    # alternative
library("clusterProfiler")
library( "DESeq2" )
library("ggplot2")
library("dplyr")
library("htmltools")
# BiocManager::install("DESeq2")

# 加载featureCount 结果
data <- read.delim("/Users/jeffersonchen/programming/YuceBio/Feature_Count/new.txt", sep = '\t')


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
database <- read.delim("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day4/Summarized_Data/SraRunTable9_GSE168009.txt", sep = ",")

# 重新命名数据行名
rownames(database) <- database$Run
# 将nodurable benefit和durable benefit缩写成ND和D
database$SAMPLE_TYPE <- factor(ifelse(database$SAMPLE_TYPE == "no durable benefit", "ND", "D"))
type <- c(factor(database$SAMPLE_TYPE))
type <- relevel(type, ref = "ND")
# 整理数据
database <- data.frame(name=sample_names, condition=type)
write.table(database, file = "./Original_Version_Info.txt")
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

# 储存数据
# countData <-data.frame(countData,log2FoldChange = res$log2FoldChange, padj = res$padj, pvalue = res$pvalue)

write.table(countData, "/Users/jeffersonchen/programming/YuceBio/YuceBio/Day6/Count_Data.txt", sep = '\t', row.names = F)

write.table(res, "res_output.txt", sep = '\t')

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)

write.table(resdata, "all_des_output.txt", row.names=TRUE, sep = '\t')


# plot MA Graph
ma1 <- plotMA(res, main="DESeq2", ylim=c(2, 2))
ggsave("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day4/R_Code/DeSeq2/Graphs/MA1.pdf", ma1, limitsize = F)


# 读取DESeq结果数据
newdata <- read.delim("res_output.txt", sep = '\t')

par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
# manually save


# *************************
#  Plot Volcano Graphs
p1 <- ggplot(newdata, aes(log2FoldChange, -log(padj,10))) + # -log10 conversion  s
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1
ggsave("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day4/R_Code/DeSeq2/Graphs/Volcano1.pdf", p1)

# 进行分类，根据p值和log值将数据分为三等
new<- newdata %>% 
  mutate(
    Expression = case_when(log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
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

ggsave("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day4/R_Code/DeSeq2/Graphs/Volcano2.pdf", p2, limitsize = F)

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

p3



# ********************
#  Plot HeatMap
library("vsn")
library("pheatmap")

# 取行均最高的一千行进行绘图
select <- order(rowMeans(counts(dds, normalized=T)), decreasing=T)[1:1000]
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[, c("name", "condition")])
hp1<- pheatmap(log2.norm.counts, cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col=df, fontsize=6)
ggsave("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day4/R_Code/DeSeq2/Graphs/HeatMap.pdf", hp1)


vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
ntd <- normTransform(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

# 这行出了问题，但是如果继续跑也能跑出结果（属实奇怪
df <- as.data.frame(colData(dds)[,c("condition","type")]) # cause error, keep running 
hp <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
hp
ggsave("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day4/R_Code/DeSeq2/Graphs/HeatMap2.pdf", hp)


# install.packages("corrplot")
library("reshape2")
library("corrplot")

countData <- countData[, 1:9]
V <- cor(countData, y = NULL, use = "everything", method = c("pearson"))

corr_mat <- round(V, 2)
dist <- as.dist((1-corr_mat)/2)
hc <- hclust(dist)
corr_mat <-corr_mat[hc$order, hc$order]
melted_corr_mat <- melt(corr_mat)

correlation_hm <- ggplot(data = melted_corr_mat, syms=TRUE, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + theme(axis.text.x = element_text(angle = 45, size = 8, margin = margin(t = 15))) + scale_fill_gradient('performance', limits=c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),  low = "black", high = "#FE9000") + ggtitle("Correlation _HeatMap") + xlab("Genes") + ylab("Genes") 

correlation_hm
ggsave("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day4/R_Code/DeSeq2/Graphs/Correlation_Heat_Map.pdf", correlation_hm)

col <- colorRampPalette(c("#4D9DE0", "#E15554"))(20)
heatmap(V, col=col, symm=TRUE,  xlab = "Genes", ylab = "Genes", main = "Correlation_HeatMap(Hie)", margins = c(10, 10))

# ***********************
# 尝试中。。。
rownames(resdata) <- resdata$Row.names
top_genes <- resdata[(resdata$baseMean > 50) & (abs(resdata$log2FoldChange) > 2), ]
top_genes <- top_genes[order(top_genes$log2FoldChange, decreasing = TRUE),]
rld <- rlog(dds, blind=FALSE)
mat <- assay(rld)[rownames(top_genes), sample_names]

colnames(mat) <- sample_names
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale))
colnames(mat.scaled) <- colnames(mat)

num_keep <- 25
row_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)))

l2_val <- as.matrix(top_genes[row_keep,]$log2FoldChange)
colnames(l2_val) <- "logFC"

mean <- as.matrix(top_genes[row_keep,]$baseMean)
colnames(mean) <- "AveExpr"


# BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")

col_logFC <- colorRamp2(c(min(l2_val), 0, max(l2_val)), c("blue", "white", "red"))

col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[row_keep,], cluster_rows = F, column_labels = sample_names, name = "Z-score", cluster_columns = T, row_names_gp = gpar(fontsize = 0) )

# row_labels = top_genes$Row.names[row_keep]
#  cell_fun = function(j, i, x, y, w, h, col){grid.text(round(l2_val[i, j], 2), x, y)})
# h2 <- Heatmap(l2_val, cluster_rows = F, name = "logFC", top_annotation = ha, col = col_logFC)

#  row_labels = top_genes$Row.names[row_keep]
# cell_fun = function(j, i, x, y, w, h, col){grid.text(round(l2_val[i, j], 2), x, y)})
# h3 <- Heatmap(mean, cluster_rows = F, name = "AveExpr", col = col_AveExpr)

# h <- h1 + h2 +h3
h1
# *************************

#  Plot PCA Graph

library("vsn")
library("ggplot2")
vsdata <- vst(dds, blind=FALSE)
pc <- plotPCA(vsdata, intgroup="condition") #using the DESEQ2 plotPCA fxn we can
ggsave("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day4/R_Code/DeSeq2/Graphs/PC.pdf", pc)




