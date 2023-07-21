setwd("/Users/jeffersonchen/programming/YuceBio/YuceBio/Organized_R_Code/R_Code")
library("devtools")
library("magrittr") # needs to be run every time you start R and want to use %>%
library("dplyr")    # alternative
library("ggplot2")
library("dplyr")

# BiocManager::install("DESeq2")

# 加载featureCount 结果

# featureCounts -f -a /home/wangwei/lncRNA_usage_script/gencode.v38.annotation.gtf -t exon -o ./STAR_Comparsion/Feature_Count/without_p.txt --largestOverlap -T 32 ./STAR_Comparsion/*Aligned.sortedByCoord.out.bam
with_p_exon <- read.delim("/Users/jeffersonchen/programming/YuceBio/Feature_Count/With_P_Exon.txt", sep = '\t')

# featureCounts -p -a /home/wangwei/lncRNA_usage_script/gencode.v38.annotation.gtf -o ./STAR_Comparsion/Feature_Count/without_p.txt --largestOverlap -T 32 ./STAR_Comparsion/*Aligned.sortedByCoord.out.bam
with_p_gene <- read.delim("/Users/jeffersonchen/programming/YuceBio/Feature_Count/With_P_Gene.txt", sep = '\t')

# 重新命名数据列名
sample_names <- c("SRR1382274", "SRR1382275", "SRR1382276", "SRR1382277", "SRR1382278", "SRR1382279", "SRR1382280", "SRR1382281", "SRR1382282")
names(with_p_exon)[7:15] <- sample_names
with_p_exon <- with_p_exon[, c(1, 7:15)]
with_p_exon <- aggregate(with_p_exon, .~ Geneid, FUN = "sum")
rownames(with_p_exon) <- with_p_exon$Geneid
with_p_exon <- with_p_exon[, -1]
# with_p_exon <- with_p_exon[ rowSums(with_p_exon) > 1, ]
# 取出所有count值
cd_with_p_exon <- as.matrix(with_p_exon)
cd_with_p_exon <- log(cd_with_p_exon + 1)

names(with_p_gene)[7:15] <- sample_names
with_p_gene <- with_p_gene[, c(1, 7:15)]
with_p_gene <- aggregate(with_p_gene, .~ Geneid, FUN = "sum")
rownames(with_p_gene) <- with_p_gene$Geneid
with_p_gene <- with_p_gene[, -1]
# with_p_gene <- with_p_gene[ rowSums(with_p_gene) > 1, ]
# 取出所有count值
cd_with_p_gene <- as.matrix(with_p_gene)
cd_with_p_gene <- log(cd_with_p_gene + 1)
# write.table(countData, "/Users/jeffersonchen/programming/YuceBio/YuceBio/Day6/Count_Data.txt", sep = '\t', row.names = TRUE)

# 略微调整数据
# 
pdf("./Graphs/FC_Cor_Without_Filter.pdf", width = 15, height = 15)
par(mfrow=c(3,3))
for (i in 1:ncol(cd_with_p_exon)){
  correlation <- cor(cd_with_p_exon[, i], cd_with_p_gene[, i])
  correlation
  # Create a dot plot
  plot(cd_with_p_exon[, i], cd_with_p_gene[, i], pch = 16, xlab = "with_p_exon", ylab = "with_p_gene", main = paste("Correlation Dot Plot", colnames(cd_with_p_exon)[i]))
  
  # Add a correlation line
  abline(lm(cd_with_p_gene[, i] ~ cd_with_p_exon[, i]), col = "red")
  print(paste("Correlation For Sample ", colnames(cd_with_p_exon)[i]))
  print(correlation)
  # Add text for the correlation coefficient
  text(
    x = min(cd_with_p_exon[, i]), y = max(cd_with_p_gene[, i]),
    labels = paste("Correlation:", round(correlation, 2)),
    pos = 4, col = "blue"
  )
}
dev.off()