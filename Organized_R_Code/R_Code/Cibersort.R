setwd("/Users/jeffersonchen/programming/YuceBio/YuceBio/Organized_R_Code/R_Code")
source("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day6/Assist/CiberSort.R")
# install_github("omnideconv/immunedeconv")
library(e1071)
library(parallel)
library(preprocessCore)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(tidyverse)
library(tinyarray)

# exp <- read.delim("Count_Data.txt", sep = '\t')
# exp <- trans_exp(exp,mrna_only = T)

LM22.file <- read.delim("./../Data/LM22.txt", sep = '\t')

# new version
data <- read.delim("./../Data/GSE135222.hg38.genename.tpm.txt", sep = '\t')
rownames(data) <- data$X1
colnames(data)[1] <- "Gene Symbol"
# data <- data[, -1]
data
write.table(data, "./../Results/GEO/filtered.txt", sep = '\t', row.names = F)

# old version

# exon--gene
# Gene: 14.54% of input gene IDs are fail to map...
# file <- file[rowSums(file) > 1,]
# 18706
# 15975
# Exon: 25.01% of input gene IDs are fail to map...
# file <- file[rowSums(file) > 1,]
# 30086
# 22530
file <- read.delim("data_outfpkm.txt", sep = '\t')

sample_names <- c("SRR1382274", "SRR1382275", "SRR1382276", "SRR1382277", "SRR1382278", "SRR1382279", "SRR1382280", "SRR1382281", "SRR1382282")
file <- aggregate(file, .~ Geneid, FUN = "sum")
rownames(file) <- file$Geneid
file <- file[, -1]
file <- file[rowSums(file) > 1,]
all_genes <- as.character(rownames(file))
all_genes <- sub("\\..*", "", all_genes)
rownames(file) <- all_genes
colnames(file) <- sample_names

keytypes(org.Hs.eg.db) 
all_genes <-  bitr(all_genes, #数据集
                   fromType="ENSEMBL", #输入为SYMBOL格式
                   toType=c("SYMBOL", "ENTREZID"),  # 转为ENTERZID格式
                   OrgDb="org.Hs.eg.db") #人类 数

filtered_file <- file %>%
  filter(rownames(.) %in% all_genes$ENSEMBL)

filtered_file <- filtered_file %>%
  mutate(SYMBOL = all_genes$SYMBOL[match(rownames(filtered_file), all_genes$ENSEMBL)])
# need to manually change to Gene Symbol
nrow(filtered_file)
filtered_file <- replace(filtered_file, is.na(filtered_file), 0)
filtered_file <- filtered_file[, c(10, 1:9)]
colnames(filtered_file)[1] <- "Gene Symbol"
filtered_file <- aggregate(filtered_file, .~ `Gene Symbol`, FUN = mean)
write.table(filtered_file, file = "ciber_input.txt", sep = "\t", quote = FALSE, row.names = FALSE)

res_cibersort <- CIBERSORT('../Data/LM22.txt','ciber_input.txt', perm = 100, QN = F)

write.table(res_cibersort, file="res_cibersort.txt", sep="\t", col.names=T, row.names=T, quote=F)
save(res_cibersort, file="res_cibersort.Rdata")
# write.table(res_cibersort, file="res_cibersort.txt", sep="\t", col.names=T, row.names=T, quote=F)
# save(res_cibersort, file="res_cibersort.Rdata")

# 画图
library(ggplot2)
library(immunedeconv)
library(tidyverse)
library(dplyr)

colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

my_theme <- function(){
  theme(panel.grid = element_blank(),       # 网格线
        panel.border = element_blank(),     # 面板边框
        legend.position="right",            # legend位置
        legend.text = element_text(size=8), # legend内容大小
        legend.title = element_text(size=8),# legend标题大小
        axis.line = element_line(size=1),   # 坐标轴线
        text = element_text(family="Times"),# 文本字体
        axis.text.y = element_text(size = 8,face='bold',color='black'),# y轴标签样式
        axis.text.x = element_text(size = 8,face='bold',color='black',angle=90,hjust=1),        # x轴标签样式，angle=45 倾斜 45 度
        axis.title = element_text(size=10,face="bold"),  # 轴标题
        plot.title = element_text(hjust=0.5,size=10))    # 距，设置绘图区域距离边的据类，上、右、下、左
}  

p1 <- res_cibersort[,1:22] %>% reshape2::melt() %>%
  ggplot(aes(x=Var1,y=value,fill=Var2)) +
  geom_bar(stat='identity') +
  coord_flip()  +
  scale_fill_manual(values =colour ) +
  theme_bw()+ theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+my_theme()

p1
pdf("cibersort.pdf")
p1
dev.off()


library(pheatmap)
# res_cibersort <- read.delim("CIBERSORT-Results.txt", sep = '\t')
# row.names(res_cibersort) <- res_cibersort$Mixture

res_cibersort <- as.matrix(res_cibersort)
re <- res_cibersort[, -(23:25)]
k <- apply(re,2,function(x) {sum(x == 0) < nrow(res_cibersort)/2})
table(k)

re2 <- as.data.frame(t(re[,k]))

# new version
sample_id <- read.delim("./../Data/Sample_Names.txt", sep = '\t', header = F)
colnames(sample_id) <- c("Samples", "Condition", "Dataset")
rownames(sample_id) <- sample_id$Samples
sample_id$Condition <- ifelse(sample_id$Condition == "Nonresponder", "NR", "R") %>% factor( levels = c("NR", "R")) %>% relevel(ref = "NR")
an = data.frame(group = sample_id$Condition,
                row.names = sample_id$Samples)


# old version
sample_id <- read.delim("./../Data/Original_Version_Info.txt", sep = ' ', header = T)
colnames(sample_id) <- c("Samples", "Condition")
sample_id$Condition <- sample_id$Condition %>% factor( levels = c("ND", "D")) %>% relevel(ref = "ND")
an = data.frame(group = sample_id$Condition,
                row.names = sample_id$Samples)

pdf("HeatMap_C.pdf")
hm <- pheatmap(re2,scale = "row",
         show_colnames = T,
         annotation_col = an,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_cols = FALSE,
         main = "CiberSort VS. Sample Groups"
         )
hm
dev.off()

library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

hg <- ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))

hg
ggsave("Histogram.pdf", hg, width = 10, height = 8)


nbg <- ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))

nbg
ggsave("Unorganized_Boxplot.pdf", nbg, width = 10, height = 8)

a = dat %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)

dat$Cell_type = factor(dat$Cell_type,levels = a)

lg <- ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))

lg
ggsave("Ordered_boxplot.pdf", lg, width = 11, height = 8)

# new version
dat$Group = ifelse(sample_id[dat$Sample, 2] == "NR", "NR", "R")
dat$Group <- factor(dat$Group, levels = c("NR", "R")) %>% relevel(ref = "NR")

# old version
rownames(sample_id) <- sample_id$Samples
dat$Group = ifelse(sample_id[dat$Sample, 2] == "ND", "ND", "D")
dat$Group <- factor(dat$Group, levels = c("ND", "D")) %>% relevel(ref = "ND")

library(ggpubr)
cpg <- ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")


ggsave("Comparsion_Graph.pdf", cpg,  width = 11, height = 8)



V <- cor(res_cibersort, y = NULL, use = "everything", method = c("pearson"))

corr_mat <- round(V, 2)
dist <- as.dist((1-corr_mat)/2)
hc <- hclust(dist)
corr_mat <-corr_mat[hc$order, hc$order]
melted_corr_mat <- melt(corr_mat)

col <- colorRampPalette(c("#4D9DE0", "#E15554"))(20)
pdf("Testing.pdf")
heatmap(V, col=col, symm=TRUE,  xlab = "Genes", ylab = "Genes", main = "Correlation_HeatMap(Hie)", margins = c(10, 10))
dev.off()


data <- data[, -1]
V <- cor(data, y = NULL, use = "everything", method = c("pearson"))

corr_mat <- round(V, 2)
dist <- as.dist((1-corr_mat)/2)
hc <- hclust(dist)
corr_mat <-corr_mat[hc$order, hc$order]
melted_corr_mat <- melt(corr_mat)


correlation_hm <- ggplot(data = melted_corr_mat, syms=TRUE, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + theme(axis.text.x = element_text(angle = 45, size = 8, margin = margin(t = 15))) + scale_fill_gradient('performance', limits=c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),  low = "black", high = "#FE9000") + ggtitle("Correlation _HeatMap") + xlab("Genes") + ylab("Genes") 

ggsave("Correlation_HM.pdf", correlation_hm, width = 11, height = 8)

