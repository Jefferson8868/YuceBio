setwd("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day6")
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
# BiocManager::install("tinyarray")
library(tinyarray)
library(tidyverse)

# exp <- read.delim("Count_Data.txt", sep = '\t')
# exp <- trans_exp(exp,mrna_only = T)

LM22.file <- read.delim("./LM22.txt", sep = '\t')
# Changing ID
gene_expression <- read.table("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day6/data_outfpkm.txt", sep = '\t')
sample_names <- c("Gene Symbol", "SRR1382274", "SRR1382275", "SRR1382276", "SRR1382277", "SRR1382278", "SRR1382279", "SRR1382280", "SRR1382281", "SRR1382282")
gene_expression <- gene_expression[-1, ]
colnames(gene_expression) <- sample_names
rownames(gene_expression) <- gene_expression$`Gene Symbol`

selecting_data <- read.delim("Count_Data.txt", sep = '\t')
signif_res <- selecting_data[abs(selecting_data$padj < 0.05) & !is.na(selecting_data$padj) & (abs(selecting_data$log2FoldChange) > 2), ]
signif_genes <- as.character(rownames(signif_res))




gene_expression <- gene_expression %>%
  filter( gene_expression$`Gene Symbol` %in% rownames(signif_res))

signif_genes <-  sub("\\..*", "", signif_genes)

rownames(signif_res) <- signif_genes
rownames(gene_expression) <- signif_genes
# 
# signif_res

# rm(all_genes)


keytypes(org.Hs.eg.db) 

gene_list <- bitr(signif_genes, #数据集
                  fromType="ENSEMBL", #输入为SYMBOL格式
                  toType=c("SYMBOL", "ENTREZID"),  # 转为ENTERZID格式
                  OrgDb="org.Hs.eg.db") #人类 数据库

# all_genes <-  bitr(all_genes, #数据集
#                    fromType="ENSEMBL", #输入为SYMBOL格式
#                    toType=c("SYMBOL", "ENTREZID"),  # 转为ENTERZID格式
#                    OrgDb="org.Hs.eg.db") #人类 数据库
gene_list$SYMBOL

filtered_file <- gene_expression %>%
  filter(rownames(.) %in% gene_list$ENSEMBL)

filtered_file <- filtered_file %>%
  mutate(SYMBOL = gene_list$SYMBOL[match(rownames(filtered_file), gene_list$ENSEMBL)])
# need to manually change to Gene Symbol

filtered_file <- replace(filtered_file, is.na(filtered_file), 0)
filtered_file <- filtered_file[, -1]
names(filtered_file)[names(filtered_file) == "SYMBOL"] <- "Gene Symbol"

filtered_file <- filtered_file[, c('Gene Symbol', 'SRR1382274', 'SRR1382275', 'SRR1382276', 'SRR1382277', 'SRR1382278', 'SRR1382279', 'SRR1382280', 'SRR1382281', 'SRR1382282')]

write.table(filtered_file, "filtered_file.txt", row.names=F, sep = '\t')

Res.results <- CIBERSORT("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day6/LM22.txt" , "filtered_file.txt", perm =100, Q = F)  

# BiocManager::install("immunedeconv")
# install.packages("remotes")
# remotes::install_github("icbi-lab/immunedeconv")
library(ggplot2)
library(immunedeconv)
library(tidyverse)

# write.table(res_cibersort, file="res_cibersort.txt", sep="\t", col.names=T, row.names=F, quote=F)
# save(res_cibersort, file="res_cibersort.Rdata")

# 画图
library(ggplot2)
library(tidyverse)
library(dplyr)

# 其中 nperm 指定的是置换的次数，QN 分位数归一化,如果是芯片设置为 T，如果是测序就设置为 F，测序数据最好是TPM
res_cibersort <- CIBERSORT('LM22.txt','filtered_file.txt', perm = 1000, QN = F)

write.table(res_cibersort, file="res_cibersort.txt", sep="\t", col.names=T, row.names=F, quote=F)
save(res_cibersort, file="res_cibersort.Rdata")

# 画图
library(ggplot2)
library(tidyverse)
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
re <- res_cibersort[, -(23:25)]
k <- apply(re,2,function(x) {sum(x == 0) < nrow(res_cibersort)/2})
table(k)

re2 <- as.data.frame(t(re[,k]))

an = data.frame(group = c("ND", "ND", "ND", "ND", "D", "D", "D", "D", "D"),
                row.names = colnames(gene_expression)[2:10])

pdf("HeatMap.pdf")
hm <- pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
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

ggsave("Histogram.pdf", hg)


nbg <- ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))

ggsave("Unorganized_Boxplot.pdf", nbg)

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

ggsave("Ordered_boxplot.pdf", lg)

dat$Group = ifelse(as.numeric(str_sub(dat$Sample,9,10))<78,"ND","D")
library(ggpubr)
cpg <- ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")

ggsave("Comparsion_Graph.pdf", cpg)



V <- cor(res_cibersort, y = NULL, use = "everything", method = c("pearson"))

corr_mat <- round(V, 2)
dist <- as.dist((1-corr_mat)/2)
hc <- hclust(dist)
corr_mat <-corr_mat[hc$order, hc$order]
melted_corr_mat <- melt(corr_mat)

col <- colorRampPalette(c("#4D9DE0", "#E15554"))(20)
heatmap(V, col=col, symm=TRUE,  xlab = "Genes", ylab = "Genes", main = "Correlation_HeatMap(Hie)", margins = c(10, 10))
