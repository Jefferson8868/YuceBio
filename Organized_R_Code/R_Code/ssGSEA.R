setwd("/Users/jeffersonchen/programming/YuceBio/YuceBio/Organized_R_Code/R_Code")
library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)
library(rio)


# **** older version
mygroup <- read.delim("./../Data/Original_Version_Info.txt", header = T, sep =' ')
mygroup <- mygroup[, -1]
class(mygroup)
mygroup <- factor(mygroup, levels = c("ND", "D")) %>% relevel(ref = "ND")
# **** newer version
mygroup <- read.delim("./../Data/Sample_Names.txt", header = F)
mygroup <- mygroup[, 1:2]
colnames(mygroup) <- c("Samples", "Condition")
mygroup$Condition <- ifelse(mygroup$Condition == "Nonresponder", "NR", "R")
group <- mygroup[, 2]
row.names(mygroup) <- mygroup[,1]
mygroup <- mygroup[,-1]
mygroup <- factor(dat$Group, levels = c("NR", "R")) %>% relevel(ref = "NR")
class(mygroup)
# ****

# **** newer version
exp <- read.table("./../Results/GEO/filtered.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

# **** older version
exp <- read.table("./../Results/With_P_Exon/ciber_input.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp <- read.table("./../Results/With_P_Gene/ciber_input.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
# 
# exp <- read.table("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day6/ciber_input.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
# 
# 
# exp <- read.table("./../Data/Original_Version_Filtered_Genes.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

# *****************
exp<- as.matrix(exp)


geneset = rio::import("./../Data/gene_set_1.Rdata") 
geneset = geneset[-1]

re <- gsva(exp, geneset, method="ssgsea", mx.diff=FALSE, verbose=FALSE, min.sz > 1) 


ssgsea_score = gsva(exp, geneset, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'

pdf("GSEA_Result.pdf")
draw_boxplot(re,mygroup,color = c("#1d4a9b","#e5171a"))
dev.off()


re
nc[,11:146]
ncol(nc) - length(mygene)+ 1
ncol(nc)
# ***********************
# 需要改变
mygene <- gene_list$SYMBOL  #定义你的目的基因
nc = t(rbind(re,exp[mygene,]))  ;#将你的目的基因匹配到表达矩阵---行名匹配--注意大小写
m = rcorr(nc)$r[1:nrow(re),(ncol(nc)-length(mygene)+1):ncol(nc)]

##计算p值
p = rcorr(nc)$P[1:nrow(re),(ncol(nc)-length(mygene)+1):ncol(nc)]
head(p)


tmp <- matrix(case_when(as.vector(p) < 0.01 ~ "**",
                        as.vector(p) < 0.05 ~ "*",
                        TRUE ~ ""), nrow = nrow(p))

##绘制热图
pdf("HeatMap_ssGSEA.pdf", width = 20, height = 40)
p1 <- pheatmap(t(m),
               display_numbers =t(tmp),
               angle_col =45,
               color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
               border_color = "white",
               cellwidth = 20, 
               cellheight = 20,
               width = 7, 
               height=9.1,
               treeheight_col = 0,
               treeheight_row = 0)
dev.off()





# limma 差异基因
library(limma)
library(edgeR)
design <- model.matrix(~0+group)
rownames(design) = mygroup$Samples
colnames(design) <- c("NR", "R")

DGElist <- DGEList(counts = exp, group = group)


keep_gene <- rowSums(cpm(DGElist) > 1) >= 2
table(keep_gene)
DGElist <- DGElist[keep_gene, ,keep.lib.sizes =FALSE]

DGElist <- calcNormFactors( DGElist )

v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")

fit <- lmFit(v, design)

cont.matrix <- makeContrasts(contrasts = c('NR-R'), levels = design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

nrDEG_limma_voom = topTable(fit2, coef = 'NR-R', n = Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
library(dplyr)
res<-cbind(rownames(nrDEG_limma_voom),nrDEG_limma_voom)
res_1<-res %>% dplyr::filter((logFC>1 | logFC < (-1)) & adj.P.Val < 0.05)
colnames(res_1)[1]<-"Symbol"
