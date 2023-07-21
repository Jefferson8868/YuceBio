setwd("/Users/jeffersonchen/programming/YuceBio/YuceBio/Organized_R_Code/R_Code")
# Gene ID transfer
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(tidyverse)

res <- read.delim("./../Results/With_P_Exon/res_output.txt", sep = '\t')
# exon: 348
# up 159
# down 189
# gene: 347
# up 151
# down 196
signif_res <- res[res$padj < 0.1 & !is.na(res$padj) & (res$log2FoldChange) > 1,]
# signif_res <- res[res$padj < 0.1 & !is.na(res$padj) & (res$log2FoldChange) < 1,]
signif_genes <- as.character(rownames(signif_res))
signif_genes <-  sub("\\..*", "", signif_genes)
all_genes <- as.character(rownames(res))
all_genes <- sub("\\..*", "", all_genes)
rownames(res) <- all_genes

keytypes(org.Hs.eg.db) 

gene_list <- bitr(signif_genes, #数据集
                  fromType="ENSEMBL", #输入为SYMBOL格式
                  toType=c("SYMBOL", "ENTREZID"),  # 转为ENTERZID格式
                  OrgDb="org.Hs.eg.db") #人类 数据库
all_genes <-  bitr(all_genes, #数据集
                   fromType="ENSEMBL", #输入为SYMBOL格式
                   toType=c("SYMBOL", "ENTREZID"),  # 转为ENTERZID格式
                   OrgDb="org.Hs.eg.db") #人类 数据库

all_genes
rownames(res)

filtered_file <- res %>%
  filter(rownames(.) %in% all_genes$ENSEMBL)

filtered_file <- filtered_file %>%
  mutate(ENTREZID = all_genes$ENTREZID[match(rownames(filtered_file), all_genes$ENSEMBL)])
# need to manually change to Gene Symbol

filtered_file <- replace(filtered_file, is.na(filtered_file), 0)

filtered_file

write.table(filtered_file, file = "./../Results/With_P_Gene/filtered_file.txt", sep = "\t", quote = FALSE, row.names = FALSE)

ego <- enrichGO(gene = gene_list$ENTREZID,
                # universe = all_genes$ENTREZID,
                keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                readable = TRUE)





cluster_summary <- data.frame(ego)


pdf("Dotplot.pdf")
dotplot(ego, showCategory=50)
dev.off()
library(enrichplot)
x2 <- pairwise_termsim(ego)

pdf("Emapplot.pdf")
emapplot(x2, showCategory=50)
dev.off()
signif_res_lFC <- signif_res$log2FoldChange
pdf("Centplot.pdf")
cnetplot(ego,
         categorySize="pvalue",
         showCategory = 5,
         color.params= list(signif_res_lFC),
         vertex.label.font=6,
         max.overlaps = 6
)
dev.off() 

pdf("GOgraph.pdf")
plotGOgraph(ego, firstSigNodes=20)
dev.off()

# ego_MF.fil <- simplify(ego_MF)

# ego_ALL.sig <- ego_ALL[ego_ALL$pvalue <= 0.01]
# 过滤后为数据框，不能用自带的参数直接绘制，可以使用ggplot2进行绘制。（暂略）
gene_list$ENTREZID



kk <- enrichKEGG(gene = gene_list$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 0.05)




head(kk,2)
pdf("KEGG_Dotplot.pdf")
dotplot(kk,title="Enrichment KEGG_dot")
dev.off()
kk <- data.frame(kk)
write.csv(summary(kk),"KEGG-enrich.csv",row.names =FALSE)


hsa04974 <- pathview(gene.data = geneList,
                     
                     pathway.id = "hsa04974", #上述结果中的hsa04750通路
                     
                     species = "hsa",
                     
                     limit = list(gene=max(abs(geneList)), cpd=1))

# **************************************
# Perform the GSEA using KEGG gene sets:
# 
# mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
# genes <- getBM(filters="ensembl_gene_id",
#                attributes=c("ensembl_gene_id", "entrezgene_id"),
#                values= all_genes,
#                mart=mart)
# indNA = which(is.na(genes$entrezgene_id))
# genes_noNA <- genes[-indNA,]
# indnodup = which(duplicated(genes_noNA$ entrezgene_id) == F)
# genes_noNA_nodup <- genes_noNA[indnodup,]
# lFC <- res$log2FoldChange[-indNA]
# lFC <- lFC[indnodup]
# names(lFC) <- genes_noNA_nodup$entrezgene_id
# # Sort fold changes in decreasing order
# lFC <- sort(lFC, decreasing = TRUE)
# 
# 
# gseaKEGG <- gseKEGG(geneList = lFC,
#                     organism = "mmu",
#                     nPerm = 1000, # default number permutations
#                     minGSSize = 5, # minimum gene set size
#                     pvalueCutoff = 0.1, # padj cutoff value
#                     verbose = FALSE)
# # Extract the GSEA results
# gseaKEGG_results <- gseaKEGG@result
