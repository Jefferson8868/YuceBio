setwd("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day7/R_Code")
library(dplyr)


filtered_file <- read.delim("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day6/data_outfpkm.txt", sep = '\t' )  
sample_names <- c("SRR1382274", "SRR1382275", "SRR1382276", "SRR1382277", "SRR1382278", "SRR1382279", "SRR1382280", "SRR1382281", "SRR1382282")

rownames(filtered_file) <- filtered_file[, 1]
filtered_file <- filtered_file[, -1]
filtered_file<- filtered_file[ rowSums(filtered_file) > 1, ]
rownames(filtered_file) <-  sub("\\..*", "", rownames(filtered_file))
colnames(filtered_file) <- sample_names

Gene_Names = read.delim("Gene_Names.txt", sep = '\t')

filtered_file <- filtered_file %>% filter(rownames(.) %in% Gene_Names$ENSEMBL)

filtered_file <- filtered_file %>%
  mutate(SYMBOL = Gene_Names$SYMBOL[match(rownames(filtered_file), Gene_Names$ENSEMBL)])

filtered_file <- replace(filtered_file, is.na(filtered_file), 0)
rownames(filtered_file) <- filtered_file$SYMBOL
filtered_file <- filtered_file[,1:9]



library(ggpubr)
# for (i in 1:ncol(filtered_file))
# {
#   colname <- colnames(filtered_file)[i]
#   hist(filtered_file[, i], prob = TRUE, main = paste("Density plot of ", colname))
#   xfit <- seq(min(filtered_file[, i]), max(filtered_file[, i]), length = 20)
#   yfit <- dnorm(xfit, mean(filtered_file[, i]), sd(filtered_file[, i]))
#   lines(xfit, yfit, col = "red", lwd = 2)
#   lines(density(filtered_file[, i]), col = "blue", lwd = 2)
#  # print( shapiro.test((filtered_file[, i])))
# }

# ggqqplot(filtered_file[, 6], color = "blue", main = paste("Normal Q-Q Plot", colname), ylim = c(0, 5000))

# for (i in 1:ncol(filtered_file))
# {
#   for (j in 1:ncol(filtered_file))
#   {
#     x <- t.test(filtered_file[, i], y = c(filtered_file[, j]), alternative = "less")
#   if(x$p.value <= 0.05)
#   {
#     print(i)
#     print(j)
#     print(x)
#   }
#   }
#   
# }

# for (i in 1:ncol(filtered_file))
# {
# print(shapiro.test(filtered_file[, i]))
# }
# 
# qqnorm(filtered_file[, i])
# qqline(filtered_file[, i])
# 



for (i in 1:ncol(filtered_file))
{
  for (j in 1:ncol(filtered_file))
  {
    x <- wilcox.test(filtered_file[, i], y = c(filtered_file[, j]), alternative = "two.sided", paired = FALSE)
    if(x$p.value <= 0.05)
    {
    
      print(i)
      print(j)
      print(x$alternative)
      print(x)
      
    }
  }
}

library(tidyverse)

# Assuming your data is stored in a data frame called 'original_data'
filtered_file$Gene <- rownames(filtered_file)

new_data <- filtered_file %>% 
  pivot_longer(!Gene, names_to = "Samples", values_to = "Counts")


# new_data$Samples_In_Num <- as.numeric(str_sub(new_data$Samples, 9, 10))-74
new_data$Group <- ifelse(as.numeric(str_sub(new_data$Samples,9,10))<78,"ND","D")
# group_by(new_data, Group)
new_data$Group <- relevel(new_data$Group, ref = "ND")
levels(new_data$Group)
wilcox.test(Counts ~ Group, data = new_data, alternative = "greater")
# wilcox.test(new_data$Samples_In_Num, new_data$Counts, alternative = "two.sided")

library("ggpubr")
ggboxplot(new_data, x = "Group", y = "Counts", 
          color = "Group", palette = c("#00AFBB", "#E7B800"),
          ylab = "Counts", xlab = "Groups", ylim = c(0, 200))


ggboxplot(new_data,x = "Group",y = "Counts",color = "Group",add = "jitter") + stat_compare_means()
group=levels(factor(new_data$Group)) 
#将Treatment转换成因子型变量
new_data$Group =factor(new_data$Group, levels=group)
#获得Treatment中元素之间的组合，即：设置比较组（将所有实验组分成两两一组进行后续比较）
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

# group=levels(factor(new_data$Group))
# new_data$Group=factor(new_data$Group, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
mytheme <- theme(axis.title=element_text(size=30),
                 axis.title.x=element_text(size=15),
                 axis.title.y=element_text(size=15),
                 legend.title=element_text(size=15),
                 legend.text=element_text(size=15),
                 axis.text.x=element_text(size=15))
#绘图
#stat_compare_means(comparisons = my_comparisons):指定需要进行比较以及添加p-value、显著性标记的组
boxplot=ggboxplot(new_data, 
                  x="Group",   
                  y="Counts", 
                  ylab="Counts",
                  legend.title=x,
                  color="Group",
                  palette = c("blue","red"),
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons) +
  scale_y_continuous(limits=c(0,1.05)) +
  mytheme

boxplot

library("ggpubr")
data(ToothGrowth)

df3 <- ToothGrowth
ggboxplot(new_data, x = "Group", y = "Counts")+
  stat_compare_means(ref.group = "0.5", 
                     method = "wilcox.test",
                     method.args = list(alternative = "less"))

#绘制boxplot
pdf(file="Boxplot.pdf",width=5,height=4.5)
boxplot <- ggboxplot(new_data,
                  x="Group",
                  y="Counts",
                  ylim = c(0, 100),
                  color="Group",
                  legend.title=x,
                  palette = c("#00AFBB", "#E7B800"),
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "more")) #默认是wilcox.test，可换成t.test,kruskal.test，anova，下面同理
boxplot
print(boxplot)
dev.off()



