setwd("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day7/R_Code")
library(dplyr)
library(ggpubr)
library(tidyverse)

ciber_result <- read.delim("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day7/Data/res_cibersort.txt", sep = '\t')
ciber_result <- ciber_result[, 1:22]

sample_names <- c("SRR1382274", "SRR1382275", "SRR1382276", "SRR1382277", "SRR1382278", "SRR1382279", "SRR1382280", "SRR1382281", "SRR1382282")


rownames(ciber_result) <- sample_names
ciber_result$Samples <- rownames(ciber_result)

new_data <- ciber_result %>% 
  pivot_longer(!Samples, names_to = "Cells", values_to = "Proportion")
new_data$Group <- ifelse(as.numeric(str_sub(new_data$Samples,9,10))<78,"ND","D")

a = new_data %>% 
  group_by(Cells) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cells)

new_data$Cells = factor(new_data$Cells,levels = a)



pdf("plot.pdf", width = 12, height = 16)
p <- ggboxplot(new_data, x = "Group", y = "Proportion",
               color = "Group", palette = "jco",
               add = "jitter",
               facet.by = "Cells", short.panel.labs = FALSE)


# Use only p.format as label. Remove method name.
p + stat_compare_means(method = "wilcox.test", label = "p.format")

dev.off()

library(RColorBrewer)
# install.packages("ggdist")
library(ggdist)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

pdf("Comparsion.pdf", width = 12, height = 8)
p <- ggboxplot(new_data, x = "Cells", y = "Proportion",
               color = "Group", palette = "jco",
               add = "jitter")
p + stat_compare_means(aes(group = Group), label = "p.signif") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom") 
  # ggdist::stat_halfeye(
  #   adjust = 0.5, 
  #   justification = .2,
  #   .width = 0,
  #   point_colour = NA
  # ) + 
  # ggdist::stat_dots(
  #   side = "left",
  #   justification = 1.1,
  #   
  # )

dev.off()

