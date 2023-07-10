library(ggplot2)
setwd("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day4/R_Code/Bar_Chart")
# install.packages(c("DESeq2", "Biobase", "gplots", "ggplot2", "RColorBrewer", "pheatmap", "matrixStats", "lattice", "gtools" ,"dplyr"))
my_data <- read.delim("STAR_DATA.txt")
# normalize the data
data <- data.frame( mappedreads=c(my_data[2]), multimapping_reads=c(my_data[3]), unmapped_data=c(my_data[4]))
rownames(data)<- c(my_data$X)

genetype <- c(rep(rownames(data), each = ncol(data)))

Read_type <- rep(c(colnames(data)), time = nrow(data))

Values <- c()
for(i in 1:(ncol(t(data))* ncol(data))) {       # for-loop over columns
  Values <- c(Values, t(data)[i])
}

data <- data.frame(Genetype = genetype, Read_type= Read_type, Values = Values)

bar_chart <- ggplot(data, aes(fill=Read_type, y=Values, x=Genetype)) + 
  geom_bar(position="fill", stat="identity", width =0.9) + scale_y_reverse() + scale_fill_brewer(palette="Paired")+
  theme_minimal() + expand_limits(y=c(0,1)) + theme(axis.text.x = element_text(angle = 45, size = 5))

ggsave("/Users/jeffersonchen/programming/YuceBio/YuceBio/Day4/R_Code/Bar_Chart/Bar_Chart.pdf", bar_chart)
bar_chart
dev.off()

