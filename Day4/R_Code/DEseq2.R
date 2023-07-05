#install.packages("htmltools")
library(htmltools)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

library( "DESeq2" )
library(ggplot2)

countsName <- "http://bioconnector.org/workshops/data/airway_scaledcounts.csv"
download.file(countsName, destfile = "airway_scaledcounts.csv", method = "auto")

countData <- read.csv('airway_scaledcounts.csv', header = TRUE, sep = ",")
head(countData)
