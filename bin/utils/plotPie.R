#!/usr/bin/env Rscript
# Barplot
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
pdf(args[3], onefile=TRUE)
myMatrix <- read.table(args[1], sep="\t", header=TRUE)
myMatrix
myPie <- ggplot(myMatrix , aes(x= args[2], y=cnt, fill=names))+ 
geom_bar(width = 1, stat = "identity")
myPie
