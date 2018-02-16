#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

myObj <- read.delim(args[1], header = F, sep = "\t", fill = TRUE)
myObj

library("ggpubr")

ggscatter(myObj, x = "V2", y = "V4", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", label = "V1", repel = TRUE, 
          xlab = "Breaks Chicken", ylab = "Breaks Finch")


#add = c("none", "reg.line", "loess"), font.label = c(12, "plain"), font.family = "", label = "V1", repel = TRUE, color="V1",
#png(file="mygraphic.png",width=1000,height=850)


myObj <- read.delim(args[2], header = F, sep = "\t", fill = TRUE)
myObj

library("ggpubr")

ggscatter(myObj, x = "V2", y = "V4", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", label = "V1", repel = TRUE, 
          xlab = "Breaks Chicken", ylab = "Breaks Finch")


#add = c("none", "reg.line", "loess"), font.label = c(12, "plain"), font.family = "", label = "V1", repel = TRUE, color="V1",
#png(file="mygraphic.png",width=1000,height=850)

dev.off()
