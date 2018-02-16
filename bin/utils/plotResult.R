#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

#REF1 here -----------------------------------------------------
myObj <- read.delim(args[1], header = F, sep = "\t", fill = TRUE)

library("ggpubr")

ggscatter(myObj, x = "V2", y = "V5", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", label = "V1", repel = TRUE, 
          xlab = "Breaks in Chicken Reference", ylab = "Breaks in Finch Reference")


#add = c("none", "reg.line", "loess"), font.label = c(12, "plain"), font.family = "", label = "V1", repel = TRUE, color="V1",
#png(file="mygraphic.png",width=1000,height=850)

#REF2 here ------------------------------------------------------
myObj <- read.delim(args[2], header = F, sep = "\t", fill = TRUE)

library("ggpubr")

ggscatter(myObj, x = "V2", y = "V5", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", label = "V1", repel = TRUE, 
          xlab = "Breaks in Chicken Reference", ylab = "Breaks in Finch Reference")


#add = c("none", "reg.line", "loess"), font.label = c(12, "plain"), font.family = "", label = "V1", repel = TRUE, color="V1",
#png(file="mygraphic.png",width=1000,height=850)

# library
library(ggplot2)

myData <- read.delim(args[3], header = F, sep = "\t", fill = TRUE) 

#                      V1       V2  V3
#1          columba_livia    Break  77
#2          columba_livia Detected  77
#3          columba_livia     Real  73
#4  corvus_brachyrhynchos    Break  28
#5  corvus_brachyrhynchos Detected  24
#6  corvus_brachyrhynchos     Real  24

# Grouped
p1<-ggplot(myData, aes(fill=V2, y=V3, x=V1)) + 
    geom_bar(position="dodge", stat="identity")
p1 + coord_flip()

# Stacked
p2<-ggplot(myData, aes(fill=V2, y=V3, x=V1)) + 
    geom_bar( stat="identity")
 
p2 + coord_flip() 

# Stacked Percent
p3<-ggplot(myData, aes(fill=V2, y=V3, x=V1)) + 
    geom_bar( stat="identity", position="fill")

p3 + coord_flip()


myData <- read.delim(args[4], header = F, sep = "\t", fill = TRUE) 

#                      V1       V2  V3
#1          columba_livia    Break  77
#2          columba_livia Detected  77
#3          columba_livia     Real  73
#4  corvus_brachyrhynchos    Break  28
#5  corvus_brachyrhynchos Detected  24
#6  corvus_brachyrhynchos     Real  24

# Grouped
p1<-ggplot(myData, aes(fill=V2, y=V3, x=V1)) + 
    geom_bar(position="dodge", stat="identity")
p1 + coord_flip()

# Stacked
p2<-ggplot(myData, aes(fill=V2, y=V3, x=V1)) + 
    geom_bar( stat="identity")
 
p2 + coord_flip() 

# Stacked Percent
p3<-ggplot(myData, aes(fill=V2, y=V3, x=V1)) + 
    geom_bar( stat="identity", position="fill")

p3 + coord_flip()
dev.off()
