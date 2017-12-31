#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
# library
library(ggplot2)

myData <- read.delim(args[1], header = F, sep = "\t", fill = TRUE) 

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
