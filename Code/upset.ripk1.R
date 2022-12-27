library(data.table)
library(UpSetR)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)


args <- c("SideStories/RIPK1.upsetR.dat.txt")
dat <- fread(args[1])

upset(dat,sets=c("Basal","0","1","2","-1","High","Low","Med","MedHigh","MedLow"),boxplot.summary = c("Protein"))
upset(dat,nsets = 15,boxplot.summary = c("Protein"))


dat2 <- dat[which(dat$Basal == 1),]
upset(dat2,nsets = 15,boxplot.summary = c("Protein"))
