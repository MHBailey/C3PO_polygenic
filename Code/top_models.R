##################
#Author Matt Bailey
#Last updated 2022-Jun-28
##################


#Objective
#Evaluate whether the top predicted genes are also the genes that performed well (better) in the validation 

library(data.table)
library(reshape2)
library(ggplot2)
library(stringr)
library(viridis)
library(dplyr)


'%!in%' <- function(x,y)!('%in%'(x,y))
args = commandArgs(trailingOnly=TRUE)

#Breast DNA:
args = c("Validation/BRCA.SampleCorrelations.dna.v1.txt","Processed_data/BRCA.SampleCorrelations.dna.v1.txt")

valid <- fread(args[1])
orig <- fread(args[2])

vquant <- quantile(valid$Coef)[4]
oquant <- quantile(orig$Coef)[4]

vgenes <- valid[which(valid$Coef >= vquant),]$Protein
ogenes <- orig[which(orig$Coef >= oquant),]$Protein

#Shrink this list to the same size and same overlapiping genes
all_genes <- intersect(valid$Protein,orig$Protein)

new_valid <- valid[which(valid$Protein %in% all_genes),]
new_orig <- orig[which(orig$Protein %in% all_genes),]

new_valid$rank = rank(new_valid$Coef)
new_orig$rank <- rank(new_orig$Coef)


new_valid$Protein == new_orig$Protein

o_valid <- new_valid[order(new_valid$Protein),]
o_orig <- new_orig[order(new_orig$Protein),]

o_valid$rank = rank(o_valid$Coef)
o_orig$rank <- rank(o_orig$Coef)


res = cor.test(o_valid$rank,o_orig$rank,method="kendall")

