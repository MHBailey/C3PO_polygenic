library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)


'%!in%' <- function(x,y)!('%in%'(x,y))
args = commandArgs(trailingOnly=TRUE)

#This is just going to be some scratch code for one examples before 
#args = c("Processed_data/LUAD.PolyRisk.dna.v2.txt","Processed_data/LUAD.PolyRisk.amp.v2.txt","Processed_data/LUAD.PolyRisk.del.v2.txt","LUAD",'Processed_data/LUAD.C3PO.combined.v1.txt')


prs_dna = data.frame(fread(args[1]))
rownames(prs_dna) <- prs_dna$Proteome_Sample_ID 
prs_dna$Proteome_Sample_ID <- NULL 
prs_amp = data.frame(fread(args[2]))
rownames(prs_amp) <- prs_amp$Proteome_Sample_ID
prs_amp$Proteome_Sample_ID <- NULL 
prs_del = data.frame(fread(args[3]))
rownames(prs_del) <- prs_del$Proteome_Sample_ID
prs_del$Proteome_Sample_ID <- NULL
cancer = args[4]

u1 <- intersect(colnames(prs_dna),colnames(prs_amp))
u2 <- intersect(u1,colnames(prs_del))

dna2 <- prs_dna %>% select(all_of(u2))
amp2 <- prs_amp %>% select(all_of(u2))
del2 <- prs_del %>% select(all_of(u2))

outdf <- data.frame(as.matrix(dna2) + as.matrix(amp2) + as.matrix(del2))

write.table(outdf,args[5],sep="\t",quote=F,row.names=T)
