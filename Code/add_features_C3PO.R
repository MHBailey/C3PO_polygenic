library(reshape2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(viridis)
library(data.table)



'%!in%' <- function(x,y)!('%in%'(x,y))
args = commandArgs(trailingOnly=TRUE)

 
#args <- c("Processed_data/{can}.C3PO.combined.v1.txt","Processed_data/CPTAC.{can}.hallmark.scores.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/Meta_table/count.Neoantigen.Mutation.AllCancerType.v2.smoking.tsv")

#args <- c("Processed_data/LUAD.C3PO.combined.v1.txt","Processed_data/CPTAC.LUAD.hallmark.scores.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/Meta_table/count.Neoantigen.Mutation.AllCancerType.v2.smoking.tsv","LUAD","Figures/C3PO.LUAD.Smoking.features.heatmap.pdf","Processed_data/LUAD.Notch.highlights.txt","Processed_data/LUAD.DNARepair.highlights.txt","Processed_data/LUAD.AdaptiveImmune.highlights.txt","Processed_data/LUAD.ChromMod.highlights.txt")

#args <- c("Processed_data/PDAC.C3PO.combined.v1.txt","Processed_data/CPTAC.PDAC.hallmark.scores.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/Meta_table/count.Neoantigen.Mutation.AllCancerType.v2.smoking.tsv","PDAC","Figures/C3PO.PDAC.Smoking.features.heatmap.pdf","Processed_data/PDAC.Notch.highlights.txt","Processed_data/PDAC.DNARepair.highlights.txt","Processed_data/PDAC.AdaptiveImmune.highlights.txt","Processed_data/PDAC.ChromMod.highlights.txt")

#args <- c("Processed_data/LSCC.C3PO.combined.v1.txt","Processed_data/CPTAC.LSCC.hallmark.scores.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/Meta_table/count.Neoantigen.Mutation.AllCancerType.v2.smoking.tsv","LSCC","Figures/C3PO.LSCC.Smoking.features.heatmap.pdf","Processed_data/LSCC.Notch.highlights.txt","Processed_data/LSCC.DNARepair.highlights.txt","Processed_data/LSCC.AdaptiveImmune.highlights.txt","Processed_data/LSCC.ChromMod.highlights.txt")

#args <- c("Processed_data/HNSCC.C3PO.combined.v1.txt","Processed_data/CPTAC.HNSCC.hallmark.scores.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/Meta_table/count.Neoantigen.Mutation.AllCancerType.v2.smoking.tsv","HNSCC","Figures/C3PO.HNSCC.Smoking.features.heatmap.pdf","Processed_data/HNSCC.Notch.highlights.txt","Processed_data/HNSCC.DNARepair.highlights.txt","Processed_data/HNSCC.AdaptiveImmune.highlights.txt","Processed_data/HNSCC.ChromMod.highlights.txt")

#args <- c("Processed_data/UCEC.C3PO.combined.v1.txt","Processed_data/CPTAC.UCEC.hallmark.scores.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/Meta_table/count.Neoantigen.Mutation.AllCancerType.v2.smoking.tsv","UCEC","Figures/C3PO.UCEC.Smoking.features.heatmap.pdf","Processed_data/UCEC.Notch.highlights.txt","Processed_data/UCEC.DNARepair.highlights.txt","Processed_data/UCEC.AdaptiveImmune.highlights.txt","Processed_data/UCEC.ChromMod.highlights.txt")

#args <- c("Processed_data/COAD.C3PO.combined.v1.txt","Processed_data/CPTAC.COAD.hallmark.scores.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/Meta_table/count.Neoantigen.Mutation.AllCancerType.v2.smoking.tsv","COAD","Figures/C3PO.COAD.Smoking.features.heatmap.pdf","Processed_data/COAD.Notch.highlights.txt","Processed_data/COAD.DNARepair.highlights.txt","Processed_data/COAD.AdaptiveImmune.highlights.txt","Processed_data/COAD.ChromMod.highlights.txt")

dat <- fread(args[1])
halls <- fread(args[2])
meta <- fread(args[3])
smoking <- fread(args[4]) 

###So now I want to merge these data together to build some features
meta_lim <- select(meta,c("CASE_ID","Sample","Proteome_Sample_ID"))
smokingm <- merge(smoking,meta_lim,by.x="Sample",by.y="CASE_ID",all.x=T)
cansmoke <- smokingm[which(smokingm$Proteome_Sample_ID %in% halls$ID),]
os <- cansmoke[order(cansmoke$Proteome_Sample_ID),]
baha = columnAnnotation(Count=os$Count,Mutation=os$Mutation,Smoking=os$Smoking)

hs <- dcast(halls,ID~sectors,value.var="absvalue")
hs2 <- hs[which(hs$ID %in% os$Proteome_Sample_ID),]

rownames(hs2) = hs2$ID
hs2$ID <- NULL
mh <- as.matrix(hs2,rownames=rownames(hs2))
#Heatmap(mh)
zmh <- scale(t(mh))


os$Proteome_Sample_ID == colnames(zmh)

pdf(args[6],height=7,width=7,useDingbats=F)
Heatmap(zmh,bottom_annotation=baha,show_row_names = TRUE,show_column_names=FALSE)
dev.off()

#A quick dive into Notch 
myNotch <- hs %>% select(ID,Notch_pathway)
on <- myNotch[order(myNotch$Notch_pathway),]
ons <- merge(on,os,by.x="ID",by.y="Proteome_Sample_ID")
write.table(ons,args[7],sep="\t",quote=F,row.names=F)

#A quick dive into DNA_Repair
myDNA <- hs %>% select(ID,DNA_repair)
on <- myDNA[order(myDNA$DNA_repair),]
ons <- merge(on,os,by.x="ID",by.y="Proteome_Sample_ID")
write.table(ons,args[8],sep="\t",quote=F,row.names=F)

#A quick dive into AdaptiveImmune
myAImmune <- hs %>% select(ID,Adaptive_immunity)
on <- myAImmune[order(myAImmune$Adaptive_immunity),]
ons <- merge(on,os,by.x="ID",by.y="Proteome_Sample_ID")
write.table(ons,args[9],sep="\t",quote=F,row.names=F)


#Chromatin_modification
myChromMod <- hs %>% select(ID,Chromatin_modification)
on <- myChromMod[order(myChromMod$Chromatin_modification),]
ons <- merge(on,os,by.x="ID",by.y="Proteome_Sample_ID")
write.table(ons,args[10],sep="\t",quote=F,row.names=F)


#Now I need an enrichment statistic for tp Chromatin modification analysis. 
plot(density(ons$Smoking)) #This said that there are two major groups above and below .2 

plot(density(ons$Chromatin_modification))
ons$scaleHall <- scale(ons$Chromatin_modification)
plot(density(ons$scaleHall))
ons$scaleSmoke <- scale(ons$Smoking)
ons$highHall <- ifelse(ons$scaleHall > 0, 1, 0)
ons$highSmoke <- ifelse(ons$scaleSmoke > 0, 1, 0)


etest <- chisq.test(ons$highHall,ons$highSmoke)







