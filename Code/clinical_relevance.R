library(data.table)
library(reshape2)
library(ggplot2)
library(stringr)
library(viridis)
library(dplyr)



'%!in%' <- function(x,y)!('%in%'(x,y))
args = commandArgs(trailingOnly=TRUE)

#args = c("Processed_data/CPTAC.BRCA.hallmark.scores.txt","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","BRCA",'Processed_data/immune_correlations.BRCA.txt','Figures/C3POxImmune.BRCA.pdf')

#args = c("Processed_data/CPTAC.LUAD.hallmark.scores.txt","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","LUAD",'Processed_data/immune_correlations.LUAD.txt','Figures/C3POxImmune.LUAD.pdf')

#args = c("Processed_data/CPTAC.CCRCC.hallmark.scores.txt","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","CCRCC",'Processed_data/immune_correlations.CCRCC.txt','Figures/C3POxImmune.CCRCC.pdf')

 
#args = c("Processed_data/CPTAC.OV.hallmark.scores.txt","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","OV",'Processed_data/immune_correlations.OV.txt','Figures/C3POxImmune.OV.pdf')

#args = c("Processed_data/CPTAC.PDAC.hallmark.scores.txt","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","PDAC",'Processed_data/immune_correlations.PDAC.txt','Figures/C3POxImmune.PDAC.pdf')

#args = c("Processed_data/CPTAC.COAD.hallmark.scores.txt","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","COAD",'Processed_data/immune_correlations.COAD.txt','Figures/C3POxImmune.COAD.pdf')

#args = c("Processed_data/CPTAC.PANCAN.hallmark.scores.txt","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","PANCAN",'Processed_data/immune_correlations.PANCAN.txt','Figures/C3POxImmune.PANCAN.pdf')

#args = c("Processed_data/CPTAC.UCEC.hallmark.scores.txt","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","UCEC",'Processed_data/immune_correlations.UCEC.txt','Figures/C3POxImmune.UCEC.pdf')
 
#args = c("Processed_data/CPTAC.GBM.hallmark.scores.txt","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","GBM",'Processed_data/immune_correlations.GBM.txt','Figures/C3POxImmune.GBM.pdf')


hall = fread(args[1])
clust = fread(args[2])
can = args[3]
outf = args[4]
ofig = args[5]

dh <- dcast(hall,value.var="value",ID~sectors)
hc <- merge(dh,clust,by.x="ID",by.y="Proteome_Sample_ID")


#CLEAN UP THE NAs
hc[hc=="n/a"] <- NA
hc[hc=="Not Applicable"] <- NA
hc[hc=="NA"] <- NA


hcc <- hc %>% select(where(~mean(is.na(.))< 0.1))

mine = c("ID", "Adaptive_immunity"  ,"Angiogenesis" ,"Cancer_driver_genes"  ,"Cell_cycle_and_apoptosis" ,"Chromatin_modification" ,"DNA_repair" ,"Epithelial.mesenchymal_transition"  ,"Extracellular_matrix" ,"Hedgehog_pathway" ,"Humoral_immunity" ,"Inflammation" ,"Innate_immunity"  ,"JAK.STAT_pathway" ,"MAPK_pathway" ,"Metastasis" ,"Notch_pathway"  ,"PI3K_pathway" ,"RAS_pathway"  ,"TGFB_pathway" ,"Transcriptional_misregulation"  ,"Wnt_pathway", "Immune_subtype_pancan" )

immune =  c("ID", "Adaptive_immunity","Humoral_immunity" ,"Inflammation", "Innate_immunity","Immune_subtype_pancan")

getmine <- data.frame(hcc %>% select(all_of(mine)))
mh <- melt(getmine,variable.var=c("ID","Immune_subtype_pancan"))

p <- ggplot(mh[which(mh$variable %in% immune),],aes(x=variable,y=value,color=Immune_subtype_pancan))
p <- p + geom_boxplot()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
p


toremove = c("order" ,"CASE_ID" ,"Sample" ,"Other.label" ,"Proteome_Normal_Sample_ID" ,"Phosphoproteome_Sample_ID" ,"Phosphoproteome_Normal_Sample_ID" ,"Acetylome_Sample_ID" ,"Acetylome_Normal_Sample_ID" ,"RNA_Tumor" ,"RNA_Normal" ,"Germline_WGS" ,"WXS" ,"tumor_code")

tokeep = colnames(hcc)[which(colnames(hcc) %!in% toremove)]

hccc <- data.frame(hcc %>% select(all_of(tokeep)))


rownames(hccc) <- hccc$ID
numhc <- dplyr::select_if(hccc, is.numeric)
mn <- as.matrix(numhc)

nona <- na.omit(mn)
mcor <- cor(nona)

pdf(ofig,height=10,width=10,useDingbats=F)
corrplot(mcor, method="circle",tl.cex =.3)
dev.off()


yo = data.frame(nona)
plot(yo$)





