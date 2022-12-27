library(reshape2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(viridis)
library(data.table)
library(ggplot2)

uniqString <- function(x){
    yo <- unique(unlist(strsplit(x, ";")))
    return(paste(yo,collapse=";"))
}
'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)

#args <- c("Processed_data/CPTAC.BRCA.hallmark.scores.txt","BRCA",'Figures/CPTAC.BRCA.hall_heatmap.v2.pdf')
#args <- c("Processed_data/CPTAC.LSCC.hallmark.scores.txt","LSCC",'Figures/CPTAC.LSCC.hall_heatmap.v2.pdf')
#args <- c("Processed_data/CPTAC.UCEC.hallmark.scores.txt","UCEC",'Figures/CPTAC.UCEC.hall_heatmap.v2.pdf')
#args <- c("Processed_data/CPTAC.GBM.hallmark.scores.txt","GBM",'Figures/CPTAC.GBM.hall_heatmap.v2.pdf')

hallscore <- fread(args[1])
cancer = args[2]


hs <- dcast(hallscore,ID~sectors,value.var="absvalue")
rownames(hs) = hs$ID
hs$ID <- NULL
mh <- as.matrix(hs,rownames=rownames(hs))
#Heatmap(mh)
zmh <- scale(t(mh))

#write.table(mh,"ForCollin/LSCC.raw.hallmarks.txt",sep="\t",quote=F)
#write.table(zmh,"ForCollin/LSCC.scaled.hallmarks.txt",sep="\t",quote=F)
#write.table(mh,"ForCollin/UCEC.raw.hallmarks.txt",sep="\t",quote=F)
#write.table(zmh,"ForCollin/UCEC.scaled.hallmarks.txt",sep="\t",quote=F)

pdf(args[3],height=10,width=10,useDingbats=F)
Heatmap(zmh)
dev.off()

mean(zmh)
sd(zmh)

#Annotations to add 
#Multiomic_subtype_original_annotation Immune_subtype_cohort sex age ESTIMATE_ImmuneScore predicted_ancestry medical_history Overall survival, days
smeta <- meta %>% select("Proteome_Sample_ID","cohort","Multiomic_subtype_original_annotation","Immune_subtype_cohort","sex","age","ESTIMATE_ImmuneScore","predicted_ancestry","medical_history/bmi","Overall survival, days")

bsm <- smeta[which(smeta$cohort == "GBM"),]
baha=columnAnnotation(subtype=bsm$Multiomic_subtype_original_annotation,immune_cohort=bsm$Immune_subtype_cohort,sex=bsm$sex,age=bsm$age,estimateImmune=bsm$ESTIMATE_ImmuneScore,ancestry=bsm$predicted_ancestry,history=as.numeric(bsm$"medical_history/bmi"),survival=as.numeric(bsm$"Overall survival, days"),na_col = "black")

Heatmap(zmh,bottom_annotation=baha,show_row_names = FALSE)



