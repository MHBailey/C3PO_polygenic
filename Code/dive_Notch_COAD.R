library(data.table)
library(stringr)
library(reshape2)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyr)
library(ComplexHeatmap)

args = commandArgs(trailingOnly=TRUE)

#args = c("Data/Hallmarks/nanostring.gl.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/NESTs/TheNEST.csv","Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Somatic_cnv/Broad_pipeline_wxs_formatted/GISTIC_all_thresholded.by_genes.tsv","Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_raw_gene_level_v1.3.tsv","Processed_data/COAD.PolyRisk.dna.v2.txt","Processed_data/COAD.PolyRisk.amp.v2.txt","Processed_data/COAD.PolyRisk.del.v2.txt","Processed_data/COAD.cnv.p.effect.mannu.amplification.txt","Processed_data/COAD.cnv.p.effect.mannu.deletion.txt","Processed_data/COAD.dna.p.effect.mannu.txt","X01BR043","Notch_pathway","Processed_data/COAD.Notch.highlights.txt","Processed_data/COAD.C3PO.combined.v1.txt")


halls = fread(args[1])
myhall = c("Genes",args[14])

myh = halls %>% select(all_of(myhall))
genesInHall = myh[which(myh[,2] == "Y"),]$Genes

notch <- fread(args[15])


#Get Proteins in the list
prot <- fread(args[6])
protgenes <- prot$external_gene_name

#Intersect hallmarks with Proteins
all_genes <- intersect(protgenes,genesInHall)


ph = prot %>% select("external_gene_name",all_of(notch$ID))
phg <- ph[which(ph$external_gene_name %in% all_genes),]
mp <- na.omit(melt(phg))
dp <- na.omit(dcast(mp,external_gene_name~variable))

mm <- merge(mp,notch,by.x="variable",by.y="ID",all.x=T)


p <- ggplot(mm,aes(external_gene_name,value,fill=Signature))
p <- p + geom_boxplot()
p


plot(log10(mm[which(mm$external_gene_name == "NOTCH2"),]$Mutation),mm[which(mm$external_gene_name == "NOTCH2"),]$value)

#Now lets see which gene is pushing it in the scores not just the proteins... 
comb <- fread(args[16])

cg <- comb %>% select("V1",any_of(all_genes))
cgm <- merge(cg,notch,by.x="V1",by.y="ID")

plot(log10(cgm$Mutation),cgm$HDAC2)
plot(log10(cgm$Mutation),cgm$NOTCH2)
plot(log10(cgm$Mutation),cgm$NOTCH3)

plot(log10(mm[which(mm$external_gene_name == "HDAC2"),]$Mutation),mm[which(mm$external_gene_name == "HDAC2"),]$value)
plot(log10(mm[which(mm$external_gene_name == "NOTCH2"),]$Mutation),mm[which(mm$external_gene_name == "NOTCH2"),]$value)
plot(log10(mm[which(mm$external_gene_name == "NOTCH3"),]$Mutation),mm[which(mm$external_gene_name == "NOTCH3"),]$value)

summary(lm(log10(mm[which(mm$external_gene_name == "NOTCH3"),]$Mutation)~mm[which(mm$external_gene_name == "NOTCH3"),]$value))
summary(lm(log10(mm[which(mm$external_gene_name == "NOTCH2"),]$Mutation)~mm[which(mm$external_gene_name == "NOTCH2"),]$value
rn = dp$external_gene_name
rownames(dp) = dp$external_gene_nam
dm <- as.matrix(dp[-1],rownames.force=T)

os <- unique(data.frame(mm %>% select(variable,Mutation)))

os$variable == colnames(dm)
os$m200 <- ifelse(os$Mutation > 200,1,0)
os$m100 <- ifelse(os$Mutation > 100,1,0)
baha = columnAnnotation(Mutation=log10(os$Mutation), Morethan200=os$m200,Morethan100=os$m100 )
Heatmap(dm,bottom_annotation=baha,show_row_names = TRUE,show_column_names=FALSE)

#NOW I WANT TO ADD A CNV TRACT to maybe caputre some of the CNV HEAVY NOT EVENTS 



meta = fread(args[2])
ucec <- meta[which(meta$cohort == "UCEC"),]
cnvid = ucec$CASE_ID
cnv <- fread(args[5])
cnv_keep = c("Gene Symbol",cnvid)
subcnv <- cnv %>% select(all_of(cnv_keep))

sc <- subcnv[-1]
sc[sc > 0] <- 1
sc[sc < 0] <- 1
sc[sc == 0] <- 0 
df2 <- mutate_all(sc[-1], function(x) as.numeric(as.character(x)))
cnvcnt <- data.frame(cnvcnt=colSums(df2[,-1]))
cnvcnt$CASE_ID = rownames(cnvcnt)
ucec_keep <- ucec %>% select(CASE_ID, Proteome_Sample_ID)
cm <- merge(cnvcnt,ucec_keep,by="CASE_ID")

cm$Proteome_Sample_ID == colnames(dm)
cm$cnv7500 <- ifelse(cm$cnvcnt > 7500, 1,0)

baha = columnAnnotation(Mutation=log10(os$Mutation), Morethan200=os$m200,Morethan100=os$m100,CNVmuts=cm$cnvcnt,CNV7500=cm$cnv7500)
Heatmap(dm,bottom_annotation=baha,show_row_names = TRUE,show_column_names=FALSE)

#Now lets come up with a test
tmm <- mm[which(mm$external_gene_name == "NOTCH2"),]
plot(density(log10(tmm$Mutation))) #this provided two nice groups that I could split into to based on counts.
abline(v=275)
t.test(x=tmm[which(tmm$Mutation >= 275),]$value,y=tmm[which(tmm$Mutation<275),]$value)

