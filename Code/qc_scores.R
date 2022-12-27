library(ggplot2)
library(data.table)
library(viridis)
library(dplyr)
library(reshape2)
library(viridis)


'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)
 
#args <- c('Processed_data/CCRCC.tn.prots.txt','Processed_data/CCRCC.pairs.prots.txt','Processed_data/CCRCC.dna.p.effect.mannu.txt','Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_imputed_gene_level.tsv','Data/GeneLists/299.genes.txt','Data/GeneLists/CPTAC-Pancan-summary.txt','Processed_data/CCRCC.cnv.p.effect.mannu.amplification.txt','Processed_data/CCRCC.cnv.p.effect.mannu.deletion.txt', 'CCRCC')


tn = fread(args[1])
pairs = fread(args[2])
dna = fread(args[3])
dna$logP <- -log10(dna$TPVALUE)
prot = fread(args[4])
g299 = fread(args[5])
cptacDs = fread(args[6])
cnva = fread(args[7])
cnva$logP <- -log10(cnva$TPVALUE)
cnvd = fread(args[8])
cnvd$logP <- -log10(cnvd$TPVALUE)
cancertype = args[9]

#work with the mean differences as the metrics because the p-values are off the charts
tn$diff = tn$meanTumor-tn$meanNormal 
tn_sub = tn[which(tn$diff > 3 | tn$diff < -3),]
pdf(args[12],height=6,width=20)
plot(tn[order(tn$diff),]$diff)
dev.off()

#Not work with the pairs to make sure that these differences are across the board
top10 = round(length(unique(pairs$TumorID)) * .1)
nPairs = length(unique(pairs$TumorID))
bot10 = nPairs - top10
pairs_prot <- pairs %>% group_by(Protein) %>% summarize(MeanDelta=mean(DeltaTN),numPos=length(DeltaTN[DeltaTN > 0]),numNeg=length(DeltaTN[ DeltaTN < 0]))
pdf(args[13],height=6,width=20)
plot(pairs_prot$numPos,pairs_prot$MeanDelta)
dev.off()

pairs_sub <- pairs_prot[which(pairs_prot$MeanDelta > 2 | pairs_prot$MeanDelta < -2 | pairs_prot$numPos > bot10 | pairs_prot$numPos < top10),]
 
#NOTE: Keep an eye on these. I also think that there is something to be said about the genes that don't differ at all between tumor and normal (These hiding genes or ultra-immune survailance or fundemental for biology or ???)  
pairs_look <- pairs_prot[which(pairs_prot$MeanDelta < 2 & pairs_prot$numPos > bot10),]
pairs_look <- pairs_prot[which(pairs_prot$MeanDelta > -2 & pairs_prot$numPos < top10),]

#What at about these genes that "cluster tumors into sub populations" That that hover around 40/30 

#NOTE: So the write up for this, get Mean TN differences that are greater than 2,-2 and also pull in all of the genes that move in the same direction for 90% of sample pairs. This hopefully means that they are moving in the same direction just not as much. 


#Now look at the intersection of pair differences and TN differce as "Possible Driver Proteins" and their overlaps with g299
genes <- intersect(tn_sub$Protein, pairs_sub$Protein)
can_genes <- intersect(genes, g299$'Approved symbol')
cptac_genes <- intersect(genes,cptacDs$'Gene Symbol')
#The MAJOR CIS HITTERS - Mutations chage proteins
#"CBFB"  "EGFR"  "HLA-A" "HLA-B" "IDH2"  "PTPRC" "SPTA1"

#Which of the remaining Proteins can be explained via the genotype (My stab at PRS analyis using the DNA NESTs analysis. 

#Now I'm going to build my own little PRS calculator with the NEST DNA data to try to figure it out.

NEST_PROT = NULL
for(g in genes){
    dna_sub = dna[which(dna$Protein == g),]
    cnva_sub = cnva[which(cnva$Protein == g),]
    cnvd_sub = cnvd[which(cnvd$Protein == g),]
    append_dna = dna_sub[which(dna_sub$TPVALUE < 0.5),] %>% select(Protein,NESTv1,COHORT,logP,COHEN_D)


    append_dna$SUBSTRATE = "DNA"
    append_cnva = cnva_sub[which(cnva_sub$TPVALUE < 0.5),] %>% select(Protein,NESTv1,COHORT,logP,COHEN_D)
    append_cnva$SUBSTRATE = "CNVa"
    append_cnvd = cnvd_sub[which(cnvd_sub$TPVALUE < 0.5),] %>% select(Protein,NESTv1,COHORT,logP,COHEN_D)
    append_cnvd$SUBSTRATE = "CNVd"

    NEST_PROT = rbind(NEST_PROT,append_dna)
    NEST_PROT = rbind(NEST_PROT,append_cnva)
    NEST_PROT = rbind(NEST_PROT,append_cnvd)
}

#p <- ggplot(NEST_PROT[NEST_PROT$Protein == g], aes(x=NESTv1,y=Protein,size=logP,color=SUBSTRATE))
#p <- p + geom_point(alpha=0.33)
#p <- p + theme_bw()
#p <- p + theme(axis.text=element_text(size=4),axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
#p <- p + scale_color_viridis_d()
#p 


#p <- ggplot(NEST_PROT[NEST_PROT$Protein == "PLOD2"], aes(x=NESTv1,y=Protein,size=logP,color=SUBSTRATE))
#p <- p + geom_point(alpha=0.33)
#p <- p + theme_bw()
#p <- p + theme(axis.text=element_text(size=4),axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
#p <- p + scale_color_viridis_d()
#p

#p <- ggplot(NEST_PROT, aes(x=NESTv1,y=Protein,size=logP,color=SUBSTRATE))
#p <- p + geom_point(alpha=0.33)
#p <- p + theme_bw()
#p <- p + theme(axis.text=element_text(size=4),axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
#p <- p + scale_color_viridis_d()
#p


p <- ggplot(NEST_PROT[which(NEST_PROT$logP > 2),], aes(x=NESTv1,y=Protein,size=logP,color=SUBSTRATE))
p <- p + geom_point(alpha=0.33)
p <- p + theme_bw()
p <- p + theme(axis.text=element_text(size=4),axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
p <- p + scale_color_viridis_d()
p <- p + ggtitle(cancertype)

pdf(args[10],height=20,width=20)
print(p)
dev.off()

write.table(NEST_PROT,args[11],sep="\t",row.names=F,quote=F) #I could also turn this table into a heatmap of some kind. 


#Figure to build r2 by many figures correleation of the 
#CCRCC = PLOD2 
### https://www.spandidos-publications.com/ijo/48/5/1837
### https://pubmed.ncbi.nlm.nih.gov/26983694/ 
### Cancer status- independant weight. 

#OV =
### Cologen finding finding from Nest
### alcohol dehydrogenase 1B https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5982754/
### HMCN2 


##HEAD and NECK 
#PLOD2x


####UCEC
#LAMB2
#GSTM2
#DCLK1 
#SNCG
#SNCB
#TMOD2
#VWA1
#CAPS
#APCS
#ANGPTL2
#ADIRF























 


