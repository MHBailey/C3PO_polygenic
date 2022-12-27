library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)
library(viridis)

'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)


#args <- c("Data/NESTs/TheNEST.csv","Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_imputed_gene_level.tsv","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt")



nest = fread(args[1],header=F,sep=",")
maf = fread(args[2])
prot = fread(args[3])
meta = fread(args[4])

#STEP1 (reduces the maf to missense, nonsense, frameshifts)
muts2keep <- c("Missense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation")
cmaf <- maf[which(maf$Variant_Classification %in% muts2keep),]


#STEP2 (make a list of samples with and without a mutation in a NEST)
NEST = "NEST:53"
PROT = "NNMT"
nestgs <-  unlist(strsplit(nest[which(nest$V1 == NEST),]$V3," "))

mutsamps = cmaf[which(cmaf$Hugo_Symbol %in% nestgs),]$Tumor_Sample_Barcode
wtsamps = unique(cmaf$Tumor_Sample_Barcode[which(cmaf$Tumor_Sample_Barcode %!in% mutsamps)])

protmuts <- meta[which(meta$WXS %in% mutsamps),]$Proteome_Sample_ID
protwts <- meta[which(meta$WXS %in% wtsamps),]$Proteome_Sample_ID
normals <- meta$Proteome_Normal_Sample_ID[meta$Proteome_Normal_Sample_ID != ""]

#STEP4 (PLOT normal, PLOTmut, PLOTwt x Protein)
prott <- prot
GN <- prot$external_gene_name
prott$external_gene_name = NULL
prottt <- data.frame(t(prott))
colnames(prottt) <- GN
prottt$PSAMP <- rownames(prottt)
prottt$Status <- ifelse(rownames(prottt) %in% normals,"Normal",NA)
prottt$Status <- ifelse(rownames(prottt) %in% protmuts,"NEST_mut",prottt$Status)
prottt$Status <- ifelse(rownames(prottt) %in% protwts,"NEST_wt",prottt$Status)


plotp <- prottt %>% select(PROT,Status)
colnames(plotp) <- c("GENE","STATUS")

p <- ggplot(plotp,aes(x=STATUS,y=GENE))
p <- p + geom_boxplot()
p <- p + ggtitle(paste(NEST,PROT,sep="  "))
p

#Things that I want to plot are there 
#Normal, NESTmut, NESTwt
#Test comparison 

totest1 <- plotp[which(plotp$STATUS != "Normal"),]
totest2 <- totest1[which(!is.na(totest1$STATUS)),]
wilcox.test(GENE~STATUS,data=totest2)

mixed <- c("C3N-01998","C3N-02275","C3N-00498","C3N-03042","C3L-00079")
nogo <- meta[which(meta$CASE_ID %in% mixed),]$Proteome_Normal_Sample_ID

prottt[which(rownames(prottt) %in% nogo),]$TOP1
prottt[which(rownames(prottt) %in% nogo),]$BOD1L1
prottt[which(rownames(prottt) %in% nogo),]$ARHGAP31
