library(data.table)
library(reshape2)
library(ggplot2)
library(stringr)
library(dplyr)

'%!in%' <- function(x,y)!('%in%'(x,y))

args = commandArgs(trailingOnly=TRUE)

#Rscript --vanilla --quiet code.R file1.txt file2.txt 

args <- c("Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Gene_expression/ALL_RNA-Seq_Expr_WashU_FPKM_UQ_annotation.tsv","Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_raw_gene_level_v1.3.tsv","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt")


maf=fread(args[1])

maf1 <- maf[,c(1,4,5,6,8,13,25)] #essentials

#Restrict to my deep dive
indeltypes <- c("DE_NOVO_START_IN_FRAME","DE_NOVO_START_OUT_FRAME","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Nonsense_Mutation","Nonstop_Mutation","START_CODON_INS","START_CODON_SNP","Splice_Site")
maf2 <- maf1[which(maf1$Variant_Classification %in% indeltypes ),]
#45006 events 

#I think I want to use splice-sites as a special case. removed in maf3 

maf3 <- maf2[which(maf2$Variant_Classification != "Splice_Site"),]
#34747 events 
#Count events by genes
m3cnt <- data.frame(t(sort(table(maf3$Hugo_Symbol))))
m3cnt <- m3cnt[,-1]

#Count events by samples
maf4 <- data.frame(unique(maf3[,c(1,6)]))
#33726 events
#This means that there were 1021 events with multiple indels in the same sample and gene 
m4cnt <- data.frame(t(sort(table(maf4$Hugo_Symbol))))
m4cnt <- m4cnt[,-1]

m3.m4 <- merge(m3cnt,m4cnt,by="Var2")
m34order <- m3.m4[order(m3.m4$Freq.x, decreasing = T),]

m34top50 <- head(m34order,50)

m34top50$ordGene <- factor(x=m34top50$Var2,levels=m34top50$Var2)

mm34 <- melt(m34top50)
mm34$CountType = ifelse(mm34$variable == "Freq.x","All","Unique")



#Now what 
p <- ggplot(mm34,aes(x=ordGene,y=value,fill=CountType))
p <- p + geom_bar(stat="identity",position="dodge")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
p

#NOW I want to bring in proteomics and RNAseq data to look at the effect of mutations on 

rna = fread(args[2])
prot = fread(args[3])
meta = fread(args[4])

#Clean up RNA names 
newnames <- colnames(rna)
nn <- str_replace_all(newnames,"-A","_N")
nnn <- str_replace_all(nn,"-T","_T")
colnames(rna) <- nnn
tokeep <- intersect(nnn,meta$Sample)

tokeep2 <- c("gene_name",tokeep)
rna2 <- select(rna,tokeep2)
tomatch <- meta[which(meta$Sample %in% tokeep),][,c(2,5)]
setnames(rna2,tomatch$Sample,tomatch$Proteome_Sample_ID)
mr <- melt(rna2)

#Clean up Protein data 
mp <- melt(prot)

#META_WES_PROT 
mwp = meta[,c(14,5)]

#Get the top 50genes to test
top50 <- m34top50$Var2
#g="APC" #I might want to swap out the maf3 for maf And I may want to revers this so that all RNA and Protein is shown and mutation is layed ontop. 

for(g in top50){
    muts <- maf[which(maf$Hugo_Symbol == g),]
    mutsmeta <- merge(muts,mwp,by.x="Tumor_Sample_Barcode",by.y="WXS",all.x=T)
    gr <- mr[which(mr$gene_name==g),]
    gp <- mp[which(mp$external_gene_name ==g),]
    mutr <- merge(mutsmeta,gr,by.x="Proteome_Sample_ID",by.y="variable",all.y=T)
    mutrp <- merge(mutr,gp,by.x="Proteome_Sample_ID",by.y="variable",all.y=T)

    scattername = paste("Figures/Scatter_",g,".20220930.pdf",sep="")
    RNAdensename = paste("Figures/RNAdensity_",g,".20220930.pdf",sep="")
    PROTdensename = paste("Figures/PROTdensity_",g,".20220930.pdf",sep="")  

    pdf(scattername,height=8,width=8)
    p <- ggplot(mutrp,aes(x=value.x,y=value.y,color=Variant_Classification))
    p <- p + geom_point()
    print(p)
    dev.off()
   
    pdf(RNAdensename,height=8,width=8)
    p <- ggplot(mutrp,aes(x=value.x,group=Variant_Classification,fill=Variant_Classification))
    p <- p + geom_boxplot(alpha=.3)
    p <- p + ggtitle("RNAseq breakdown")
    print(p)
    dev.off()

    pdf(PROTdensename,height=8,width=8)
    p <- ggplot(mutrp,aes(x=value.y,group=Variant_Classification,fill=Variant_Classification))
    p <- p + geom_boxplot(alpha=.3)
    p <- p + ggtitle("Preotome breakdown")
    print(p)
    dev.off()

}
#mutationdata maf3 
#rna data = mr
#protein data = mp 




