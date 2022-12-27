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

'%!in%' <- function(x,y)!('%in%'(x,y))
args = commandArgs(trailingOnly=TRUE)


#args <- c("Data/ssGSEA/Pancan_driver_DNA_ssGSEA-combined.gct","Data/ssGSEA/Pancan_driver_RNA_ssGSEA-combined.gct","Data/ssGSEA/Pancan_driver_protein_ssGSEA-combined.gct","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_raw_gene_level_v1.3.tsv","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt")


dna <- fread(args[1])
rna <- fread(args[2])
protein <- fread(args[3])
meta <- fread(args[4])

#can and id 
meta_lim <- select(meta,c(Proteome_Sample_ID,cohort,Sample,CASE_ID))
meta_lim$variable <- sub("^(\\d)", "X\\1",meta_lim$CASE_ID)
meta_lim$variable2 <- str_replace(meta_lim$variable,"-","\\.")


#So here is where it is again
rownames(dna) <- dna$id
dna$id <- NULL
toremove <- c("Signature","pvalue","fdr","No.columns.scored",".A")

dna_c <- select(dna,!contains(toremove))
rna_c <- select(rna,!contains(toremove))
protein_c <- select(protein,!contains(toremove))





#######################################DNA#######################
#Now I need to melt, merge, split, plot
dna_c$Hallmarks = rownames(dna_c)
dna_m <- melt(dna_c)
dna_mc <- merge(dna_m,meta_lim,by.x ="variable", by.y="Proteome_Sample_ID",all.x=T)

#Now I need to grab the cancertypes 
cancers = unique(meta_lim$cohort)

for(i in cancers){
    new_dat <- dna_mc[which(dna_mc$cohort == i),]
    dd <- dcast(new_dat,Hallmarks~variable,value.var="value")
    rownames(dd) <- dd$Hallmarks
    halls <- dd$Hallmarks
    dd$Hallmarks = NULL
    mm <- as.matrix(dd)
    rownames(mm) <- halls 
    Heatmap(mm)
    smm <- as.matrix(scale(mm))
        ofname <- paste("Processed_data/",i,".ssGSEA.sample.matrix.DNA.txt",sep="")
    figname <- paste("Figures/",i,".ssGSEA.heatmap.DNA.pdf",sep="")
    write.table(smm,ofname,row.names=T, sep="\t",quote=F)
    
    pdf(figname,height=10,width=18,useDingbats=F)
    print(Heatmap(smm))
    dev.off()

}


#####################################RNA############################
#Now I need to melt, merge, split, plot
rna_m <- melt(rna_c)
#These names seem quite odd. I'm going to convert to the . id and go from there: here is my attempt at this: 
rna_m$variable2 <- str_replace(rna_m$variable,"\\.T","") 


rna_mc <- merge(rna_m,meta_lim,by="variable2",all.x=T)

#Now I need to grab the cancertypes 
cancers = unique(meta_lim$cohort)

for(i in cancers){
    new_dat <- rna_mc[which(rna_mc$cohort == i),]
    dd <- dcast(new_dat,id~Proteome_Sample_ID,value.var="value")
    halls <- dd$id
    rownames(dd) <- dd$id
    dd$id = NULL
    mm <- as.matrix(dd)
    rownames(mm) <- halls
    mss <- as.matrix(scale(mm))
    ofname <- paste("Processed_data/",i,".ssGSEA.sample.matrix.RNA.txt",sep="")
    figname <- paste("Figures/",i,".ssGSEA.heatmap.RNA.pdf",sep="")
    write.table(mss,ofname,row.names=T, sep="\t",quote=F)
    
    pdf(figname,height=10,width=18,useDingbats=F)
    print(Heatmap(mss))
    dev.off()
}




###########################PROTEIN##################################
#Now I need to melt, merge, split, plot
protein_m <- melt(protein_c)
#These names seem quite odd. I'm going to convert to the . id and go from there: here is my attempt at this: 


protein_mc <- merge(protein_m,meta_lim,by.x="variable",by.y="Proteome_Sample_ID")

#Now I need to grab the cancertypes 
cancers = unique(meta_lim$cohort)

for(i in cancers){
    new_dat <- protein_mc[which(protein_mc$cohort == i),]
    dd <- dcast(new_dat,id~variable,value.var="value")
    halls <- dd$id
    rownames(dd) <- dd$id
    dd$id = NULL
    mm <- as.matrix(dd)
    rownames(mm) <- halls
    Heatmap(mm)
    mms <- as.matrix(scale(mm))

        ofname <- paste("Processed_data/",i,".ssGSEA.sample.matrix.Protein.txt",sep="")
    figname <- paste("Figures/",i,".ssGSEA.heatmap.Protein.pdf",sep="")
    write.table(mms,ofname,row.names=T, sep="\t",quote=F)
    
    pdf(figname,height=10,width=18,useDingbats=F)
    print(Heatmap(mms))
    dev.off()
}


##### NOW LET ME ADD THE SMOKING DATA ######## 
smoking = fread("Data/Meta_table/count.Neoantigen.Mutation.AllCancerType.v2.smoking.tsv")

meta_lim <- select(meta,c("CASE_ID","Sample","Proteome_Sample_ID","cohort"))
smokingm <- merge(smoking,meta_lim,by.x="Sample",by.y="CASE_ID",all.x=T)
cansmoke <- smokingm[which(smokingm$Proteome_Sample_ID %in% unique(protein_mc$variable)),]
os <- cansmoke[order(cansmoke$Proteome_Sample_ID),]
#Now I think that I want to do all of this by cancertype noow 

cans <- unique(os$cohort)

for(i in cans){
    os_can = os[which(os$cohort == i),]
    halls = protein_mc[which(protein_mc$cohort == i),]
    #set up the annotations 
    baha = columnAnnotation(Count=os_can$Count,Mutation=os_can$Mutation,Smoking=os_can$Smoking)
    #set up the matrix
    hs <- dcast(halls,variable~id,value.var="value")
    hs2 <- hs[which(hs$variable %in% os_can$Proteome_Sample_ID),]
    rownames(hs2) = hs2$variable
    hs2$variable <- NULL
    mh <- as.matrix(hs2,rownames=rownames(hs2))
    #Heatmap(mh)
    zmh <- scale(t(mh))
    
    ofname = paste("Figures/C3PO.",i,".Smoking.ssGSEA.heatmap.pdf",sep="")
    pdf(ofname,height=7,width=7,useDingbats=F)
    print(Heatmap(zmh,bottom_annotation=baha,show_row_names = TRUE,show_column_names=FALSE,name=i))
    dev.off()
}






#Now I'm going to look into the protein data specifically for LUAD and UCEC and the EMT pathway that is realy different in a few of these samples. 

ucec <- fread("Processed_data/UCEC.ssGSEA.sample.matrix.Protein.txt")
luad <- fread("Processed_data/LUAD.ssGSEA.sample.matrix.Protein.txt")

mu <- melt(ucec)
ml <- melt(luad)


emtu <- mu[which(mu$V1 == "Epithelial-mesenchymal_transition"),]
emtl <- ml[which(ml$V1 == "Epithelial-mesenchymal_transition"),]

plot(density(emtu$value))
plot(density(emtl$value))


highu <- emtu[which(emtu$value > 0),]$variable
highl <- emtl[which(emtl$value > 0),]$variable

#NOW GET THE GENE LISTS AND PROTEINS NEEDED
halls = fread("Data/Hallmarks/nanostring.gl.txt")
myhall = c("Genes","Epithelial-mesenchymal_transition")

myh = halls %>% select(all_of(myhall))
genesInHall = myh[which(myh[,2] == "Y"),]$Genes

#Get Proteins in the list
prot <- fread("Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_raw_gene_level_v1.3.tsv")
protgenes <- prot$external_gene_name

#Intersect hallmarks with Proteins
all_genes <- intersect(protgenes,genesInHall)

#all of the UCEC and LUAD sample ids 
ucec_samp = colnames(ucec[,-1])
luad_samp = colnames(luad[,-1])


ph_ucec = prot %>% select("external_gene_name",all_of(ucec_samp))
ph_luad = prot %>% select("external_gene_name",all_of(luad_samp))

phg_ucec <- ph_ucec[which(ph_ucec$external_gene_name %in% all_genes),]
phg_luad <- ph_luad[which(ph_luad$external_gene_name %in% all_genes),]

mp_ucec <- na.omit(melt(phg_ucec))
mp_luad <- na.omit(melt(phg_luad))

dp_ucec <- na.omit(dcast(mp_ucec,external_gene_name~variable))
dp_luad <- na.omit(dcast(mp_luad,external_gene_name~variable))

mp_ucec$EMT <- ifelse(mp_ucec$variable %in% highu,"High","Low")
mp_luad$EMT <- ifelse(mp_luad$variable %in% highl,"High","Low")


#There are too many to look by eye, so I need to look at it systematically. 
emtgenes <- unique(mp_ucec$external_gene_name)
COMPU = NULL
i="AKT3"
for(i in emtgenes){
    dat <- mp_ucec[which(mp_ucec$external_gene_name == i),]
    a <- dat[which(dat$EMT == "High"),]$value
    b <- dat[which(dat$EMT == "Low"),]$value
    if(length(a) > 5 & length(b) > 5){
        myt <- t.test(a,b)
        out = data.frame("Gene"=i,"pvalue"=myt$p.value,CIlow=myt$conf.int[1],CIhigh=myt$conf.int[2])
        COMPU = rbind(COMPU,out)
    }
}
COMPU[order(COMPU$pvalue),]

emtgenes <- unique(mp_luad$external_gene_name)
COMPL = NULL
i="AKT3"
for(i in emtgenes){
    dat <- mp_luad[which(mp_luad$external_gene_name == i),]
    a <- dat[which(dat$EMT == "High"),]$value
    b <- dat[which(dat$EMT == "Low"),]$value
    if(length(a) > 5 & length(b) > 5){
        myt <- t.test(a,b)
        out = data.frame("Gene"=i,"pvalue"=myt$p.value,CIlow=myt$conf.int[1],CIhigh=myt$conf.int[2])
        COMPL = rbind(COMPL,out)
    }
}
COMPL[order(COMPL$pvalue),]


mykeep <- c("VPS13A","FSTL1","TMEM30B")
p <- ggplot(mp_ucec[which(mp_ucec$external_gene_name %in% mykeep),],aes(external_gene_name,value,fill=factor(EMT)))
p <- p + geom_boxplot()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
p


p <- ggplot(mp_luad,aes(external_gene_name,value,fill=factor(EMT)))
p <- p + geom_boxplot()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
p












