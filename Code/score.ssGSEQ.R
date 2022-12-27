library(data.table)
library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)


'%!in%' <- function(x,y)!('%in%'(x,y))
args = commandArgs(trailingOnly=TRUE)

#args = c("Processed_data/{can}.ssGSEA.sample.matrix.DNA.txt","Processed_data/{can}.ssGSEA.sample.matrix.RNA.txt","Processed_data/{can}.ssGSEA.sample.matrix.Protein.txt","{can}")

#args = c("Processed_data/HNSCC.ssGSEA.sample.matrix.DNA.txt","Processed_data/HNSCC.ssGSEA.sample.matrix.RNA.txt","Processed_data/HNSCC.ssGSEA.sample.matrix.Protein.txt","HNSCC")

#args = c("Processed_data/UCEC.ssGSEA.sample.matrix.DNA.txt","Processed_data/UCEC.ssGSEA.sample.matrix.RNA.txt","Processed_data/UCEC.ssGSEA.sample.matrix.Protein.txt","UCEC")

#args = c("Processed_data/BRCA.ssGSEA.sample.matrix.DNA.txt","Processed_data/BRCA.ssGSEA.sample.matrix.RNA.txt","Processed_data/BRCA.ssGSEA.sample.matrix.Protein.txt","BRCA")

#args = c("Processed_data/OV.ssGSEA.sample.matrix.DNA.txt","Processed_data/OV.ssGSEA.sample.matrix.RNA.txt","Processed_data/OV.ssGSEA.sample.matrix.Protein.txt","OV")


dna <- fread(args[1])
rna <- fread(args[2])
protein <- fread(args[3])
can = args[4]

#Rank by Sample
dr <- data.frame(apply(dna[,-1], 2, rank))
rownames(dr) <- dna$V1
rr <- data.frame(apply(rna[,-1], 2, rank))
rownames(rr) <- rna$V1
pr <- data.frame(apply(protein[,-1], 2, rank))
rownames(pr) <- protein$V1

tmp <- intersect(colnames(dr),colnames(rr))
allsamples <- intersect(tmp,colnames(pr))

drs <- select(dr,all_of(allsamples))
rrs <- select(rr,all_of(allsamples))
prs <- select(pr,all_of(allsamples))


COMPARISONS = NULL
for(i in 1:dim(rrs)[2]){
    #Do the drs v rrs
    name = "DNAvRNA"
    sample = colnames(drs)[i]
    wtd <- wilcox.test(drs[,i],rrs[,i],paired=T)
    out=data.frame(Comparison = name, sample=sample, pval=wtd$p.value,stat="SignedRank",cancer=can)
    COMPARISONS = rbind(COMPARISONS,out)

    #Do the drs v prs
    name = "DNAvProtein"
    wtd <- wilcox.test(drs[,i],prs[,i],paired=T)
    out=data.frame(Comparison = name,sample=sample, pval=wtd$p.value,stat="SignedRank",cancer=can)
    COMPARISONS = rbind(COMPARISONS,out)

    #Do the rrs v prs 
    name = "RNAvProtein"
    wtd <- wilcox.test(rrs[,i],prs[,i],paired=T)
    out=data.frame(Comparison = name,sample=sample, pval=wtd$p.value,stat="SignedRank",cancer=can)
    COMPARISONS = rbind(COMPARISONS,out)
}

#head(COMPARISONS[order(COMPARISONS$pval),],100)

write.table(COMPARISONS,args[5],sep="\t",quote=F,row.names=F)

gtitle = paste(can,"Sample Ranks",sep=" ")
p <- ggplot(COMPARISONS,aes(x=pval,fill=Comparison))
p <- p + geom_density(alpha=.4)
p <- p + ggtitle(gtitle)

pdf(args[6],height=5,width=5,useDingbats=F)
print(p)
dev.off()















######################CALC RANK BY HALLMARK####################
#Rank by Hallmark
dr <- data.frame(t(apply(dna[,-1], 1, rank)))
rownames(dr) <- dna$V1
rr <- data.frame(t(apply(rna[,-1], 1, rank)))
rownames(rr) <- rna$V1
pr <- data.frame(t(apply(protein[,-1], 1, rank)))
rownames(pr) <- protein$V1

tmp <- intersect(colnames(dr),colnames(rr))
allsamples <- intersect(tmp,colnames(pr))

drs <- select(dr,all_of(allsamples))
rrs <- select(rr,all_of(allsamples))
prs <- select(pr,all_of(allsamples))


COMPARISONS = NULL
for(i in 1:dim(rrs)[2]){
    #Do the drs v rrs
    name = "DNAvRNA"
    sample = colnames(drs)[i]
    wtd <- wilcox.test(drs[,i],rrs[,i],paired=T)
    out=data.frame(Comparison = name, sample=sample, pval=wtd$p.value,stat="SignedRank",cancer=can)
    COMPARISONS = rbind(COMPARISONS,out)

    #Do the drs v prs
    name = "DNAvProtein"
    wtd <- wilcox.test(drs[,i],prs[,i],paired=T)
    out=data.frame(Comparison = name,sample=sample, pval=wtd$p.value,stat="SignedRank",cancer=can)
    COMPARISONS = rbind(COMPARISONS,out)

    #Do the rrs v prs 
    name = "RNAvProtein"
    wtd <- wilcox.test(rrs[,i],prs[,i],paired=T)
    out=data.frame(Comparison = name,sample=sample, pval=wtd$p.value,stat="SignedRank",cancer=can)
    COMPARISONS = rbind(COMPARISONS,out)
}

#head(COMPARISONS[order(COMPARISONS$pval),],100)

write.table(COMPARISONS,args[7],sep="\t",quote=F,row.names=F)

gtitle = paste(can,"Hallmark Ranks",sep=" ")
p <- ggplot(COMPARISONS,aes(x=pval,fill=Comparison))
p <- p + geom_density(alpha=.4)
p <- p + ggtitle(gtitle)

pdf(args[8],height=5,width=5,useDingbats=F)
print(p)
dev.off()
