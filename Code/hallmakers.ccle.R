library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)
library(stringr)
library(circlize)


'%!in%' <- function(x,y)!('%in%'(x,y))

args = commandArgs(trailingOnly=TRUE)

#args = c('Validation/BRCA.cis.PolyRisk.dna.v2.txt','Validation/BRCA.cis.PolyRisk.amp.v2.txt','Validation/BRCA.cis.PolyRisk.del.v2.txt',"BRCA",'Validation/BRCA.pancan.PolyRisk.dna.v2.txt','Validation/BRCA.pancan.PolyRisk.amp.v2.txt','Validation/BRCA.pancan.PolyRisk.del.v2.txt','Data/Hallmarks/nanostring.gl.txt')

#args = c('Validation/OV.cis.PolyRisk.dna.v2.txt','Validation/OV.cis.PolyRisk.amp.v2.txt','Validation/OV.cis.PolyRisk.del.v2.txt',"OV",'Validation/OV.pancan.PolyRisk.dna.v2.txt','Validation/OV.pancan.PolyRisk.amp.v2.txt','Validation/OV.pancan.PolyRisk.del.v2.txt','Data/Hallmarks/nanostring.gl.txt')

#args = c('Validation/COAD.cis.PolyRisk.dna.v2.txt','Validation/COAD.cis.PolyRisk.amp.v2.txt','Validation/COAD.cis.PolyRisk.del.v2.txt',"COAD",'Validation/COAD.pancan.PolyRisk.dna.v2.txt','Validation/COAD.pancan.PolyRisk.amp.v2.txt','Validation/COAD.pancan.PolyRisk.del.v2.txt','Data/Hallmarks/nanostring.gl.txt')


prs_dna = data.frame(fread(args[1]))
rownames(prs_dna) <- prs_dna$DepMap_ID
prs_dna$Proteome_Sample_ID <- NULL
prs_amp = data.frame(fread(args[2]))
rownames(prs_amp) <- prs_amp$DepMap_ID
prs_amp$Proteome_Sample_ID <- NULL
prs_del = data.frame(fread(args[3]))
rownames(prs_del) <- prs_del$DepMap_ID
prs_del$Proteome_Sample_ID <- NULL
cancer = args[4]
halls = data.frame(fread(args[5]))

h = colnames(halls)[-1] #Removes genes 

HALLMARKS = NULL
for(i in h){
    hg = halls %>% select("Genes",all_of(i))
    gs = hg[which(hg[,2]=="Y"),]$Genes
    dnas = prs_dna %>% select(any_of(gs))
    amps = prs_amp %>% select(any_of(gs))
    dels = prs_del %>% select(any_of(gs))
    dnas$ID = rownames(dnas)
    amps$ID = rownames(amps)
    dels$ID = rownames(dels)
    
    yo <- data.frame(bind_rows(dnas,amps,dels) %>%
    group_by(ID) %>% 
    summarise_all(sum, na.rm = T))
    rownames(yo) <- yo$ID
    yo$ID = NULL
    yo2 = data.frame(rowSums(yo))
    colnames(yo2) = i
    totg = max(c(dim(dnas)[2],dim(amps)[2],dim(dels)[2]))
    yo3 = yo2/totg
    if(is.null(HALLMARKS)){
        HALLMARKS = yo3
    }else{
        HALLMARKS = cbind(HALLMARKS,yo3)
    }
}

mhall <- melt(as.matrix(HALLMARKS))

p <- ggplot(mhall,aes(x=Var2,y=value))
p <- p + geom_boxplot()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))


pdf(args[6],height=5,width=8,useDingbats=F)
print(p)
dev.off()

##############################PANCAN###########################################33
prs_dna = data.frame(fread(args[7]))
rownames(prs_dna) <- prs_dna$DepMap_ID
prs_dna$Proteome_Sample_ID <- NULL
prs_amp = data.frame(fread(args[8]))
rownames(prs_amp) <- prs_amp$DepMap_ID
prs_amp$Proteome_Sample_ID <- NULL
prs_del = data.frame(fread(args[9]))
rownames(prs_del) <- prs_del$DepMap_ID
prs_del$Proteome_Sample_ID <- NULL

h = colnames(halls)[-1] #Removes genes 

HALLMARKS = NULL
for(i in h){
    hg = halls %>% select("Genes",all_of(i))
    gs = hg[which(hg[,2]=="Y"),]$Genes
    dnas = prs_dna %>% select(any_of(gs))
    amps = prs_amp %>% select(any_of(gs))
    dels = prs_del %>% select(any_of(gs))
    dnas$ID = rownames(dnas)
    amps$ID = rownames(amps)
    dels$ID = rownames(dels)

    yo <- data.frame(bind_rows(dnas,amps,dels) %>%
    group_by(ID) %>%
    summarise_all(sum, na.rm = T))
    rownames(yo) <- yo$ID
    yo$ID = NULL
    yo2 = data.frame(rowSums(yo))
    colnames(yo2) = i
    totg = max(c(dim(dnas)[2],dim(amps)[2],dim(dels)[2]))
    yo3 = yo2/totg
    if(is.null(HALLMARKS)){
        HALLMARKS = yo3
    }else{
        HALLMARKS = cbind(HALLMARKS,yo3)
    }
}

mhall <- melt(as.matrix(HALLMARKS))

p <- ggplot(mhall,aes(x=Var2,y=value))
p <- p + geom_boxplot()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))


pdf(args[10],height=5,width=8,useDingbats=F)
print(p)
dev.off()

#Not I just need to find the top 3 pathways and see if they are matching across these groups
