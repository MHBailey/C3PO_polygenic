library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)


'%!in%' <- function(x,y)!('%in%'(x,y))
args = commandArgs(trailingOnly=TRUE)

#This is just going to be some scratch code for one examples before 
#args = c("Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_imputed_gene_level.tsv","Processed_data/LUAD.PolyRisk.dna.v2.txt","Processed_data/LUAD.PolyRisk.amp.v2.txt","Processed_data/LUAD.PolyRisk.del.v2.txt","LUAD",'Processed_data/LUAD.SampleCorrelations.dna.v1.txt','Processed_data/LUAD.SampleCorrelations.amp.v1.txt','Processed_data/LUAD.SampleCorrelations.del.v1.txt','Processed_data/LUAD.SampleCorrelations.combined.v1.txt')

#args = c("Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_imputed_gene_level.tsv", "Processed_data/COAD.PolyRisk.dna.v2.txt", "Processed_data/COAD.PolyRisk.amp.v2.txt", "Processed_data/COAD.PolyRisk.del.v2.txt", "COAD", "Processed_data/COAD.SampleCorrelations.dna.v1.txt", "Processed_data/COAD.SampleCorrelations.amp.v1.txt", "Processed_data/COAD.SampleCorrelations.del.v1.txt", "Processed_data/COAD.SampleCorrelations.combined.v1.txt")

prot = fread(args[1])
prs_dna = data.frame(fread(args[2]))
rownames(prs_dna) <- prs_dna$Proteome_Sample_ID 
prs_dna$Proteome_Sample_ID <- NULL 
prs_amp = data.frame(fread(args[3]))
rownames(prs_amp) <- prs_amp$Proteome_Sample_ID
prs_amp$Proteome_Sample_ID <- NULL 
prs_del = data.frame(fread(args[4]))
rownames(prs_del) <- prs_del$Proteome_Sample_ID
prs_del$Proteome_Sample_ID <- NULL
cancer = args[5]

#Clean up PROT 
GN <- prot$external_gene_name
prot$external_gene_name = NULL
prott <- data.frame(t(prot))
colnames(prott) <- GN

prottt <- prott[which(rownames(prott) %in% rownames(prs_dna)),]

genes = colnames(prott)


CORRS_dna = NULL
g = "EGFR"
for(g in genes){
    if(g %in% colnames(prs_dna) & g %in% colnames(prott)){
        tmp_mut = prs_dna %>% select(mutg=all_of(g))
        tmp_prot = prottt %>% select(protg=all_of(g))
        pm = merge(tmp_mut,tmp_prot,by=0)
        rownames(pm) <- pm$"Row.names"
        pm$"Row.names" <- NULL
        corrs = cor.test(pm[,1],pm[,2])
        #plot(pm[,1],pm[,2])
        out = data.frame("Protein" = g, "Coef"=corrs$estimate ,"P.value"=corrs$p.value,"ConfHigh"=corrs$conf.int[2],"ConfLow"=corrs$conf.int[1],"Substrate"="DNA")
        CORRS_dna = rbind(CORRS_dna,out)
    }
}
write.table(CORRS_dna,args[6],sep="\t",quote=F,row.names=F)
#head(CORRS_dna[order(CORRS_dna$P.value),],75)
#tail(CORRS_dna[order(CORRS_dna$P.value),],75)


CORRS_amp = NULL
g = "EGFR"
for(g in genes){
    if(g %in% colnames(prs_amp) & g %in% colnames(prott)){
        tmp_mut = prs_amp %>% select(mutg=all_of(g))
        tmp_prot = prottt %>% select(protg=all_of(g))
        pm = merge(tmp_mut,tmp_prot,by=0)
        rownames(pm) <- pm$"Row.names"
        pm$"Row.names" <- NULL
        corrs = cor.test(pm[,1],pm[,2]) 
        plot(pm[,1],pm[,2])
        out = data.frame("Protein" = g, "Coef"=corrs$estimate ,"P.value"=corrs$p.value,"ConfHigh"=corrs$conf.int[2],"ConfLow"=corrs$conf.int[1],"Substrate"="AMP")
        CORRS_amp = rbind(CORRS_amp,out)
    }
}
write.table(CORRS_amp,args[7],sep="\t",quote=F,row.names=F)
#head(CORRS_amp[order(CORRS_amp$P.value),],75)
#tail(CORRS_amp[order(CORRS_amp$P.value),],75)

CORRS_del = NULL
g = "EGFR"
for(g in genes){
    if(g %in% colnames(prs_del) & g %in% colnames(prott)){
        tmp_mut = prs_del %>% select(mutg=all_of(g))
        tmp_prot = prottt %>% select(protg=all_of(g))
        pm = merge(tmp_mut,tmp_prot,by=0)
        rownames(pm) <- pm$"Row.names"
        pm$"Row.names" <- NULL
        corrs = cor.test(pm[,1],pm[,2])
        plot(pm[,1],pm[,2])
        out = data.frame("Protein" = g, "Coef"=corrs$estimate ,"P.value"=corrs$p.value,"ConfHigh"=corrs$conf.int[2],"ConfLow"=corrs$conf.int[1],"Substrate"="AMP")
        CORRS_del = rbind(CORRS_del,out)
    }
}
write.table(CORRS_del,args[8],sep="\t",quote=F,row.names=F)
#head(CORRS_del[order(CORRS_del$P.value),],75)
#tail(CORRS_del[order(CORRS_del$P.value),],75)




CORRS_combined = NULL
g = "EGFR"
g = "GGA2"
for(g in genes){
    if(g %in% colnames(prs_dna) &  g %in% colnames(prs_amp) & g %in% colnames(prs_del) & g %in% colnames(prott)){
        tmp_mut = prs_dna %>% select(mut_dna=all_of(g))
        tmp_amp = prs_amp %>% select(mut_amp=all_of(g))
        tmp_del = prs_del %>% select(mut_amp=all_of(g))
        tmp_prot = prottt %>% select(protg=all_of(g))
        ma = merge(tmp_mut,tmp_amp,by=0)
        rownames(ma) = ma$"Row.names"
        ma$"Row.names" <- NULL
        mad = merge(ma,tmp_del,by=0)
        rownames(mad) <- mad$"Row.names"
        mad$"Row.names" <- NULL
        madp = merge(mad,tmp_prot,by=0)
        rownames(madp) <- madp$"Row.names"
        madp$"Row.names" <- NULL

        mylm = lm(madp[,4] ~ madp[,1]+madp[,2]+madp[,3])
        myg = glance(mylm)
        #plot(madp[,1]+madp[,2]+madp[,3],madp[,4],main=g)

        out = data.frame("Protein" = g, "Coef"=myg$adj.r.squared,"P.value"=myg$p.value,"Substrate"="COMBINED")
        CORRS_combined = rbind(CORRS_combined,out)
    }
}
write.table(CORRS_combined,args[9],sep="\t",quote=F,row.names=F)
#head(CORRS_combined[order(CORRS_combined$P.value),],75)
#tail(CORRS_combined[order(CORRS_combined$P.value),],75)






###############PHOSPHO 

