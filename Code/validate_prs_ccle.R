library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)
library(stringr)

'%!in%' <- function(x,y)!('%in%'(x,y))



args = commandArgs(trailingOnly=TRUE)

#This is just going to be some scratch code for one examples before 
#args = c("Validation/CCLE.OV.updated.protein.table.tsv", "Validation/CCLE.OV.cis.PolyRisk.dna.v2.txt", "Validation/CCLE.OV.cis.PolyRisk.amp.v2.txt", "Validation/CCLE.OV.cis.PolyRisk.del.v2.txt", "OV", "Validation/CCLE.OV.pancan.PolyRisk.dna.v2.txt", "Validation/CCLE.OV.pancan.PolyRisk.amp.v2.txt", "Validation/CCLE.OV.pancan.PolyRisk.del.v2.txt", "Validation/CCLE.OV.SampleCorrelations.dna.v1.txt", "Validation/CCLE.OV.SampleCorrelations.amp.v1.txt", "Validation/CCLE.OV.SampleCorrelations.del.v1.txt", "Validation/CCLE.OV.SampleCorrelations.combined.v1.txt", "Validation/CCLE.OV.pancan.SampleCorrelations.dna.v1.txt", "Validation/CCLE.OV.pancan.SampleCorrelations.amp.v1.txt", "Validation/CCLE.OV.pancan.SampleCorrelations.del.v1.txt", "Validation/CCLE.OV.pancan.SampleCorrelations.combined.v1.txt", "Figures/CCLE.validation.OV.dna.qqplot.pdf", "Figures/CCLE.validation.OV.cnva.qqplot.pdf", "Figures/CCLE.validation.OV.cnvd.qqplot.pdf", "Figures/CCLE.validation.OV.combined.qqplot.pdf", "Figures/CCLE.validation.OV.pdna.qqplot.pdf", "Figures/CCLE.validation.OV.pcnva.qqplot.pdf", "Figures/CCLE.validation.OV.pcnvd.qqplot.pdf", "Figures/CCLE.validation.OV.pcombined.qqplot.pdf")

#args = c("Validation/CCLE.BRCA.updated.protein.table.tsv", "Validation/CCLE.BRCA.cis.PolyRisk.dna.v2.txt", "Validation/CCLE.BRCA.cis.PolyRisk.amp.v2.txt", "Validation/CCLE.BRCA.cis.PolyRisk.del.v2.txt", "BRCA", "Validation/CCLE.BRCA.pancan.PolyRisk.dna.v2.txt", "Validation/CCLE.BRCA.pancan.PolyRisk.amp.v2.txt", "Validation/CCLE.BRCA.pancan.PolyRisk.del.v2.txt", "Validation/CCLE.BRCA.SampleCorrelations.dna.v1.txt", "Validation/CCLE.BRCA.SampleCorrelations.amp.v1.txt", "Validation/CCLE.BRCA.SampleCorrelations.del.v1.txt", "Validation/CCLE.BRCA.SampleCorrelations.combined.v1.txt", "Validation/CCLE.BRCA.pancan.SampleCorrelations.dna.v1.txt", "Validation/CCLE.BRCA.pancan.SampleCorrelations.amp.v1.txt", "Validation/CCLE.BRCA.pancan.SampleCorrelations.del.v1.txt", "Validation/CCLE.BRCA.pancan.SampleCorrelations.combined.v1.txt", "Figures/CCLE.validation.BRCA.dna.qqplot.pdf", "Figures/CCLE.validation.BRCA.cnva.qqplot.pdf", "Figures/CCLE.validation.BRCA.cnvd.qqplot.pdf", "Figures/CCLE.validation.BRCA.combined.qqplot.pdf", "Figures/CCLE.validation.BRCA.pdna.qqplot.pdf", "Figures/CCLE.validation.BRCA.pcnva.qqplot.pdf", "Figures/CCLE.validation.BRCA.pcnvd.qqplot.pdf", "Figures/CCLE.validation.BRCA.pcombined.qqplot.pdf")

#args = c("Validation/CCLE.COAD.updated.protein.table.tsv", "Validation/CCLE.COAD.cis.PolyRisk.dna.v2.txt", "Validation/CCLE.COAD.cis.PolyRisk.amp.v2.txt", "Validation/CCLE.COAD.cis.PolyRisk.del.v2.txt", "COAD", "Validation/CCLE.COAD.pancan.PolyRisk.dna.v2.txt", "Validation/CCLE.COAD.pancan.PolyRisk.amp.v2.txt", "Validation/CCLE.COAD.pancan.PolyRisk.del.v2.txt", "Validation/CCLE.COAD.SampleCorrelations.dna.v1.txt", "Validation/CCLE.COAD.SampleCorrelations.amp.v1.txt", "Validation/CCLE.COAD.SampleCorrelations.del.v1.txt", "Validation/CCLE.COAD.SampleCorrelations.combined.v1.txt", "Validation/CCLE.COAD.pancan.SampleCorrelations.dna.v1.txt", "Validation/CCLE.COAD.pancan.SampleCorrelations.amp.v1.txt", "Validation/CCLE.COAD.pancan.SampleCorrelations.del.v1.txt", "Validation/CCLE.COAD.pancan.SampleCorrelations.combined.v1.txt", "Figures/CCLE.validation.COAD.dna.qqplot.pdf", "Figures/CCLE.validation.COAD.cnva.qqplot.pdf", "Figures/CCLE.validation.COAD.cnvd.qqplot.pdf", "Figures/CCLE.validation.COAD.combined.qqplot.pdf", "Figures/CCLE.validation.COAD.pdna.qqplot.pdf", "Figures/CCLE.validation.COAD.pcnva.qqplot.pdf", "Figures/CCLE.validation.COAD.pcnvd.qqplot.pdf", "Figures/CCLE.validation.COAD.pcombined.qqplot.pdf")

#SOME NECESSARY FUNCTIONS: #https://slowkow.com/notes/ggplot2-qqplot/ 

inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  return(lambda)
}

gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  lambda = inflation(ps)
  mylabel = sprintf("λ = %.2f", lambda)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(length(-log10(sort(ps))))),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:length(-log10(sort(ps))), shape2 = length(-log10(sort(ps))):1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:length(-log10(sort(ps))), shape2 = length(-log10(sort(ps))):1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po) +
    annotate("text", x = 1, y = 0, label = mylabel,size = 8)
}


prot = fread(args[1])
prs_dna = data.frame(fread(args[2]))
rownames(prs_dna) <- prs_dna$DepMap_ID
prs_dna$Proteome_Sample_ID <- NULL
prs_amp = data.frame(fread(args[3]))
rownames(prs_amp) <- prs_amp$DepMap_ID
prs_amp$Proteome_Sample_ID <- NULL
prs_del = data.frame(fread(args[4]))
rownames(prs_del) <- prs_del$DepMap_ID
prs_del$Proteome_Sample_ID <- NULL
cancer = args[5]

prott <- data.frame(prot)
SN <- prot$Sample
prott$Sample = NULL
colnames(prott) <- make.unique(colnames(prott), sep="_")
prottt <- prott %>% select(where(~mean(is.na(.))< 0.5))
genes = colnames(prottt)
rownames(prottt) = make.unique(SN)


CORRS_dna = NULL
g = "XIAP"
for(g in genes){
    g2 = str_split_fixed(g,"_",2)[,1]
    if(g2 %in% colnames(prs_dna) & g %in% colnames(prottt)){
        tmp_mut = prs_dna %>% select(mutg=all_of(g2))
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
write.table(CORRS_dna,args[9],sep="\t",quote=F,row.names=F)
#head(CORRS_dna[order(CORRS_dna$P.value),],75)
#tail(CORRS_dna[order(CORRS_dna$P.value),],75)



CORRS_amp = NULL
g = "EGFR"
for(g in genes){
    g2 = str_split_fixed(g,"_",2)[,1]
    if(g2 %in% colnames(prs_amp) & g %in% colnames(prott)){
        tmp_mut = prs_amp %>% select(mutg=all_of(g2))
        tmp_prot = prottt %>% select(protg=all_of(g))
        pm = merge(tmp_mut,tmp_prot,by=0)
        rownames(pm) <- pm$"Row.names"
        pm$"Row.names" <- NULL
        corrs = cor.test(pm[,1],pm[,2])
        #plot(pm[,1],pm[,2])
        out = data.frame("Protein" = g, "Coef"=corrs$estimate ,"P.value"=corrs$p.value,"ConfHigh"=corrs$conf.int[2],"ConfLow"=corrs$conf.int[1],"Substrate"="AMP")
        CORRS_amp = rbind(CORRS_amp,out)
    }
}
write.table(CORRS_amp,args[10],sep="\t",quote=F,row.names=F)
#head(CORRS_amp[order(CORRS_amp$P.value),],75)
#tail(CORRS_amp[order(CORRS_amp$P.value),],75)

CORRS_del = NULL
g = "EGFR"
for(g in genes){
    g2 = str_split_fixed(g,"_",2)[,1]
    if(g2 %in% colnames(prs_del) & g %in% colnames(prott)){
        tmp_mut = prs_del %>% select(mutg=all_of(g2))
        tmp_prot = prottt %>% select(protg=all_of(g))
        pm = merge(tmp_mut,tmp_prot,by=0)
        rownames(pm) <- pm$"Row.names"
        pm$"Row.names" <- NULL
        corrs = cor.test(pm[,1],pm[,2])
        #plot(pm[,1],pm[,2])
        out = data.frame("Protein" = g, "Coef"=corrs$estimate ,"P.value"=corrs$p.value,"ConfHigh"=corrs$conf.int[2],"ConfLow"=corrs$conf.int[1],"Substrate"="DEL")
        CORRS_del = rbind(CORRS_del,out)
    }
}
write.table(CORRS_del,args[11],sep="\t",quote=F,row.names=F)
#head(CORRS_del[order(CORRS_del$P.value),],75)
#tail(CORRS_del[order(CORRS_del$P.value),],75)



CORRS_combined = NULL
g = "EGFR"
g = "GGA2"
g = "S100A3"
for(g in genes){
    g2 = str_split_fixed(g,"_",2)[,1]
    if(g2 %in% colnames(prs_dna) &  g2 %in% colnames(prs_amp) & g2 %in% colnames(prs_del) & g %in% colnames(prott)){
        tmp_mut = prs_dna %>% select(mut_dna=all_of(g2))
        tmp_amp = prs_amp %>% select(mut_amp=all_of(g2))
        tmp_del = prs_del %>% select(mut_amp=all_of(g2))
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
write.table(CORRS_combined,args[12],sep="\t",quote=F,row.names=F)
#head(CORRS_combined[order(CORRS_combined$P.value),],75)
#tail(CORRS_combined[order(CORRS_combined$P.value),],75)


pdf(args[17],height=5,width=5,useDingbats=F)
gg_qqplot(CORRS_dna$P.value) +
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )
dev.off()

pdf(args[18],height=5,width=5,useDingbats=F)
gg_qqplot(CORRS_amp$P.value) +
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )
dev.off()

pdf(args[19],height=5,width=5,useDingbats=F)
gg_qqplot(CORRS_del$P.value) +
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )
dev.off()

pdf(args[20],height=5,width=5,useDingbats=F)
gg_qqplot(CORRS_combined$P.value) +
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )
dev.off()




###########################################PANCAN##############################################
prs_dna = data.frame(fread(args[6]))
rownames(prs_dna) <- prs_dna$DepMap_ID
prs_dna$Proteome_Sample_ID <- NULL
prs_amp = data.frame(fread(args[7]))
rownames(prs_amp) <- prs_amp$DepMap_ID
prs_amp$Proteome_Sample_ID <- NULL
prs_del = data.frame(fread(args[8]))
rownames(prs_del) <- prs_del$DepMap_ID
prs_del$Proteome_Sample_ID <- NULL
cancer = args[5]

CORRS_dna = NULL
g = "DNAJC17"
for(g in genes){
    g2 = str_split_fixed(g,"_",2)[,1]
    if(g2 %in% colnames(prs_dna) & g %in% colnames(prott)){
        tmp_mut = prs_dna %>% select(mutg=all_of(g2))
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
write.table(CORRS_dna,args[13],sep="\t",quote=F,row.names=F)
#head(CORRS_dna[order(CORRS_dna$P.value),],75)
#tail(CORRS_dna[order(CORRS_dna$P.value),],75)


CORRS_amp = NULL
g = "EGFR"
for(g in genes){
    g2 = str_split_fixed(g,"_",2)[,1]
    if(g2 %in% colnames(prs_amp) & g %in% colnames(prott)){
        tmp_mut = prs_amp %>% select(mutg=all_of(g2))
        tmp_prot = prottt %>% select(protg=all_of(g))
        pm = merge(tmp_mut,tmp_prot,by=0)
        rownames(pm) <- pm$"Row.names"
        pm$"Row.names" <- NULL
        corrs = cor.test(pm[,1],pm[,2])
        #plot(pm[,1],pm[,2])
        out = data.frame("Protein" = g, "Coef"=corrs$estimate ,"P.value"=corrs$p.value,"ConfHigh"=corrs$conf.int[2],"ConfLow"=corrs$conf.int[1],"Substrate"="AMP")
        CORRS_amp = rbind(CORRS_amp,out)
    }
}
write.table(CORRS_amp,args[14],sep="\t",quote=F,row.names=F)
#head(CORRS_amp[order(CORRS_amp$P.value),],75)
#tail(CORRS_amp[order(CORRS_amp$P.value),],75)

CORRS_del = NULL
g = "EGFR"
for(g in genes){
    g2 = str_split_fixed(g,"_",2)[,1]
    if(g2 %in% colnames(prs_del) & g %in% colnames(prott)){
        tmp_mut = prs_del %>% select(mutg=all_of(g2))
        tmp_prot = prottt %>% select(protg=all_of(g))
        pm = merge(tmp_mut,tmp_prot,by=0)
        rownames(pm) <- pm$"Row.names"
        pm$"Row.names" <- NULL
        corrs = cor.test(pm[,1],pm[,2])
        out = data.frame("Protein" = g, "Coef"=corrs$estimate ,"P.value"=corrs$p.value,"ConfHigh"=corrs$conf.int[2],"ConfLow"=corrs$conf.int[1],"Substrate"="DEL")
        CORRS_del = rbind(CORRS_del,out)
    }
}
write.table(CORRS_del,args[15],sep="\t",quote=F,row.names=F)
#head(CORRS_del[order(CORRS_del$P.value),],75)
#tail(CORRS_del[order(CORRS_del$P.value),],75)


CORRS_combined = NULL
g = "EGFR"
g = "GGA2"
g = "S100A3"
for(g in genes){
    g2 = str_split_fixed(g,"_",2)[,1]
    if(g2 %in% colnames(prs_dna) &  g2 %in% colnames(prs_amp) & g2 %in% colnames(prs_del) & g %in% colnames(prott)){
        tmp_mut = prs_dna %>% select(mut_dna=all_of(g2))
        tmp_amp = prs_amp %>% select(mut_amp=all_of(g2))
        tmp_del = prs_del %>% select(mut_amp=all_of(g2))
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
write.table(CORRS_combined,args[16],sep="\t",quote=F,row.names=F)
#head(CORRS_combined[order(CORRS_combined$P.value),],75)
#tail(CORRS_combined[order(CORRS_combined$P.value),],75)


pdf(args[21],height=5,width=5,useDingbats=F)
gg_qqplot(CORRS_dna$P.value) +
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )
dev.off()

pdf(args[22],height=5,width=5,useDingbats=F)
gg_qqplot(CORRS_amp$P.value) +
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )
dev.off()

pdf(args[23],height=5,width=5,useDingbats=F)
gg_qqplot(CORRS_del$P.value) +
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )
dev.off()

pdf(args[24],height=5,width=5,useDingbats=F)
gg_qqplot(CORRS_combined$P.value) +
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )
dev.off()

       

