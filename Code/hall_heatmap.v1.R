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

#args <- c("Processed_data/CPTAC.PDAC.hallmark.scores.txt","Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","PDAC","KRAS TP53 CDKN2A")


#args <- c("Processed_data/CPTAC.BRCA.hallmark.scores.txt","Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","BRCA","PIK3CA TP53 PTEN","Processed_data/BRCA.PolyRisk.amp.v2.txt","Processed_data/BRCA.PolyRisk.del.v2.txt","Processed_data/BRCA.PolyRisk.dna.v2.txt","Data/Hallmarks/nanostring.gl.txt")

#args <- c("Processed_data/CPTAC.LUAD.hallmark.scores.txt","Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","LUAD","EGFR TP53 KRAS","Processed_data/LUAD.PolyRisk.amp.v2.txt","Processed_data/LUAD.PolyRisk.del.v2.txt","Processed_data/LUAD.PolyRisk.dna.v2.txt","Data/Hallmarks/nanostring.gl.txt")

#args <- c("Processed_data/CPTAC.HNSCC.hallmark.scores.txt","Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","HNSCC","EGFR TP53 KRAS","Processed_data/HNSCC.PolyRisk.amp.v2.txt","Processed_data/HNSCC.PolyRisk.del.v2.txt","Processed_data/HNSCC.PolyRisk.dna.v2.txt","Data/Hallmarks/nanostring.gl.txt")

#args <- c("Processed_data/CPTAC.GBM.hallmark.scores.txt","Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","GBM","EGFR TP53 KRAS","Processed_data/GBM.PolyRisk.amp.v2.txt","Processed_data/GBM.PolyRisk.del.v2.txt","Processed_data/GBM.PolyRisk.dna.v2.txt","Data/Hallmarks/nanostring.gl.txt")

#Bring in the data
hallscore <- fread(args[1])
maf <- fread(args[2])
meta <- fread(args[3])
cancer = args[4]


#Subset to cancer (NOTE: put in a line for pancan) 
myids = meta[which(meta$cohort == cancer),]
mutid = myids$WXS
cnvid = myids$CASE_ID


#Subset the maf to get keep genes in each cancer type
keep_muts <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins")
maf_cnt <- maf[which(maf$Variant_Classification %in% keep_muts),]
submaf <- maf_cnt[which(maf_cnt$Tumor_Sample_Barcode %in% mutid),]
topgenes <- unlist(str_split(args[5]," "))

#Get Hypermutators 
qrange = quantile(table(submaf$Tumor_Sample_Barcode))
IQR = qrange[4]-qrange[2]
hyperthresh = IQR*1.5+qrange[4]
hypers = names(table(submaf$Tumor_Sample_Barcode)[table(submaf$Tumor_Sample_Barcode) > hyperthresh ])

#Get Tumor Sample Barcodes of mutational hits
maf1 <- unique(submaf[which(submaf$Hugo_Symbol == topgenes[1]),]$Tumor_Sample_Barcode)
maf2 <- unique(submaf[which(submaf$Hugo_Symbol == topgenes[2]),]$Tumor_Sample_Barcode)
maf3 <- unique(submaf[which(submaf$Hugo_Symbol == topgenes[3]),]$Tumor_Sample_Barcode)

#Convert tumore sample barcodes to protein ones
maf1p = meta[which(meta$WXS %in% maf1),]$Proteome_Sample_ID
maf2p = meta[which(meta$WXS %in% maf2),]$Proteome_Sample_ID
maf3p = meta[which(meta$WXS %in% maf3),]$Proteome_Sample_ID
hyperp = meta[which(meta$WXS %in% hypers),]$Proteome_Sample_ID


#now for the NMF clusters and see how they pan out. 
nmf1s = meta[which(meta$Multiomic_subtype_new_annotation == "NMF1/Estrogen_immune_cold"),]$Proteome_Sample_ID
nmf2s = meta[which(meta$Multiomic_subtype_new_annotation == "NMF2/Immune_hot"),]$Proteome_Sample_ID
nmf3s = meta[which(meta$Multiomic_subtype_new_annotation == "NMF3/Proliferation"),]$Proteome_Sample_ID
nmf4s = meta[which(meta$Multiomic_subtype_new_annotation == "NMF4/EMT"),]$Proteome_Sample_ID

#Now I want to pull out some information on the BRCA subtypes to see what falls out. 
pam50 <- fread("Data/Clinical_data/prosp-brca-v5.3-sample-annotation.csv")
Basal <- pam50[which(pam50$PAM50 == "Basal"),]$Sample.ID
Her2 <- pam50[which(pam50$PAM50 == "Her2"),]$Sample.ID
LumA <- pam50[which(pam50$PAM50 == "LumA"),]$Sample.ID
LumB <- pam50[which(pam50$PAM50 == "LumB"),]$Sample.ID
Normal <- pam50[which(pam50$PAM50 == "Normal-like"),]$Sample.ID

#Samps with KRAS (example)
hallscore$HIT1 = ifelse(hallscore$ID %in% maf1p,topgenes[1],"WT")
hallscore$HIT2 = ifelse(hallscore$ID %in% maf2p,topgenes[2],"WT")
hallscore$HIT3 = ifelse(hallscore$ID %in% maf3p,topgenes[3],"WT")
hallscore$NMFS = ifelse(hallscore$ID %in% nmf1s,"NMF1",NA)
hallscore$NMFS = ifelse(hallscore$ID %in% nmf2s,"NMF2",hallscore$NMFS)
hallscore$NMFS = ifelse(hallscore$ID %in% nmf3s,"NMF3",hallscore$NMFS)
hallscore$NMFS = ifelse(hallscore$ID %in% nmf4s,"NMF4",hallscore$NMFS)
hallscore$HYPE = ifelse(hallscore$ID %in% hyperp, "Hypermut","NonHyper")
hallscore$HIST = ifelse(hallscore$ID %in% Basal,"Basal",NA) 
hallscore$HIST = ifelse(hallscore$ID %in% Her2,"Her2",hallscore$HIST) 
hallscore$HIST = ifelse(hallscore$ID %in% LumA,"LumA",hallscore$HIST) 
hallscore$HIST = ifelse(hallscore$ID %in% LumB,"LumB",hallscore$HIST) 
hallscore$HIST = ifelse(hallscore$ID %in% Normal,"Normal",hallscore$HIST) 
 
#Let's try a heatmap now 
h = dcast(hallscore,sectors~HIT1,fun.aggregate=mean, na.rm = TRUE,value.var="value")
h = dcast(hallscore,sectors~HIT1+HIT2+HIT3+HYPE,fun.aggregate=mean, na.rm = TRUE,value.var="value")
h = dcast(hallscore,sectors~NMFS+HYPE,fun.aggregate=mean, na.rm = TRUE,value.var="value")
h = dcast(hallscore,sectors~HIST+HIT1+HIT2+HIT3,fun.aggregate=mean, na.rm = TRUE,value.var="value")

#Lets add some counts to this 
h = dcast(hallscore,sectors~HIST+HYPE,fun.aggregate=mean, na.rm = TRUE,value.var="value")
hc = dcast(hallscore,sectors~HIST+HYPE,fun=length)[,-1]
#hc = dcast(hallscore,sectors~HIT1,fun=length)[,-1]
hcnames <- paste(names(hc),hc[1,],sep="_n=")

hm <- as.matrix(h,rownames="sectors")
colnames(hm) <- hcnames
zhm <- scale((hm))
Heatmap(hm)
Heatmap(zhm)


tophall <-  data.frame(hallscore %>%
  arrange(desc(value)) %>%
  group_by(ID) %>%
  slice(1:3))

##### I just had another idea. It may be interesting: 
#DNA Repair is the top hit for Breast Cancer (duh!) BRCA
#What genes are in the DNA repair and being jacked up. 
amp <- fread(args[6])
del <- fread(args[7])
dna <- fread(args[8])
hallgenes <- fread(args[9])

#NOTE: somehow I need to find the top hits 
#BELOW IS FOR BRCA AFTER THAT IS LUAD with similar code 
th = "Humoral_immunity"
tg = hallgenes %>% select("Genes",any_of(th))
tgs = tg[which(tg[,2] == "Y"),]$Genes

amp_t = amp %>% select("Proteome_Sample_ID",any_of(tgs))
del_t = del %>% select("Proteome_Sample_ID",any_of(tgs))
dna_t = dna %>% select("Proteome_Sample_ID",any_of(tgs))

amp_t$HIST = ifelse(amp_t$Proteome_Sample_ID%in% Basal,"Basal",NA)
amp_t$HIST = ifelse(amp_t$Proteome_Sample_ID%in% Her2,"Her2",amp_t$HIST)
amp_t$HIST = ifelse(amp_t$Proteome_Sample_ID%in% LumA,"LumA",amp_t$HIST)
amp_t$HIST = ifelse(amp_t$Proteome_Sample_ID%in% LumB,"LumB",amp_t$HIST)
amp_t$HIST = ifelse(amp_t$Proteome_Sample_ID%in% Normal,"Normal",amp_t$HIST)
mat <- as.matrix(amp_t[,-1],rownames="HIST")

m_amp <- melt(amp_t)
dat = dcast(m_amp,HIST~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")
mdat <- as.matrix(dat,rownames="HIST")
zmdat <- scale((mdat))
Heatmap(mdat)
Heatmap(zmdat)


del_t$HIST = ifelse(del_t$Proteome_Sample_ID%in% Basal,"Basal",NA)
del_t$HIST = ifelse(del_t$Proteome_Sample_ID%in% Her2,"Her2",del_t$HIST)
del_t$HIST = ifelse(del_t$Proteome_Sample_ID%in% LumA,"LumA",del_t$HIST)
del_t$HIST = ifelse(del_t$Proteome_Sample_ID%in% LumB,"LumB",del_t$HIST)
del_t$HIST = ifelse(del_t$Proteome_Sample_ID%in% Normal,"Normal",del_t$HIST)
mdt <- as.matrix(del_t[,-1],rownames="HIST")

m_del <- melt(del_t)
ddt = dcast(m_del,HIST~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")
mddt <- as.matrix(ddt,rownames="HIST")
zmddt <- scale((mddt))
Heatmap(mddt)
Heatmap(zmddt)



dna_t$HIST = ifelse(dna_t$Proteome_Sample_ID%in% Basal,"Basal",NA)
dna_t$HIST = ifelse(dna_t$Proteome_Sample_ID%in% Her2,"Her2",dna_t$HIST)
dna_t$HIST = ifelse(dna_t$Proteome_Sample_ID%in% LumA,"LumA",dna_t$HIST)
dna_t$HIST = ifelse(dna_t$Proteome_Sample_ID%in% LumB,"LumB",dna_t$HIST)
dna_t$HIST = ifelse(dna_t$Proteome_Sample_ID%in% Normal,"Normal",dna_t$HIST)
mmt <- as.matrix(dna_t[,-1],rownames="HIST")

m_dna <- melt(dna_t)
mmt = dcast(m_dna,HIST~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")
mmmt <- as.matrix(mmt,rownames="HIST")
zmmmt <- scale((mmmt))
Heatmap(mmmt)
Heatmap(zmmmt)


p <- ggplot(m_amp,aes(x=variable,y=HIST,fill=abs(value)))
p <- p + geom_tile()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
p


#Let's look at all of the samples 
sat = dcast(m_amp,Proteome_Sample_ID~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")
sdt = dcast(m_del,Proteome_Sample_ID~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")
smt = dcast(m_dna,Proteome_Sample_ID~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")
 
msat <- as.matrix(sat,rownames="Proteome_Sample_ID")
msdt <- as.matrix(sdt,rownames="Proteome_Sample_ID")
msmt <- as.matrix(smt,rownames="Proteome_Sample_ID")


samps <- unique(data.frame(m_dna %>% select(Proteome_Sample_ID,HIST)))
Heatmap(msmt,row_split=samps$HIST)
zmsmt <- scale((msmt))
Heatmap(zmsmt,row_split=samps$HIST)



comb <- msat+msdt+msmt
Heatmap(comb,row_split=samps$HIST)
zcomb <- scale((comb))


#NOW I WANT TO ADD SOME ANNOTATIONS 
psmaf <- merge(submaf,meta,by.x="Tumor_Sample_Barcode",by.y="WXS",all.x=T)
brcaanno <- c("PIK3CA","TP53","MAP3K1","CEP250","KMT2C","CDH1")
bannomaf <- psmaf[which(psmaf$Hugo_Symbol %in% brcaanno),]
ba_dna <- dcast(bannomaf,Proteome_Sample_ID~Hugo_Symbol,fun=length)
zcomb_ids <- data.frame("Proteome_Sample_ID"=as.vector(rownames(zcomb)))
ba <- merge(zcomb_ids,ba_dna,all.x=T)
ba[is.na(ba)] <- 0 

baha = rowAnnotation(TP53=ba$TP53,PIK3CA=ba$PIK3CA,KMT2C=ba$KMT2C,CDH1=ba$CDH1,MAP3K1=ba$MAP3K1)
#NOW I NEED TO ADD SOME CNV stuffs
cnv_drivers <- fread("Data/Somatic_cnv/CPTAC_cna_driver_calls.txt")
cnv_full <- fread("Data/Somatic_cnv/Broad_pipeline_wxs_formatted/GISTIC_all_thresholded.by_genes.tsv")



pdf("Figures/DNA_damage_BRCA_genenormal_sample.pdf",height=10,width=10,useDingbats=F)
Heatmap(zcomb,row_split=samps$HIST,right_annotation = baha,show_row_names = FALSE)
dev.off()


#TRY to gather the LumA groups 
look <- Heatmap(zcomb,row_split=samps$HIST)

########################################LUAD##########################################33
##### I just had another idea. It may be interesting: 
#DNA Repair is the top hit for Breast Cancer (duh!) BRCA
#What genes are in the DNA repair and being jacked up. 
amp <- fread(args[6])
del <- fread(args[7])
dna <- fread(args[8])
hallgenes <- fread(args[9])

#NOTE: somehow I need to find the top hits 
#BELOW IS FOR BRCA AFTER THAT IS LUAD with similar code 
th = "DNA_repair"
tg = hallgenes %>% select("Genes",any_of(th))
tgs = tg[which(tg[,2] == "Y"),]$Genes

amp_t = amp %>% select("Proteome_Sample_ID",any_of(tgs))
del_t = del %>% select("Proteome_Sample_ID",any_of(tgs))
dna_t = dna %>% select("Proteome_Sample_ID",any_of(tgs))

#SET UP THE COLUMNS 
yo <- hallscore
yo$HIST = paste(yo$HIT1, yo$HIT2, yo$HIT3, yo$HYPE, sep="_")
yo$HIST = paste(yo$NMFS, yo$HYPE,sep="_")
yo_small <- unique(yo %>% select("ID","HIST"))

amp_t <- merge(amp_t,yo_small,by.x="Proteome_Sample_ID",by.y="ID")
mat <- as.matrix(amp_t[,-1],rownames="HIST")
m_amp <- melt(amp_t)
dat = dcast(m_amp,HIST~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")
mdat <- as.matrix(dat,rownames="HIST")
zmdat <- scale((mdat))
Heatmap(mdat)
Heatmap(zmdat)


del_t <- merge(del_t,yo_small,by.x="Proteome_Sample_ID",by.y="ID")
mdt <- as.matrix(del_t[,-1],rownames="HIST")
m_del <- melt(del_t)
ddt = dcast(m_del,HIST~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")
mddt <- as.matrix(ddt,rownames="HIST")
zmddt <- scale((mddt))
Heatmap(mddt)
Heatmap(zmddt)


dna_t <- merge(dna_t,yo_small,by.x="Proteome_Sample_ID",by.y="ID")
mmt <- as.matrix(dna_t[,-1],rownames="HIST")
m_dna <- melt(dna_t)
mmt = dcast(m_dna,HIST~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")
mmmt <- as.matrix(mmt,rownames="HIST")
zmmmt <- scale((mmmt))
Heatmap(mmmt)
Heatmap(zmmmt)



#Let's look at all of the samples 
sat = dcast(m_amp,Proteome_Sample_ID~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")
sdt = dcast(m_del,Proteome_Sample_ID~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")
smt = dcast(m_dna,Proteome_Sample_ID~variable,fun.aggregate=mean, na.rm = TRUE,value.var="value")

sat$E2F5 <- NULL
smt$E2F5 <- NULL

msat <- as.matrix(sat,rownames="Proteome_Sample_ID")
msdt <- as.matrix(sdt,rownames="Proteome_Sample_ID")
msmt <- as.matrix(smt,rownames="Proteome_Sample_ID")


samps <- unique(data.frame(m_dna %>% select(Proteome_Sample_ID,HIST)))
Heatmap(msmt,row_split=samps$HIST)
zmsmt <- scale((msmt))
Heatmap(zmsmt,row_split=samps$HIST)


comb <- msat+msdt+msmt
Heatmap(comb,row_split=samps$HIST)
zcomb <- scale((comb))
pdf("Figures/RAS_pathway_LUAD_genenormal_sample.pdf",height=10,width=10,useDingbats=F)
Heatmap(zcomb,row_split=samps$HIST,)
dev.off()


###### NOW I want to pull out some information on why the LumA patients look so different with DNA damage genes... generally from the genomics contributions. ###### 
luma <- fread("Data/plop",header=F)
lumaTSB <- meta[which(meta$Proteome_Sample_ID %in% luma$V1),]$WXS #this grabs the tumor samples barcodes
lumaMAF <- submaf[which(submaf$Tumor_Sample_Barcode %in% lumaTSB),]
ml <- merge(luma,meta,by.x="V1", by.y="Proteome_Sample_ID")
lumaMAFmeta = merge(lumaMAF,ml,by.x="Tumor_Sample_Barcode", by.y="WXS",all.x=T)
lumaDNA <- data.frame(lumaMAFmeta %>% group_by(V2,Hugo_Symbol) %>% tally() %>% ungroup() %>% arrange(desc(n)))
head(lumaDNA,35)

#GENE to annotate with
psmaf <- merge(submaf,meta,by.x="Tumor_Sample_Barcode",by.y="WXS",all.x=T)
brcaanno <- c("PIK3CA","TP53","MAP3K1","CEP250","KMT2C","CDH1")
bannomaf <- psmaf[which(psmaf$Hugo_Symbol %in% brcaanno),]
ba_dna <- dcast(bannomaf,Proteome_Sample_ID~Hugo_Symbol,fun=length)
zcomb_ids <- data.frame("Proteome_Sample_ID"=as.vector(rownames(zcomb)))
ba <- merge(zcomb_ids,ba_dna,all.x=T)
ba[is.na(ba)] <- 0 


#Annotate with genes
mcnv <- melt(cnv_full,id.vars=c("Gene Symbol","Gene ID","Cytoband"))
mluma <- mcnv[which(mcnv$variable %in% ml$CASE_ID),]
lmc <- merge(mluma,ml,by.x="variable",by.y="CASE_ID",all.x=T)
lmc0 <- lmc[which(lmc$value !=0 ),]
nest <- c("GOLGA2","RAB1A","RAB1B","RUSC2")
nl <- lmc0[which(lmc0$"Gene Symbol" %in% nest),]
table(nl %>% select(V2,`Gene Symbol`,value))
lumaCNV <- data.frame(lmc0 %>% group_by(V2,`Gene Symbol`) %>% tally() %>% ungroup() %>% arrange(desc(n)))
lmc_tops <- lmc %>% group_by(V2,`Gene Symbol`) %>% 


dcnv <- dcast(mcnv,`variable`~`Gene Symbol`,value.var="value")
cnv_meta <- merge(dcnv,meta,by.x="variable",by.y="CASE_ID",all.x=T)
lumaCASEID <- meta[which(meta$Proteome_Sample_ID %in% luma$V1),]$CASE_ID
lumaCNV <- dcnv[which(dcnv$variable %in% lumaCASEID),]
lmc <- merge(dcnv,ml,by.x="variable",by.y="CASE_ID",all.x=T)

lmc_tops <- lmc %>% group_by(V2,`Gene Symbol`) %>% 





############################################## LET's Go through an examples of HNSCC ########################### 
#Subset to cancer (NOTE: put in a line for pancan) 
hallscore <- fread(args[1])
myids = meta[which(meta$cohort == cancer),]
mutid = myids$WXS
cnvid = myids$CASE_ID


#Subset the maf to get keep genes in each cancer type
keep_muts <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins")
maf_cnt <- maf[which(maf$Variant_Classification %in% keep_muts),]
submaf <- maf_cnt[which(maf_cnt$Tumor_Sample_Barcode %in% mutid),]
topgenes <- unlist(str_split(args[5]," "))


tophall <-  data.frame(hallscore %>%
  arrange(desc(value)) %>%
  group_by(ID) %>%
  slice(1:3))








##############################HERE IS SOME GBM WORK 



args <- c("Processed_data/CPTAC.GBM.hallmark.scores.txt","Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Clustering/CPTAC-Pancan-Data_metatable.txt","GBM","EGFR TP53 KRAS","Processed_data/GBM.PolyRisk.amp.v2.txt","Processed_data/GBM.PolyRisk.del.v2.txt","Processed_data/GBM.PolyRisk.dna.v2.txt","Data/Hallmarks/nanostring.gl.txt")
hallscore <- fread(args[1])
maf <- fread(args[2])
meta <- fread(args[3])
cancer = args[4]



#Subset to cancer (NOTE: put in a line for pancan) 
myids = meta[which(meta$cohort == cancer),]
mutid = myids$WXS
cnvid = myids$CASE_ID
#Subset the maf to get keep genes in each cancer type
keep_muts <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins")
maf_cnt <- maf[which(maf$Variant_Classification %in% keep_muts),]
submaf <- maf_cnt[which(maf_cnt$Tumor_Sample_Barcode %in% mutid),]
topgenes <- unlist(str_split(args[5]," "))
amp <- fread(args[6])
del <- fread(args[7])
dna <- fread(args[8])
hallgenes <- fread(args[9])
unique(hallscore$sectors)
th = "Humoral_immunity"
tg = hallgenes %>% select("Genes",any_of(th))
tgs = tg[which(tg[,2] == "Y"),]$Genes
tgs
amp_t = amp %>% select("Proteome_Sample_ID",any_of(tgs))
del_t = del %>% select("Proteome_Sample_ID",any_of(tgs))
dna_t = dna %>% select("Proteome_Sample_ID",any_of(tgs))

mamp <- as.matrix(amp_t[,-1],rownames = amp_t$Proteome_Sample_ID)
mdel <- as.matrix(del_t[,-1],rownames = del_t$Proteome_Sample_ID)
mdna <- as.matrix(dna_t[,-1],rownames = dna_t$Proteome_Sample_ID)

comb <- abs(mamp)+abs(mdel)+abs(mdna) 
Heatmap(comb)
zcomb <- scale(t(comb))
Heatmap(zcomb)

