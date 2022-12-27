library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)
library(stringr)
library(circlize)


'%!in%' <- function(x,y)!('%in%'(x,y))

args = commandArgs(trailingOnly=TRUE)

#args = c('Processed_data/BRCA.PolyRisk.dna.v2.txt','Processed_data/BRCA.PolyRisk.amp.v2.txt','Processed_data/BRCA.PolyRisk.del.v2.txt',"BRCA",'Data/Hallmarks/nanostring.gl.txt','Figures/CPTAC.BRCA.hallmark.boxplots.pdf','Processed_data/CPTAC.BRCA.hallmark.scores.txt','Figures/Circos.CPTAC.BRCA.hallmarks.pdf',"Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","./hallmarks.v2.yaml","Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_raw_gene_level_v1.3.tsv")

#args = c('Processed_data/PDAC.PolyRisk.dna.v2.txt','Processed_data/PDAC.PolyRisk.amp.v2.txt','Processed_data/PDAC.PolyRisk.del.v2.txt',"PDAC",'Data/Hallmarks/nanostring.gl.txt','Figures/CPTAC.PDAC.hallmark.boxplots.pdf','Processed_data/CPTAC.PDAC.hallmark.scores.txt','Figures/Circos.CPTAC.PDAC.hallmarks.pdf',"Data/Meta_table/CPTAC-Pancan-Data_metatable.txt")

#args = c('Processed_data/{can}.PolyRisk.dna.v2.txt','Processed_data/{can}.PolyRisk.amp.v2.txt','Processed_data/{can}.PolyRisk.del.v2.txt',"{can}",'Data/Hallmarks/nanostring.gl.txt','Figures/CPTAC.{can}.hallmark.boxplots.pdf','Processed_data/CPTAC.{can}.hallmark.scores.txt','Figures/Circos.CPTAC.{can}.hallmarks.pdf',"Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","./hallmarks.v2.yaml","Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_imputed_gene_level.tsv")

prs_dna = data.frame(fread(args[1]))
rownames(prs_dna) <- prs_dna$Proteome_Sample_ID
prs_dna$Proteome_Sample_ID <- NULL
prs_amp = data.frame(fread(args[2]))
rownames(prs_amp) <- prs_amp$Proteome_Sample_ID
prs_amp$Proteome_Sample_ID <- NULL
prs_del = data.frame(fread(args[3]))
rownames(prs_del) <- prs_del$Proteome_Sample_ID
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
    yo2 = data.frame(rowSums(abs(yo))) #THIS IS THE CHANGE TO figure out direction
    #yo2 = data.frame(rowSums((yo))) #THIS IS THE CHANGE TO figure out direction
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

#Not I just need to find the top 3 pathways and see if they are matching across these groups

mhall$value2 = abs(mhall$value)

tophall <-  data.frame(mhall %>%
  arrange(desc(value2)) %>%
  group_by(Var1) %>%
  slice(1:3))


colnames(mhall) <- c("ID","sectors","value","absvalue")


#NOW ONTO THE CIRCOS
#Now I want to put it together: 
library(yaml)
data = yaml.load_file(args[10])
hallmark_list = data$symbol
hallmark_name = sapply(hallmark_list, function(x) x$name)
hallmark_src = sapply(hallmark_list, function(x) x$src)



#Set up the link colors# 
colors = fread("colors.txt")

#ADD CT to tophall 
meta = fread(args[9])
metafav = data.frame(meta %>% select("Proteome_Sample_ID","cohort"))
tophall$Var3 <- rep(c("T1","T2","T3"))
tophall_merged <- merge(tophall,metafav,by.x="Var1",by.y="Proteome_Sample_ID")
types = unique(tophall_merged$cohort)


##############NOW CALCULATE THE PROTEIN SCORES ############################
prot = fread(args[11])
metafav = data.frame(meta %>% select("Proteome_Sample_ID","cohort"))
if(cancer != "PANCAN"){
    canmeta <- metafav[which(metafav$cohort == cancer),]$Proteome_Sample_ID
}else{
    canmeta <- metafav$Proteome_Sample_ID
}

GN <- prot$external_gene_name
prot$external_gene_name = NULL
prott <- data.frame(t(prot))
colnames(prott) <- GN

prottt <- prott[which(rownames(prott) %in% canmeta),]

genes = colnames(prottt)

h = colnames(halls)[-1] #Removes genes 

HALLMARKS = NULL
for(i in h){
    hg = halls %>% select("Genes",all_of(i))
    gs = hg[which(hg[,2]=="Y"),]$Genes
    proteins = prottt %>% select(any_of(gs))
    proteins$ID = rownames(proteins)
    
    yo <- data.frame(proteins %>%
    group_by(ID) %>%
    summarise_all(sum, na.rm = T))
    rownames(yo) <- yo$ID
    yo$ID = NULL
    yo2 = data.frame(rowSums(abs(yo))) #THIS IS THE CHANGE TO figure out direction
    colnames(yo2) = i
    totg = max(c(dim(proteins)[2]))
    yo3 = yo2/totg
    if(is.null(HALLMARKS)){
        HALLMARKS = yo3
    }else{
        HALLMARKS = cbind(HALLMARKS,yo3)
    }
}

mhall_p <- melt(as.matrix(HALLMARKS))
mhall_p$value2 = abs(mhall_p$value)

tophall_p <-  data.frame(mhall_p %>%
  arrange(desc(value2)) %>%
  group_by(Var1) %>%
  slice(1:3))

colnames(mhall_p) <- c("ID","sectors","value","absvalue")

tophall_p$Var3 <- rep(c("T1","T2","T3"))
tophall_p_merged <- merge(tophall_p,metafav,by.x="Var1",by.y="Proteome_Sample_ID")
types_p = unique(tophall_p_merged$cohort)


tophall_merged$Scaled= scale(tophall_merged$value2)
tophall_p_merged$Scaled= scale(tophall_p_merged$value2)
mhall$Scaled <- scale(mhall$absvalue)
mhall_p$Scaled <- scale(mhall_p$absvalue)

#ADD colorgs
colors = fread("colors.txt")


pdf(args[8],heigh=10,width=10,useDingbats=F)
circos.par("track.height" = 0.3)
circos.initialize(mhall$sectors, xlim=c(min(mhall$Scaled),max(mhall$Scaled)))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
    image = EBImage::readImage(hallmark_src[CELL_META$sector.numeric.index])
    rasterImage(image,
        xleft = pos[1, 1] - 0.08, ybottom = pos[1, 2] - 0.08,
        xright = pos[1, 1] + 0.08, ytop = pos[1, 2]+ 0.08)
}, bg.border = 1, track.height = .25)
circos.trackHist(mhall$sectors, x = mhall$Scaled,"track.height" = 0.15)


#ADD CT to tophall 
for(t in types_p){
    hey_p <- dcast(tophall_p_merged[which(tophall_p_merged$cohort == t),],Var1~Var3,value.var="Var2")
    mycomb_p <- combn(unique(tophall_p_merged[which(tophall_p_merged$cohort == t),]$Var2),2)


    for(i in 1:dim(mycomb_p)[2]){
        cancol = colors[which(colors$Cancer == t),]$Color
        hey_p8 = NULL
        pair = mycomb_p[,i]
        p1 = pair[1]
        p2 = pair[2]
        hey_p2 = hey_p[which(hey_p$T1 == p1 & hey_p$T2 == p2),]
        hey_p3 = hey_p[which(hey_p$T2 == p1 & hey_p$T1 == p2),]

        hey_p4 = hey_p[which(hey_p$T2 == p1 & hey_p$T3 == p2),]
        hey_p5 = hey_p[which(hey_p$T3 == p1 & hey_p$T2 == p2),]

        hey_p6 = hey_p[which(hey_p$T1 == p1 & hey_p$T3 == p2),]
        hey_p7 = hey_p[which(hey_p$T3 == p1 & hey_p$T1 == p2),]

        hey_p8 = rbind(hey_p2,hey_p3,hey_p4,hey_p5,hey_p6,hey_p7)
        ids1 = tophall_p_merged[which(tophall_p_merged$cohort == t),][which(tophall_p_merged[which(tophall_p_merged$cohort == t),]$Var1 %in% hey_p8$Var1 & tophall_p_merged[which(tophall_p_merged$cohort == t),]$Var2 == p1),]$Scaled
        ids2 = tophall_p_merged[which(tophall_p_merged$cohort == t),][which(tophall_p_merged[which(tophall_p_merged$cohort == t),]$Var1 %in% hey_p8$Var1 & tophall_p_merged[which(tophall_p_merged$cohort == t),]$Var2 == p2),]$Scaled
        if(dim(hey_p8)[1] != 0){
            circos.link(p1,ids1,p2,ids2,h.ratio=.4,col=cancol)
        }
    }
}

#ADD CT to tophall 
for(t in types){
    hey <- dcast(tophall_merged[which(tophall_merged$cohort == t),],Var1~Var3,value.var="Var2")
    mycomb <- combn(unique(tophall_merged[which(tophall_merged$cohort == t),]$Var2),2)


    for(i in 1:dim(mycomb)[2]){ 
        cancol = colors[which(colors$Cancer == t),]$ColorLight
        hey8 = NULL
        pair = mycomb[,i]
        p1 = pair[1]
        p2 = pair[2]
        hey2 = hey[which(hey$T1 == p1 & hey$T2 == p2),]
        hey3 = hey[which(hey$T2 == p1 & hey$T1 == p2),]

        hey4 = hey[which(hey$T2 == p1 & hey$T3 == p2),]
        hey5 = hey[which(hey$T3 == p1 & hey$T2 == p2),]

        hey6 = hey[which(hey$T1 == p1 & hey$T3 == p2),]
        hey7 = hey[which(hey$T3 == p1 & hey$T1 == p2),]

        hey8 = rbind(hey2,hey3,hey4,hey5,hey6,hey7)
        ids1 = tophall_merged[which(tophall_merged$cohort == t),][which(tophall_merged[which(tophall_merged$cohort == t),]$Var1 %in% hey8$Var1 & tophall_merged[which(tophall_merged$cohort == t),]$Var2 == p1),]$Scaled
        ids2 = tophall_merged[which(tophall_merged$cohort == t),][which(tophall_merged[which(tophall_merged$cohort == t),]$Var1 %in% hey8$Var1 & tophall_merged[which(tophall_merged$cohort == t),]$Var2 == p2),]$Scaled
        if(dim(hey8)[1] != 0){
            circos.link(p1,ids1,p2,ids2,h.ratio=.4,col="darkgrey")
        }
    }
}

legend("topleft", pch = 1, legend = cancer)
dev.off()



pdf(args[6])
plot(1,1)
dev.off()
pdf(args[7])
plot(1,1)
dev.off()
