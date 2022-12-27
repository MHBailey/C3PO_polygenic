library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)
library(stringr)
library(circlize)


'%!in%' <- function(x,y)!('%in%'(x,y))

args = commandArgs(trailingOnly=TRUE)

#args = c("Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_imputed_gene_level.tsv","BRCA",'Data/Hallmarks/nanostring.gl.txt','Figures/ProteinOnly.CPTAC.BRCA.hallmark.boxplots.pdf','Processed_data/ProteinOnly.CPTAC.BRCA.hallmark.scores.txt','Figures/ProteinOnly.Circos.CPTAC.BRCA.hallmarks.pdf',"Data/Meta_table/CPTAC-Pancan-Data_metatable.txt")

#args = c("Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_imputed_gene_level.tsv","PDAC",'Data/Hallmarks/nanostring.gl.txt','Figures/ProteinOnly.CPTAC.PDAC.hallmark.boxplots.pdf','Processed_data/ProteinOnly.CPTAC.PDAC.hallmark.scores.txt','Figures/ProteinOnly.Circos.CPTAC.PDAC.hallmarks.pdf',"Data/Meta_table/CPTAC-Pancan-Data_metatable.txt")

#args = c("Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_imputed_gene_level.tsv","{can}",'Data/Hallmarks/nanostring.gl.txt','Figures/ProteinOnly.CPTAC.{can}.hallmark.boxplots.pdf','Processed_data/ProteinOnly.CPTAC.{can}.hallmark.scores.txt','Figures/ProteinOnly.Circos.CPTAC.{can}.hallmarks.pdf',"Data/Meta_table/CPTAC-Pancan-Data_metatable.txt")

prot = fread(args[1])
cancer = args[2]
halls = data.frame(fread(args[3]))
meta = fread(args[7])
metafav = data.frame(meta %>% select("Proteome_Sample_ID","cohort"))
if(cancer != "PANCAN"){
    canmeta <- metafav[which(metafav$cohort == cancer),]$Proteome_Sample_ID
}else{
    canmeta <- metafav$Proteome_Sample_ID
}

#Clean up PROT 
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



p <- ggplot(mhall,aes(x=Var2,y=value))
p <- p + geom_boxplot()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))


pdf(args[4],height=5,width=8,useDingbats=F)
print(p)
dev.off()


#Not I just need to find the top 3 pathways and see if they are matching across these groups

mhall$value2 = abs(mhall$value)
tophall <-  data.frame(tbl_df(mhall) %>%
  group_by(Var1) %>%
  arrange(value2, .by_group = TRUE) %>%
  top_n(3))


colnames(mhall) <- c("ID","sectors","value","absvalue")

write.table(mhall,args[5],sep="\t",quote=F,row.names=F)

#NOW ONTO THE CIRCOS

#Figure out the pictures 
library(yaml)
data = yaml.load_file("hallmarks.v2.yaml")
hallmark_list = data$symbol
hallmark_name = sapply(hallmark_list, function(x) x$name)
hallmark_src = sapply(hallmark_list, function(x) x$src)

library(EBImage)
pdf("Figures/Hallmarks.thumbnails.circos.pdf",heigh=10,width=10,useDingbats=F)
circos.par("points.overflow.warning" = FALSE)
circos.initialize(hallmark_name, xlim = c(0, 1))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
    image = EBImage::readImage(hallmark_src[CELL_META$sector.numeric.index])
    rasterImage(image,
        xleft = pos[1, 1] - 0.225, ybottom = pos[1, 2] - 0.15,
        xright = pos[1, 1] + 0.225, ytop = pos[1, 2]+ 0.15)
}, bg.border = 0, track.height = .25)
dev.off()

#Now I want to put it together: 
library(yaml)
data = yaml.load_file("hallmarks.v2.yaml")
hallmark_list = data$symbol
hallmark_name = sapply(hallmark_list, function(x) x$name)
hallmark_src = sapply(hallmark_list, function(x) x$src)



#Set up the link colors# 
colors = fread("colors.txt")

#ADD CT to tophall 
meta = fread(args[7])
metafav = data.frame(meta %>% select("Proteome_Sample_ID","cohort"))
tophall$Var3 <- rep(c("T1","T2","T3"))
tophall_merged <- merge(tophall,metafav,by.x="Var1",by.y="Proteome_Sample_ID")
types = unique(tophall_merged$cohort)


pdf(args[6],heigh=10,width=10,useDingbats=F)
circos.par("track.height" = 0.3)
circos.initialize(mhall$sectors, xlim=c(min(mhall$value),max(mhall$value)))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
    image = EBImage::readImage(hallmark_src[CELL_META$sector.numeric.index])
    rasterImage(image,
        xleft = pos[1, 1] - 0.08, ybottom = pos[1, 2] - 0.08,
        xright = pos[1, 1] + 0.08, ytop = pos[1, 2]+ 0.08)
}, bg.border = 1, track.height = .25)
circos.trackHist(mhall$sectors, x = mhall$value,"track.height" = 0.15)

#ADD CT to tophall 
for(t in types){
    hey <- dcast(tophall_merged[which(tophall_merged$cohort == t),],Var1~Var3,value.var="Var2")
    mycomb <- combn(unique(tophall_merged[which(tophall_merged$cohort == t),]$Var2),2)


    for(i in 1:dim(mycomb)[2]){ 
        cancol = colors[which(colors$Cancer == t),]$Color 
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
        ids1 = tophall_merged[which(tophall_merged$cohort == t),][which(tophall_merged[which(tophall_merged$cohort == t),]$Var1 %in% hey8$Var1 & tophall_merged[which(tophall_merged$cohort == t),]$Var2 == p1),]$value
        ids2 = tophall_merged[which(tophall_merged$cohort == t),][which(tophall_merged[which(tophall_merged$cohort == t),]$Var1 %in% hey8$Var1 & tophall_merged[which(tophall_merged$cohort == t),]$Var2 == p2),]$value
        if(dim(hey8)[1] != 0){
            circos.link(p1,ids1,p2,ids2,h.ratio=.4,col=cancol)
        }
    }
}
legend("topleft", pch = 1, legend = cancer)
dev.off()




