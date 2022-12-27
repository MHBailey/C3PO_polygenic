library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)
library(stringr)
library(circlize)


'%!in%' <- function(x,y)!('%in%'(x,y))

args = commandArgs(trailingOnly=TRUE)

#args = c('Processed_data/BRCA.PolyRisk.dna.v2.txt','Processed_data/BRCA.PolyRisk.amp.v2.txt','Processed_data/BRCA.PolyRisk.del.v2.txt',"BRCA",'Data/Hallmarks/nanostring.gl.txt','Figures/CPTAC.BRCA.hallmark.boxplots.pdf','Processed_data/CPTAC.BRCA.hallmark.scores.txt','Figures/Circos.CPTAC.BRCA.hallmarks.pdf',"Data/Meta_table/CPTAC-Pancan-Data_metatable.txt")

#args = c('Processed_data/PDAC.PolyRisk.dna.v2.txt','Processed_data/PDAC.PolyRisk.amp.v2.txt','Processed_data/PDAC.PolyRisk.del.v2.txt',"PDAC",'Data/Hallmarks/nanostring.gl.txt','Figures/CPTAC.PDAC.hallmark.boxplots.pdf','Processed_data/CPTAC.PDAC.hallmark.scores.txt','Figures/Circos.CPTAC.PDAC.hallmarks.pdf',"Data/Meta_table/CPTAC-Pancan-Data_metatable.txt")

#args = c('Processed_data/{can}.PolyRisk.dna.v2.txt','Processed_data/{can}.PolyRisk.amp.v2.txt','Processed_data/{can}.PolyRisk.del.v2.txt',"{can}",'Data/Hallmarks/nanostring.gl.txt','Figures/CPTAC.{can}.hallmark.boxplots.pdf','Processed_data/CPTAC.{can}.hallmark.scores.txt','Figures/Circos.CPTAC.{can}.hallmarks.pdf',"Data/Meta_table/CPTAC-Pancan-Data_metatable.txt")

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

p <- ggplot(mhall,aes(x=Var2,y=value))
p <- p + geom_boxplot()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))


pdf(args[6],height=5,width=8,useDingbats=F)
print(p)
dev.off()


#Not I just need to find the top 3 pathways and see if they are matching across these groups

mhall$value2 = abs(mhall$value)

tophall <-  data.frame(mhall %>%
  arrange(desc(value2)) %>%
  group_by(Var1) %>%
  slice(1:3))


colnames(mhall) <- c("ID","sectors","value","absvalue")

write.table(mhall,args[7],sep="\t",quote=F,row.names=F)

#NOW ONTO THE CIRCOS

#Figure out the pictures 
library(yaml)
data = yaml.load_file("./hallmarks.yaml")
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


pdf(args[8],heigh=10,width=10,useDingbats=F)
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




