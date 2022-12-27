library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)
library(stringr)
library(circlize)

'%!in%' <- function(x,y)!('%in%'(x,y))

args = commandArgs(trailingOnly=TRUE)

#args = c('./hallmarks.v2.yaml',"Data/Meta_table/CPTAC-Pancan-Data_metatable.txt",'Figures/top3_all_circos.pdf')


x = list.files("./Processed_data",pattern="\\.hallmark.scores.txt$")
ALL = NULL 
for(i in x){
    if(!str_detect(i,"PANCAN")){
        print(i)
        fname = paste("./Processed_data/",i,sep="")
        y = fread(fname)
        ALL = rbind(ALL,y)
    }
}

tophall <-  data.frame(ALL %>%
  arrange(desc(absvalue)) %>%
  group_by(ID) %>%
  slice(1:3))


#Now I want to put it together: 
library(yaml)
data = yaml.load_file(args[1])
hallmark_list = data$symbol
hallmark_name = sapply(hallmark_list, function(x) x$name)
hallmark_src = sapply(hallmark_list, function(x) x$src)



#Set up the link colors# 
colors = fread("colors.txt")

#ADD CT to tophall 
meta = fread(args[2])
metafav = data.frame(meta %>% select("Proteome_Sample_ID","cohort"))
tophall$Var3 <- rep(c("T1","T2","T3"))
tophall_merged <- merge(tophall,metafav,by.x="ID",by.y="Proteome_Sample_ID")
types = unique(tophall_merged$cohort)

pdf(args[3],heigh=10,width=10,useDingbats=F)
circos.par("track.height" = 0.3)
circos.initialize(ALL$sectors, xlim=c(min(ALL$value),max(ALL$value)))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
    image = EBImage::readImage(hallmark_src[CELL_META$sector.numeric.index])
    rasterImage(image,
        xleft = pos[1, 1] - 0.08, ybottom = pos[1, 2] - 0.08,
        xright = pos[1, 1] + 0.08, ytop = pos[1, 2]+ 0.08)
}, bg.border = 1, track.height = .25)
circos.trackHist(ALL$sectors, x = ALL$value,"track.height" = 0.15)

#ADD CT to tophall 
for(t in types){
    hey <- dcast(tophall_merged[which(tophall_merged$cohort == t),],ID~Var3,value.var="sectors")
    mycomb <- combn(unique(tophall_merged[which(tophall_merged$cohort == t),]$sectors),2)


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
        ids1 = tophall_merged[which(tophall_merged$cohort == t),][which(tophall_merged[which(tophall_merged$cohort == t),]$ID %in% hey8$ID & tophall_merged[which(tophall_merged$cohort == t),]$sectors == p1),]$value
        ids2 = tophall_merged[which(tophall_merged$cohort == t),][which(tophall_merged[which(tophall_merged$cohort == t),]$ID %in% hey8$ID & tophall_merged[which(tophall_merged$cohort == t),]$sectors == p2),]$value
        if(dim(hey8)[1] != 0){
            circos.link(p1,ids1,p2,ids2,h.ratio=.4,col=cancol)
        }
    }
}
legend("topleft", pch = 1, legend = "Everything")
dev.off()



