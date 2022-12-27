library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)
library(stringr)
library(circlize)
library(entropy)

'%!in%' <- function(x,y)!('%in%'(x,y))

args = commandArgs(trailingOnly=TRUE)
 
#args <- c("Data/ssGSEA/Pancan_driver_DNA_ssGSEA-combined.gct","Data/ssGSEA/Pancan_driver_RNA_ssGSEA-combined.gct","Data/ssGSEA/Pancan_driver_protein_ssGSEA-combined.gct","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","./hallmarks.v2.yaml","Figures/ssGSEA.protein.circos.all.pdf")



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

############################DNA CIRCOS#################333
dna_c$Hallmarks = rownames(dna_c)
dna_m <- melt(dna_c)
dna_mc <- merge(dna_m,meta_lim,by.x ="variable", by.y="Proteome_Sample_ID",all.x=T)

#Now I need to grab the cancertypes 
cancers = unique(meta_lim$cohort)

#########################RNA CIRCOS#######################

#########################PROTEIN CIRCOS###################
protein_m <- melt(protein_c)
#These names seem quite odd. I'm going to convert to the . id and go from there: here is my attempt at this: 


protein_mc <- merge(protein_m,meta_lim,by.x="variable",by.y="Proteome_Sample_ID")
protein_mc$abs_value <- abs(protein_mc$value)
#Needs not missing values 
protein_nona <- na.omit(protein_mc)

tophall <-  data.frame(protein_mc %>%
  arrange(desc(value)) %>%
  group_by(variable) %>%
  slice(1:2))


library(yaml)
data = yaml.load_file(args[5])
hallmark_list = data$symbol
hallmark_name = sapply(hallmark_list, function(x) x$name)
hallmark_src = sapply(hallmark_list, function(x) x$src)

#Set up the link colors# 
colors = fread("colors.txt")

#ADD CT to tophall 

colnames(protein_nona) <- c("variable","sectors","value","cohort","Sample","CASE_ID","variable.y","variable2","abs_value")

pdf(args[6],heigh=10,width=10,useDingbats=F)
circos.par("track.height" = 0.3)
circos.initialize(protein_nona$sectors, xlim=c(min(protein_nona$value),max(protein_nona$value)))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
    image = EBImage::readImage(hallmark_src[CELL_META$sector.numeric.index])
    rasterImage(image,
        xleft = pos[1, 1] - 0.08, ybottom = pos[1, 2] - 0.08,
        xright = pos[1, 1] + 0.08, ytop = pos[1, 2]+ 0.08)
}, bg.border = 1, track.height = .25)
circos.trackHist(protein_nona$sectors, x = protein_nona$value,"track.height" = 0.15)

#Now I need to get all of the cancer types 
types_p <- unique(tophall$cohort)

#Now I need to add that funny variabile T1,T2,T3 
tophall$Var3 <- rep(c("T1","T2"))

t = "BRCA"
#ADD CT to tophall 
for(t in types_p){
    hey_p <- dcast(tophall[which(tophall$cohort == t),],variable~Var3,value.var="sectors")
    mycomb_p <- combn(unique(tophall[which(tophall$cohort == t),]$sectors),2)


    for(i in 1:dim(mycomb_p)[2]){
        cancol = colors[which(colors$Cancer == t),]$Color
        hey_p8 = NULL
        pair = mycomb_p[,i]
        p1 = pair[1]
        p2 = pair[2]
        hey_p2 = hey_p[which(hey_p$T1 == p1 & hey_p$T2 == p2),]
        hey_p3 = hey_p[which(hey_p$T2 == p1 & hey_p$T1 == p2),]

#        hey_p4 = hey_p[which(hey_p$T2 == p1 & hey_p$T3 == p2),]
#        hey_p5 = hey_p[which(hey_p$T3 == p1 & hey_p$T2 == p2),]

#        hey_p6 = hey_p[which(hey_p$T1 == p1 & hey_p$T3 == p2),]
#        hey_p7 = hey_p[which(hey_p$T3 == p1 & hey_p$T1 == p2),]

#        hey_p8 = rbind(hey_p2,hey_p3,hey_p4,hey_p5,hey_p6,hey_p7)
        hey_p8 = rbind(hey_p2,hey_p3)
        ids1 = tophall[which(tophall$cohort == t),][which(tophall[which(tophall$cohort == t),]$variable %in% hey_p8$variable & tophall[which(tophall$cohort == t),]$sectors == p1),]$value
        ids2 = tophall[which(tophall$cohort == t),][which(tophall[which(tophall$cohort == t),]$variable %in% hey_p8$variable & tophall[which(tophall$cohort == t),]$sectors == p2),]$value
        if(dim(hey_p8)[1] != 0){
            circos.link(p1,ids1,p2,ids2,h.ratio=.4,col=cancol)
        }
    }
}

dev.off()



#NOW BUILD ONE FOR ALL OF THE THE CANCER TYPES 

types_p <- unique(tophall$cohort)
tophall$Var3 <- rep(c("T1","T2"))

for(t in types_p){
    ofname = paste("Figures/ssGSEA.protein.circos.",t,".pdf",sep="")
    pdf(ofname,,heigh=10,width=10,useDingbats=F)
        circos.par("track.height" = 0.3)
        circos.initialize(protein_nona$sectors, xlim=c(min(protein_nona$value),max(protein_nona$value)))
        circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
            image = EBImage::readImage(hallmark_src[CELL_META$sector.numeric.index])
            rasterImage(image,
                xleft = pos[1, 1] - 0.08, ybottom = pos[1, 2] - 0.08,
                xright = pos[1, 1] + 0.08, ytop = pos[1, 2]+ 0.08)
        }, bg.border = 1, track.height = .25)
        circos.trackHist(protein_nona$sectors, x = protein_nona$value,"track.height" = 0.15)
    
    #Now I need to get all of the cancer types 
    
    #Now I need to add that funny variabile T1,T2,T3 
    
    #ADD CT to tophall 
        hey_p <- dcast(tophall[which(tophall$cohort == t),],variable~Var3,value.var="sectors")
        mycomb_p <- combn(unique(tophall[which(tophall$cohort == t),]$sectors),2)
    
    
        for(i in 1:dim(mycomb_p)[2]){
            cancol = colors[which(colors$Cancer == t),]$Color
            hey_p8 = NULL
            pair = mycomb_p[,i]
            p1 = pair[1]
            p2 = pair[2]
            hey_p2 = hey_p[which(hey_p$T1 == p1 & hey_p$T2 == p2),]
            hey_p3 = hey_p[which(hey_p$T2 == p1 & hey_p$T1 == p2),]
    
    #        hey_p4 = hey_p[which(hey_p$T2 == p1 & hey_p$T3 == p2),]
    #        hey_p5 = hey_p[which(hey_p$T3 == p1 & hey_p$T2 == p2),]
    
    #        hey_p6 = hey_p[which(hey_p$T1 == p1 & hey_p$T3 == p2),]
    #        hey_p7 = hey_p[which(hey_p$T3 == p1 & hey_p$T1 == p2),]
    
    #        hey_p8 = rbind(hey_p2,hey_p3,hey_p4,hey_p5,hey_p6,hey_p7)
            hey_p8 = rbind(hey_p2,hey_p3)
            ids1 = tophall[which(tophall$cohort == t),][which(tophall[which(tophall$cohort == t),]$variable %in% hey_p8$variable & tophall[which(tophall$cohort == t),]$sectors == p1),]$value
            ids2 = tophall[which(tophall$cohort == t),][which(tophall[which(tophall$cohort == t),]$variable %in% hey_p8$variable & tophall[which(tophall$cohort == t),]$sectors == p2),]$value
            if(dim(hey_p8)[1] != 0){
                circos.link(p1,ids1,p2,ids2,h.ratio=.4,col=cancol)
            }
        }
    dev.off()
}





########################SHANNONS E####################################
cans <- types_p
new <- data.frame(hallmarks = unique(protein_mc$sectors))

tophall <-  data.frame(protein_mc %>%
  arrange(desc(value)) %>%
  group_by(variable) %>%
  slice(1))

counts <- tophall %>% group_by(cohort) %>% count(sectors)



SHANNON1 = NULL
for(i in cans){
    ah <- counts[which(counts$cohort == i),]
    mynew <- merge(new,ah,by.x="hallmarks",by.y="sectors",all.x=T)
    mynew$CANCER <- i
    mynew[is.na(mynew)] <- 0
    freqs = mynew$n/sum(mynew$n)
    shannon = entropy.empirical(freqs)
    out = data.frame(CANCER=i,shannons=shannon,top=1)
    SHANNON1 = rbind(SHANNON1,out)
}


tophall <-  data.frame(protein_mc %>%
  arrange(desc(value)) %>%
  group_by(variable) %>%
  slice(2))

counts <- tophall %>% group_by(cohort) %>% count(sectors)



SHANNON2 = NULL
for(i in cans){
    ah <- counts[which(counts$cohort == i),]
    mynew <- merge(new,ah,by.x="hallmarks",by.y="sectors",all.x=T)
    mynew$CANCER <- i
    mynew[is.na(mynew)] <- 0
    freqs = mynew$n/sum(mynew$n)
    shannon = entropy.empirical(freqs)
    out = data.frame(CANCER=i,shannons=shannon,top=2)
    SHANNON2 = rbind(SHANNON2,out)
}



tophall <-  data.frame(protein_mc %>%
  arrange(desc(value)) %>%
  group_by(variable) %>%
  slice(3))

counts <- tophall %>% group_by(cohort) %>% count(sectors)



SHANNON3 = NULL
for(i in cans){
    ah <- counts[which(counts$cohort == i),]
    mynew <- merge(new,ah,by.x="hallmarks",by.y="sectors",all.x=T)
    mynew$CANCER <- i
    mynew[is.na(mynew)] <- 0
    freqs = mynew$n/sum(mynew$n)
    shannon = entropy.empirical(freqs)
    out = data.frame(CANCER=i,shannons=shannon,top=3)
    SHANNON3 = rbind(SHANNON3,out)
}


#Bring them all together
SHANNON = rbind(SHANNON1,SHANNON2,SHANNON3)

#I'll wante to order shannon before I plot 
SHANNON$ord = factor(SHANNON$CANCER,levels=SHANNON1[order(SHANNON1$shannons),]$CANCER)
#and I want to add specific colors: 
colors <- fread("colors.txt")
mycol <- colors$Color
names(mycol) <- factor(colors$Cancer,levels=SHANNON1[order(SHANNON1$shannons),]$CANCER)
shanplot <- merge(SHANNON,colors,by.x="CANCER",by.y="Cancer")

pdf("Figures/Entropy.v1.ssGSEA.pdf",height=5,width=10,useDingbats=F)
p <- ggplot(shanplot,aes(x=ord,y=shannons,fill=CANCER,group=top))
p <- p + geom_bar(stat="identity",position="dodge",color="white")
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),legend.position="none")
p <- p + scale_fill_manual(values=mycol)
p <- p + scale_color_manual("white")
p <- p + ylab("Shannon entropy for top three hallmarks")
p
dev.off()


