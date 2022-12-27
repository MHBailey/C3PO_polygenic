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

halls <- list.files("./Processed_data",pattern="\\.hallmark.scores.txt$")

ALL_hall <- NULL 
for(i in halls){
    if(!str_detect(i,"PANCAN")){
        print(i)
        cancer = str_split_fixed(i,"\\.",3)[,2]
        fname = paste("./Processed_Data/",i,sep="")
        y = fread(fname)
        y$CANCER = cancer
        ALL_hall = rbind(ALL_hall,y)
    }
}

#Get the top hits for every sample 
tophall <-  data.frame(ALL_hall %>%
  arrange(desc(absvalue)) %>%
  group_by(CANCER,ID) %>%
  slice(1))

#Get the top counts and make sure that there are 0's for missing hallmarks
counts <- tophall %>% group_by(CANCER) %>% count(sectors) 

#Now I need to make sure to re-introduce some 0's 
new <- data.frame(hallmarks = unique(ALL_hall$sectors))
new$nn = 0

#Now calculate for Shannons and add to a large table 
cans <- unique(counts$CANCER)

SHANNON1 = NULL
for(i in cans){
    ah <- counts[which(counts$CANCER == i),]
    mynew <- merge(new,ah,by.x="hallmarks",by.y="sectors",all.x=T)
    mynew$CANCER <- i 
    mynew[is.na(mynew)] <- 0
    freqs = mynew$n/sum(mynew$n)
    shannon = entropy.empirical(freqs)
    out = data.frame(CANCER=i,shannons=shannon,top=1)
    SHANNON1 = rbind(SHANNON1,out)
}


tophall <-  data.frame(ALL_hall %>%
  arrange(desc(absvalue)) %>%
  group_by(CANCER,ID) %>%
  slice(2))

#Get the top counts and make sure that there are 0's for missing hallmarks
counts <- tophall %>% group_by(CANCER) %>% count(sectors)

#Now I need to make sure to re-introduce some 0's 
new <- data.frame(hallmarks = unique(ALL_hall$sectors))
new$nn = 0

#Now calculate for Shannons and add to a large table 
cans <- unique(counts$CANCER)

SHANNON2 = NULL
for(i in cans){
    ah <- counts[which(counts$CANCER == i),]
    mynew <- merge(new,ah,by.x="hallmarks",by.y="sectors",all.x=T)
    mynew$CANCER <- i
    mynew[is.na(mynew)] <- 0
    freqs = mynew$n/sum(mynew$n)
    shannon = entropy.empirical(freqs)
    out = data.frame(CANCER=i,shannons=shannon,top=2)
    SHANNON2 = rbind(SHANNON2,out)
}


tophall <-  data.frame(ALL_hall %>%
  arrange(desc(absvalue)) %>%
  group_by(CANCER,ID) %>%
  slice(3))

#Get the top counts and make sure that there are 0's for missing hallmarks
counts <- tophall %>% group_by(CANCER) %>% count(sectors)

#Now I need to make sure to re-introduce some 0's 
new <- data.frame(hallmarks = unique(ALL_hall$sectors))
new$nn = 0

#Now calculate for Shannons and add to a large table 
cans <- unique(counts$CANCER)

SHANNON3 = NULL
for(i in cans){
    ah <- counts[which(counts$CANCER == i),]
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

pdf("Figures/Entropy.v1.pdf",height=5,width=10,useDingbats=F)
p <- ggplot(shanplot,aes(x=ord,y=shannons,fill=CANCER,group=top))
p <- p + geom_bar(stat="identity",position="dodge",color="white")
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),legend.position="none")
p <- p + scale_fill_manual(values=mycol)
p <- p + scale_color_manual("white")
p <- p + ylab("Shannon entropy for top three hallmarks")
p
dev.off()

#NOW I want to plot them and pull out the variability 

tophall <-  data.frame(ALL_hall %>%
  arrange(desc(absvalue)) %>%
  group_by(CANCER,ID) %>%
  slice(1:3))

tophall$Var3 <- rep(c("T1","T2","T3"))
hallplot <- tophall %>% group_by(sectors,Var3,CANCER) %>% count() 
pdf("Figures/summed_top_hallmarks.v1.pdf",height=6,width=10,useDingbats=F)
p <- ggplot(hallplot,aes(x=(interaction(Var3,sectors)),y=n,fill=CANCER))
p <- p + geom_histogram(stat="identity")
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),legend.position="none")
p <- p + scale_fill_manual(values=mycol)
p <- p + scale_color_manual("white")
p
dev.off()




