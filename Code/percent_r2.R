library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)
library(stringr)
library(ggridges)

'%!in%' <- function(x,y)!('%in%'(x,y))

args = commandArgs(trailingOnly=TRUE)

#args = c("Data/Meta_table/CPTAC-Pancan-Data_metatable.txt",'Figures/percent_captured.pdf',"Data/Somatic_mutation_wxs/can.genes.txt")


#cancer gene list
gl <- fread(args[3],header=F)

#Start with the TCGA validation data: 
#AMPLIFICATION 
vamp = list.files("./Validation",pattern="\\.SampleCorrelations.amp.v1.txt$")
ALL_amp = NULL 
for(i in vamp){
    if(!str_detect(i,"pancan")){
        print(i)
        if(startsWith(i,"CCLE")){
            cancer = str_split_fixed(i,"\\.",3)[,2]
        }else{
            cancer = str_split_fixed(i,"\\.",2)[,1]
        }
        fname = paste("./Validation/",i,sep="")
        y = fread(fname)
        y$CANCER = cancer 
        y$COHORT = ifelse(startsWith(i,"CCLE"),"CCLE","TCGA")
        ALL_amp = rbind(ALL_amp,y)
    }
}


vdel = list.files("./Validation",pattern="\\.SampleCorrelations.del.v1.txt$")
ALL_del = NULL
for(i in vdel){
    if(!str_detect(i,"pancan")){
        print(i)
        if(startsWith(i,"CCLE")){
            cancer = str_split_fixed(i,"\\.",3)[,2]
        }else{
            cancer = str_split_fixed(i,"\\.",2)[,1]
        }
        fname = paste("./Validation/",i,sep="")
        y = fread(fname)
        y$CANCER = cancer
        y$COHORT = ifelse(startsWith(i,"CCLE"),"CCLE","TCGA")
        ALL_del = rbind(ALL_del,y)
    }
}

vdna = list.files("./Validation",pattern="\\.SampleCorrelations.dna.v1.txt$")
ALL_dna = NULL
for(i in vdna){
    if(!str_detect(i,"pancan")){
        print(i)
        if(startsWith(i,"CCLE")){
            cancer = str_split_fixed(i,"\\.",3)[,2]
        }else{
            cancer = str_split_fixed(i,"\\.",2)[,1]
        }
        fname = paste("./Validation/",i,sep="")
        y = fread(fname)
        y$CANCER = cancer
        y$COHORT = ifelse(startsWith(i,"CCLE"),"CCLE","TCGA")
        ALL_dna = rbind(ALL_dna,y)
    }
}

vcombined = list.files("./Validation",pattern="\\.SampleCorrelations.combined.v1.txt$")
ALL_combined = NULL
for(i in vcombined){
    if(!str_detect(i,"pancan")){
        print(i)
        if(startsWith(i,"CCLE")){
            cancer = str_split_fixed(i,"\\.",3)[,2]
        }else{
            cancer = str_split_fixed(i,"\\.",2)[,1]
        }
        fname = paste("./Validation/",i,sep="")
        y = fread(fname)
        y$CANCER = cancer
        y$COHORT = ifelse(startsWith(i,"CCLE"),"CCLE","TCGA")
        ALL_combined = rbind(ALL_combined,y)
    }
}

prep_amp = data.frame(ALL_amp %>% select(Protein,Coef,P.value,Substrate,CANCER,COHORT))
prep_del = data.frame(ALL_del %>% select(Protein,Coef,P.value,Substrate,CANCER,COHORT))
prep_dna = data.frame(ALL_dna %>% select(Protein,Coef,P.value,Substrate,CANCER,COHORT))

ALL <- rbind(prep_amp,prep_del,prep_dna,ALL_combined)
ALL$NEW_Y = paste(ALL$COHORT,ALL$CANCER,ALL$Substrate,sep="_")

ALL_noC <- ALL[which(ALL$Substrate %!in% c("COMBINED")),]
ALL_C <- ALL[which(ALL$Substrate %in% c("COMBINED")),]

pdf("Validation/TCGA.CCLE.coefs.pdf",height=7,width=5,useDingbats=F)
#Now to try to build the plot for it. 
p <- ggplot(ALL_noC, aes(x = Coef, y = NEW_Y, fill=Substrate))
cp <- p + theme_minimal()
p <- p + geom_density_ridges(alpha=.3)
p
dev.off()
#Now Just combined
pdf("Validation/TCGA.CCLE.combined.coefs.pdf",height=7,width=5,useDingbats=F)
p <- ggplot(ALL_C, aes(x = Coef, y = NEW_Y, fill=Substrate))
p <- p + theme_minimal()
p <- p + geom_density_ridges(alpha=.3)
p
dev.off()

#Now to summarize
sum_C <- data.frame(ALL_C %>% group_by(Substrate,CANCER,COHORT) %>% summarize(Mean=mean(Coef,na.rm=T),SD=sd(Coef,na.rm=T),MAX=max(Coef),MIN=min(Coef),Mean2=mean(Coef^2,na.rm=T),SD2=sd(Coef^2,na.rm=T),MAX2=max(Coef^2),MIN2=min(Coef^2)))
sum_noC <- data.frame(ALL_noC %>% group_by(Substrate,CANCER,COHORT) %>% summarize(Mean=mean(Coef,na.rm=T),SD=sd(Coef,na.rm=T),MAX=max(Coef,na.rm=T),MIN=min(Coef,na.rm=T),Mean2=mean(Coef^2,na.rm=T),SD2=sd(Coef^2,na.rm=T),MAX2=max(Coef^2),MIN2=min(Coef^2)))

TCGA <- sum_noC[which(sum_noC$COHORT == "TCGA"),]
TCGA_c <- sum_C[which(sum_C$COHORT == "TCGA"),]
#What are the easiest genes to predict? 

topgenes_p <-  data.frame(unique(ALL_noC) %>%
  arrange(P.value) %>%
  group_by(CANCER,COHORT,Substrate) %>%
  slice(1:100))

topgenes_p <-  data.frame(ALL_C %>%
  arrange(P.value) %>%
  group_by(CANCER,COHORT,Substrate) %>%
  slice(1:10))


#Then get the CCLE proteomics data 

#Then get into the CPTAC model (which is too circular!) but may be helpful 
#Start with the TCGA validation data: 
#AMPLIFICATION 
vamp = list.files("./Processed_data",pattern="\\.SampleCorrelations.amp.v1.txt$")
ALL_amp = NULL
for(i in vamp){
    if(!str_detect(i,"pancan")){
        print(i)
        cancer = str_split_fixed(i,"\\.",2)[,1]
        fname = paste("./Processed_data/",i,sep="")
        y = fread(fname)
        y$CANCER = cancer
        ALL_amp = rbind(ALL_amp,y)
    }
}
ALL_amp$Substrate="AMP"

vdel = list.files("./Processed_data",pattern="\\.SampleCorrelations.del.v1.txt$")
ALL_del = NULL
for(i in vdel){
    if(!str_detect(i,"pancan")){
        print(i)
        cancer = str_split_fixed(i,"\\.",2)[,1]
        fname = paste("./Processed_data/",i,sep="")
        y = fread(fname)
        y$CANCER = cancer
        ALL_del = rbind(ALL_del,y)
    }
}
ALL_del$Substrate = "DEL"


vdna = list.files("./Processed_data",pattern="\\.SampleCorrelations.dna.v1.txt$")
ALL_dna = NULL
for(i in vdna){
    if(!str_detect(i,"pancan")){
        print(i)
        cancer = str_split_fixed(i,"\\.",2)[,1]
        fname = paste("./Processed_data/",i,sep="")
        y = fread(fname)
        y$CANCER = cancer
        ALL_dna = rbind(ALL_dna,y)
    }
}
ALL_dna$Substrate = "DNA"

vcombined = list.files("./Processed_data",pattern="\\.SampleCorrelations.combined.v1.txt$")
ALL_combined = NULL
for(i in vcombined){
    if(!str_detect(i,"pancan")){
        print(i)
        cancer = str_split_fixed(i,"\\.",2)[,1]
        fname = paste("./Processed_data/",i,sep="")
        y = fread(fname)
        y$CANCER = cancer
        ALL_combined = rbind(ALL_combined,y)
    }
}
ALL_combined$Substrate ="Combined"

prep_amp = data.frame(ALL_amp %>% select(Protein,Coef,P.value,Substrate,CANCER))
prep_del = data.frame(ALL_del %>% select(Protein,Coef,P.value,Substrate,CANCER))
prep_dna = data.frame(ALL_dna %>% select(Protein,Coef,P.value,Substrate,CANCER))

ALL <- rbind(prep_amp,prep_del,prep_dna,ALL_combined)
ALL$NEW_Y = paste(ALL$CANCER,ALL$Substrate,sep="_")

ALL_noC <- ALL[which(ALL$Substrate %!in% c("Combined")),]
ALL_C <- ALL[which(ALL$Substrate %in% c("Combined")),]

pdf("Figures/CPTAC.percent.explained.pdf",height=7,width=5)
#Now to try to build the plot for it. 
p <- ggplot(ALL_noC, aes(x = Coef, y = NEW_Y, fill=Substrate))
p <- p + theme_minimal()
p <- p + geom_density_ridges(alpha=.3)
print(p)
dev.off()
#Now Just combined
pdf("Figures/CPTAC.percent.explained.combined.pdf",height=7,width=5)
p <- ggplot(ALL_C, aes(x = Coef, y = NEW_Y, fill=Substrate))
p <- p + theme_minimal()
p <- p + geom_density_ridges(alpha=.3)
print(p)
dev.off()

#Now to summarize
CPTAC_C <- data.frame(ALL_C %>% group_by(CANCER) %>% summarize(Mean=mean(Coef),SD=sd(Coef),MAX=max(Coef),MIN=min(Coef),Mean2=mean(Coef^2,na.rm=T),SD2=sd(Coef^2,na.rm=T),MAX2=max(Coef^2),MIN2=min(Coef^2)))
CPTAC_NOC <- data.frame(ALL_noC %>% group_by(Substrate,CANCER) %>% summarize(Mean=mean(Coef),SD=sd(Coef),MAX=max(Coef),MIN=min(Coef),Mean2=mean(Coef^2,na.rm=T),SD2=sd(Coef^2,na.rm=T),MAX2=max(Coef^2),MIN2=min(Coef^2)))


topgenes_p <-  data.frame(ALL_noC %>%
  arrange(P.value) %>%
  group_by(CANCER,Substrate) %>%
  slice(1:10))

topgenes_p <-  data.frame(ALL_C %>%
  arrange(P.value) %>%
  group_by(CANCER,Substrate) %>%
  slice(1:10))

select(CPTAC_C,c("CANCER","Mean2"))



#This is a graph with a different set.... 
ALL_C$Driver <- ifelse(ALL_C$Protein %in% gl$V1,1,0 )

p <- ggplot(ALL_C[which(ALL_C$CANCER == "PDAC"),],aes(x=Driver,y=Coef,fill=factor(Driver),color=factor(Driver)))
p <- p + geom_boxplot()
p

