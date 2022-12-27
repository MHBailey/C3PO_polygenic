library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)
library(viridis)
library(stringr)

'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)

#args = c('Data/Validation_CCLE/CCLE_gene_cn.csv','Figures/CCLE.validation.cnv.PTEN.pdf','Figures/CCLE.validation.cnv.MYC.pdf','Validation/CCLE.thresholded.cnv.tsv')


cnv <- fread(args[1])
ampo <- args[2]
delo <- args[3]
thresho <- args[4]



#Figure out a rough estimate for an AMP threshold 
mycname = colnames(cnv)[which(grepl("MYC ",colnames(cnv)))]
ampt = cnv %>% select(MYC=all_of(mycname))

pdf(ampo,height=5,width=5,useDingbats=F)
plot(density(ampt$MYC))
abline(v=1.2,col="red")
dev.off()

#Figure out the the rough estimate for DEL threshold 
ptenname = colnames(cnv)[which(grepl("PTEN ",colnames(cnv)))]
delt = cnv %>% select(PTEN=all_of(ptenname))

pdf(delo,height=5,width=5,useDingbats=F)
plot(density(delt$PTEN))
abline(v=0.8,col="blue")
dev.off()


####################################################################

cnvt = cnv
samples = cnv$V1
rownames(cnvt) = samples
cnvt$V1 = NULL
cnames <- str_split_fixed(colnames(cnvt)," \\(",2)[,1]
colnames(cnvt) <- cnames

cnvt[cnvt < .8] <- -1 
cnvt[cnvt >= .8 & cnvt <= 1.2] <- 0 
cnvt[cnvt > 1.2] <- 1 

write.table(cnvt,thresho,row.names=T,quote=F,sep="\t")

####################################################################






