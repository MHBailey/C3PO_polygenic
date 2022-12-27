library(data.table) 
library(reshape2)

vamp = list.files("./",pattern="\\.hallmark.scores.txt$")
ALL_c3po = NULL
for(i in vamp){
        y = fread(i)
        y$value=NULL
        z = dcast(y,ID~sectors,value.var="absvalue")
        ALL_c3po = rbind(ALL_c3po,z)
}

write.table(ALL_c3po,"Combined.C3PO.hallmarkscores.txt",sep="\t",quote=F,row.names=F)
