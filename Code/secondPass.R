library(data.table)
library(MASS)
library(ggplot2)
library(stringr)
library(reshape2)
library(viridis)
library(dplyr)
library(pROC)
library(keras)
library(tfdatasets)

'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)


args <- c("Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Gene_expression/ALL_RNA-Seq_Expr_WashU_FPKM_UQ_annotation.tsv","Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_imputed_gene_level.tsv","Data/Clinical_data/clinical_Pan-cancer.Jan2022.tsv","Data/Drivers/mutual_coocurrance.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/Hallmarks/InducingAngiogenesis.tmp.txt","Data/GeneLists/299.genes.txt","Processed_data/angio.topHITs.TODAY.txt")


maf=fread(args[1])
rna=fread(args[2])
prot=fread(args[3])
clin=fread(args[4])
mutco=fread(args[5])
meta=fread(args[6])
hall=fread(args[7],header=F)
g299=fread(args[8],header=T)
ofname=args[9]

genes <- hall$V1
g299 <- g299$"Approved symbol" 

#CHANGE THIS CODE ONCE WITH HAVE THE NEWEST VERSION
dMAF <- maf[which(maf$Hugo_Symbol %in% g299),] 
tokeep <- c("Frame_Shift_Del","Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation")
dMAF2 <- dMAF[which(dMAF$Variant_Classification %in% tokeep),] 
singles <- unique(dMAF2$Hugo_Symbol) 

g = "ANG"
g = "TOP1"
s = "TP53"


#HERE I WANT TO START TO BUILD THE THING IF I CAN (a giant 0,1) MATRIX OF DRIVER GENES and proteins with the the variables that I want, including AGE, SEX, CANCER TYPE, #AND USE GRADIENT DECENT IN PYTORCH 
#SOME ODD FEATURES
#ODD FEATURES, TOP50% for RNAseq or not 
#METHYLATION HITS 

MUTS <- data.frame(unique(dMAF2 %>% select(Tumor_Sample_Barcode,Hugo_Symbol,COHORT)))
dMUTS <- dcast(Tumor_Sample_Barcode ~ Hugo_Symbol,data=MUTS,fun.aggregate=length)

META <- data.frame(meta %>% select(Sample,Proteome_Sample_ID,Proteome_Normal_Sample_ID,WXS,cohort,age,sex))
METAT <- data.frame(META[which(META$Proteome_Sample_ID !=""),] %>% select(Proteome_Sample_ID,cohort,age,sex))
#METAT$is_tumor = 1
#colnames(METAT) <- c("PROTid","cohort","age","sex","is_tumor")
colnames(METAT) <- c("PROTid","cohort","age","sex")
METAN <- data.frame(META[which(META$Proteome_Normal_Sample_ID != ""),] %>% select(Proteome_Normal_Sample_ID,cohort,age,sex))
#METAN$is_tumor = 0
#colnames(METAN) <- c("PROTid","cohort","age","sex","is_tumor")
colnames(METAN) <- c("PROTid","cohort","age","sex")
#METAX <- rbind(METAT,METAN)
METAX <- METAT
METAM <- data.frame(META %>% select(WXS,Proteome_Sample_ID)) #TO ADD TO MAF
mm <- merge(dMUTS,METAM,by.x="Tumor_Sample_Barcode",by.y="WXS",all=T,)

prott <- prot 
GN <- prot$external_gene_name
prott$external_gene_name = NULL 
prottt <- data.frame(t(prott))
colnames(prottt) <- GN
prottt$PSAMP <- rownames(prottt)

gprottt <- data.frame(prottt %>% select(TOP1,PSAMP))
pmm <- merge(gprottt,mm,by.x="PSAMP",by.y="Proteome_Sample_ID",all.x=T)
pmmc <- merge(pmm,METAX,by.x="PSAMP",by.y="PROTid")

lung <- pmmc[which(pmmc$cohort=="LUAD"),]
lung[is.na(lung)] <- 0


lung$cohort <- NULL
lung$sex <- ifelse(lung$sex == "Male", 1, 2)
lung$PSAMP <- NULL 
lung$Tumor_Sample_Barcode <- NULL


sample = sample(nrow(lung), 75)
lung_train = lung[sample, ]
train_labels <- as.vector(lung_train$TOP1)
train_colors <- as.vector(lung_train$is_tumor)
lung_train$TOP1 <- NULL

lung_test = lung[seq(1,nrow(lung),1) %!in% sample,]
test_labels <- as.vector(lung_test$TOP1)
test_colors <- as.vector(lung_test$is_tumor)
lung_test$TOP1 <- NULL


mynames = colnames(lung_train)

train_df <- lung_train %>% as_tibble(.name_repair = "minimal") %>%  setNames(mynames)  %>% mutate(label = train_labels)
test_df <- lung_test %>% as_tibble(.name_repair = "minimal") %>%  setNames(mynames)   %>% mutate(label = test_labels)

train_dfs <- train_df %>% select(which(!colSums(train_df, na.rm=TRUE) %in% 0))

spec <- feature_spec(train_dfs, label ~ . ) %>% 
  step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>% 
  fit()

input <- layer_input_from_dataset(train_dfs %>% select(-label))

output <- input %>% 
  layer_dense_features(dense_features(spec)) %>% 
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 1) 

model <- keras_model(input, output)

summary(model)

model %>% 
  compile(
    loss = "mse",
    optimizer = optimizer_rmsprop(),
    metrics = list("mean_absolute_error")
  )


build_model <- function() {
  input <- layer_input_from_dataset(train_df %>% select(-label))
  
  output <- input %>% 
    layer_dense_features(dense_features(spec)) %>% 
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 1) 
  
  model <- keras_model(input, output)
  
  model %>% 
    compile(
      loss = "mse",
      optimizer = optimizer_rmsprop(),
      metrics = list("mean_absolute_error")
    )
  
  model
}


#THIS IS THE LONGER RUN 
print_dot_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 80 == 0) cat("\n")
    cat(".")
  }
)    

model <- build_model()

history <- model %>% fit(
  x = train_df %>% select(-label),
  y = train_df$label,
  epochs = 300,
  validation_split = 0.2,
  verbose = 0,
  callbacks = list(print_dot_callback)
)

test_predictions <- model %>% predict(test_df %>% select(-label))

plot(history)
plot(test_predictions[,1],test_labels,col=test_colors+1)
ggplot(lung,aes(x=as.factor(is_tumor),y=TOP1)) + geom_boxplot()
cor.test(test_predictions[,1],test_labels)


#THI IS THE EARLY STOP
early_stop <- callback_early_stopping(monitor = "val_loss", patience = 20)

model <- build_model()

history <- model %>% fit(
  x = train_df %>% select(-label),
  y = train_df$label,
  epochs = 500,
  validation_split = 0.2,
  verbose = 0,
  callbacks = list(early_stop)
)



test_predictions <- model %>% predict(test_df %>% select(-label))

plot(test_predictions[,1],test_labels)




#####IS THE GOAL TO PREDICT IT... YES~ 













OUTPUT = NULL
for(g in genes){
    for(s in singles){
        if(g %in% prot$external_gene_name){
            df <- prot[which(prot$external_gene_name == g),]
            df$external_gene_name <- NULL
            mdf <- melt(df)
            muts <- dMAF2[which(dMAF2$Hugo_Symbol == s),]$Tumor_Sample_Barcode
            protmuts <- meta[which(meta$Sample %in% muts),]$Proteome_Sample_ID
            mdf$MUTS <- ifelse(mdf$variable %in% protmuts,"Mut","WT")
            protmutsWT <- meta[which(meta$Sample %!in% muts),]$Proteome_Sample_ID
            WTtmpNorm <- meta[which(meta$Proteome_Sample_ID %in% protmutsWT),]$Proteome_Normal_Sample_ID
            WTnormal <- WTtmpNorm[which(WTtmpNorm != "")] 
            WTtumor <- meta[which(meta$Proteome_Sample_ID %in% protmutsWT),]$Proteome_Sample_ID
            ttest <- t.test(value~MUTS,data=mdf)
            meanWTnorm <- mean(mdf[which(mdf$variable %in% WTnormal),]$value)
            meanWTtum <- mean(mdf[which(mdf$variable %in% WTtumor),]$value)
            out <- data.frame("Driver"=s,"Hallmark"=g,"pval"=ttest$p.value,"MUTmean"=ttest$estimate[1], "WTmeanTN"=ttest$estimate[2],"WTmeanN"=meanWTnorm,"WTmeanT"=meanWTtum)
            OUTPUT = rbind(OUTPUT,out)
        }
    }
}


#IS_TUMOR Gradient decent 
#One layer at a time 
#One protein at a time

ordOUT <- OUTPUT[order(OUTPUT$pval),]
ordOUT$DIFF <- abs(ordOUT$MUTmean-ordOUT$WTmeanT)
ord2OUT <- ordOUT[order(ordOUT$DIFF, decreasing = T),]

for(g in genes){
    for(s in singles){
        print(paste(g,s,sep=" "))
        df <- prot[which(prot$external_gene_name == g),]
        df$external_gene_name <- NULL
        mdf <- melt(df)
        mdf$TP53 = ifelse(mdf$variable %in% tp53samps3, "M", NA)#This will represent the mutated group (single genes) or cooc. v wt
        mdf$TP53 = ifelse(mdf$variable %in% xtp53samps3, "W", mdf$TP53) #This will represent the more wt group 
        sing = unique(maf[which(maf$Hugo_Symbol == s),]$Tumor_Sample_Barcode)
        xsing = setdiff(allsamps,sing)
        sing3 <- meta[which(meta$Sample %in% sing),]$Proteome_Sample_ID # " 
        xsing3 <- meta[which(meta$Sample %in% xsing),]$Proteome_Sample_ID
        mdf$Grouped = ifelse(mdf$variable %in% sing3, "M", NA)#This will represent the mutated group (single genes) or cooc. v wt
        mdf$Grouped = ifelse(mdf$variable %in% xsing3, "W", mdf$Grouped) #This will represent the more wt group
        mdf$MUTCO = ifelse(mdf$TP53 == "W" & mdf$Grouped == "M", "M", NA) #This grabs ths TP53wt and calls to SINGmut
        mdf$MUTCO = ifelse(mdf$TP53 == "W" & mdf$Grouped == "W", "W", mdf$MUTCO) #This grabs the TP53wt and SINGwt
        mutco_mdf <- mdf[which(!is.na(mdf$MUTCO)),]
        mutco_triple <- lda(formula=MUTCO~value,data=mutco_mdf)
        direction <- ifelse(mutco_triple$means[2] < mutco_triple$means[1],"RespIsSmallerInMut","RespIsLargerInMut")
        out = data.frame(Genes=g,Sing=s,LDAsvd=mutco_triple$svd,direction=direction)
        CO = rbind(CO,out)
    }
}

topsSING <- head(CO[order(CO$LDAsvd, decreasing =T),],500)



