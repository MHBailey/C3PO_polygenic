library(data.table)
library(MASS)
library(ggplot2)
library(stringr)
library(reshape2)
library(viridis)
library(dplyr)
library(pROC)


'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)


args <- c("Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Gene_expression/ALL_RNA-Seq_Expr_WashU_FPKM_UQ_annotation.tsv","Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_imputed_gene_level.tsv","Data/Clinical_data/clinical_Pan-cancer.Jan2022.tsv","Data/Drivers/mutual_coocurrance.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/Hallmarks/InducingAngiogenesis.tmp.txt","Processed_data/angio.topHITs.TODAY.txt")


maf=fread(args[1])
rna=fread(args[2])
prot=fread(args[3])
clin=fread(args[4])
mutco=fread(args[5])
meta=fread(args[6])
hall=fread(args[7],header=F)
g299=fread(args[8],header=T)
ofname=args[9]


#NOTE: Step 0: Filter mutations 
#TODO: yup


#Step 1: Find the samples with somatic mutations in TP53 
#NOTE: I should/can change this to only include Collins/Eduards/Amelia's predictions
allsamps <- unique(maf$Tumor_Sample_Barcode)
tp53samps <- unique(maf[which(maf$Hugo_Symbol == "TP53"),]$Tumor_Sample_Barcode)
xtp53samps <- setdiff(allsamps,tp53samps)

#Step 2: Split the Proteomics data up into groups and look for the most differentiated
#NOTE: This can be done a couple different ways. I'm going to start with LDA 
#NOTE: Also try to build some cross validation sets in this. 70 20 10 on this front. 
#Make sure that the column label match between proteins and mutations..
#All _ to . and all - to _ 

#TODO: I need to change this code to pull in the Meta data idenfiers instead of guessing. 
allsamps3 <- meta[which(meta$Sample %in% allsamps),]$Proteome_Sample_ID #Replace -
tp53samps3 <- meta[which(meta$Sample %in% tp53samps),]$Proteome_Sample_ID # " 
xtp53samps3 <- meta[which(meta$Sample %in% xtp53samps),]$Proteome_Sample_ID # " 

#Change the column names in prot 
genes <- prot$external_gene_name

#This just goes through and compares TP53 and TP53wt 
g = "NDUFS6"
look <- NULL
for(g in genes){
    df <- prot[which(prot$external_gene_name == g),]
    df$external_gene_name <- NULL
    mdf <- melt(df)
    mdf$Grouped = ifelse(mdf$variable %in% tp53samps3, "M", NA)#This will represent the mutated group (single genes) or cooc. v wt
    mdf$Grouped = ifelse(mdf$variable %in% xtp53samps3, "W", mdf$Grouped) #This will represent the more wt group 
    clean_mdf <- mdf[which(!is.na(mdf$Grouped)),]
    this_triple <- lda(formula=Grouped~value,data=clean_mdf)
    direction <- ifelse(this_triple$means[2] < this_triple$means[1],"RepIsLarger","RespIsSmaller")
    out = data.frame(Genes=g,LDAsvd=this_triple$svd,direction=direction)
    look = rbind(look,out)
}

#This not takes the top hit from TP53 and looks at which looks best from T/N comparisons
topsTP53 <- head(look[order(look$LDAsvd, decreasing =T),],100)
tumor <- na.omit(unique(meta$Proteome_Sample_ID))
normal <- na.omit(unique(meta$Proteome_Normal_Sample_ID))

lookTN <- NULL
for(g in topsTP53$Genes){
    df <- prot[which(prot$external_gene_name == g),]
    tG <- topsTP53[which(topsTP53$Genes == g),]
    df$external_gene_name <- NULL
    mdf <- melt(df)
    mdf$Tissue = ifelse(mdf$variable %in% normal,"N","T")
    this_triple <- lda(formula=Tissue~value,data=mdf)
    direction <- ifelse(this_triple$means[2] < this_triple$means[1],"RespIsSmall","RespIsLarger")
    out = data.frame(Genes=g,LDAgene=tG$LDAsvd,LDAsvd=this_triple$svd,dirMUT=tG$direction,dirNAT=direction)
    lookTN = rbind(lookTN,out)
}

lookTN$LDAsum = lookTN$LDAgene+lookTN$LDAsvd

#Does this makes sense with a plot: 
#Think throught he co-occuring events 
ordLDA <- head(lookTN[order(lookTN$LDAsum, decreasing =T),],100)

#ID by Cohort
tcode <- na.omit(data.frame(meta %>% select(Proteome_Sample_ID,cohort)))
ncode <- na.omit(data.frame(meta %>% select(Proteome_Normal_Sample_ID,cohort)))
ncode <- ncode[which(ncode$Proteome_Normal_Sample_ID != ""),]
colnames(tcode) <- c("ID","cohort")
colnames(ncode) <- c("ID","cohort")
codes <- rbind(tcode,ncode)


favs = c("ILK","COPB2","FOXK1","HSPA5") #FIRST PASS 
favs2 <- c("NUP93","TOP1","CSE1L")

#Cell proliferation and met in TNBC NUP93 (2020)
#Promotes migration (CSE1L) 
#DNA damage: Topoisomerases (TOP1)

g = "PALM2AKAP2"
g = "SNTB2A"
g = "AK1"
g = "RPS5"
g = "PEBP1"
g = "NUP93"
g = "DAG1"
g = "TOP1"
g = "PRKDC"
g = "GGH"
g = "CSE1L"
g = "NUP155"
g = "DDB2"
g = "PDE12"

#for(g in favs){
    df <- prot[which(prot$external_gene_name == g),]
    df$external_gene_name <- NULL
    mdf <- melt(df)
    mdf$TN = ifelse(mdf$variable %in% normal,"N","T")
    mdf$MW = ifelse(mdf$variable %in% tp53samps3,"M","W")
    mdf$Group = paste(mdf$TN,mdf$MW,sep="_")
    mdf2 <- merge(mdf,codes,by.x="variable",by.y="ID")
#}

p <- ggplot(mdf2,aes(x=Group,y=value,fill=cohort))
p <- p + geom_boxplot()
p <- p + theme_classic()
p <- p + ggtitle(g)
p




this_triple <- lda(formula=Tissue~value,data=mdf)
plot(this_triple)



#Think through the Mutually exclusive events, so for the first pass I'm thinking to break this down into TP53WT v other categories. And then it would also be interesting to see Kidney data And not TP53mut 
CO <- NULL

#Think through the singles 
singles <- c("KRAS","PIK3CA","PTEN","EGFR","VHL")
#g = "A1BG"
#s = "VHL"

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




##### THIS NEXT STEP IS TO LOOK AT THIS MODEL without the TP53 CAVIATE? 
SING3 <- NULL

#Think through the singles 
singles <- c("TP53","KRAS","PIK3CA","PTEN","EGFR","VHL")
#g = "A1BG"
#s = "VHL"
#TODO: We need to split this up by cancer type 
cancers = unique(meta$cohort)

for(can in cancers){ 
    for(g in genes){
        for(s in singles){
            #print(paste(g,s,sep=" "))
            df <- prot[which(prot$external_gene_name == g),]
            df$external_gene_name <- NULL
P            mdf <- melt(df)
            mdf$TP53 = ifelse(mdf$variable %in% tp53samps3, "M", NA)#This will represent the mutated group (single genes) or cooc. v wt
            mdf$TP53 = ifelse(mdf$variable %in% xtp53samps3, "W", mdf$TP53) #This will represent the more wt group 
            sing = unique(maf[which(maf$Hugo_Symbol == s),]$Tumor_Sample_Barcode)
            xsing = setdiff(allsamps,sing)
            sing3 <- meta[which(meta$Sample %in% sing & meta$cohort == can),]$Proteome_Sample_ID # And subset by can
            xsing3 <- meta[which(meta$Sample %in% xsing & meta$cohort == can),]$Proteome_Sample_ID # And subset by can 
            mdf$Grouped = ifelse(mdf$variable %in% sing3 , "M", NA)#This will represent the mutated group (single genes) or cooc. v wt
            mdf$Grouped = ifelse(mdf$variable %in% xsing3 , "W", mdf$Grouped) #This will represent the more wt group
            #mdf$MUTCO = ifelse(mdf$TP53 == "W" & mdf$Grouped == "M", "M", NA) #This grabs ths TP53wt and calls to SINGmut
            #mdf$MUTCO = ifelse(mdf$TP53 == "W" & mdf$Grouped == "W", "W", mdf$MUTCO) #This grabs the TP53wt and SINGwt
            mutco_mdf <- mdf[which(!is.na(mdf$Grouped)),] #This step will remove unwanted cancer types
            tbl <- table(mutco_mdf$Grouped)
            yo = c(tbl[1] > 5,tbl[2] > 5)
            #see https://stackoverflow.com/questions/16822426/dealing-with-true-false-na-and-nan
            if(length(tbl) > 1 & yo[1] %in% TRUE & yo[2] %in% TRUE){ #This makes sure there are enough mutations to make a call trick 
                mutco_triple <- lda(formula=Grouped~value,data=mutco_mdf)
                direction <- ifelse(mutco_triple$means[2] < mutco_triple$means[1],"RespIsSmallerInMut","RespIsLargerInMut")
                out = data.frame(Cancer=can,Genes=g,Sing=s,LDAsvd=mutco_triple$svd,direction=direction)
                SING3 = rbind(SING3,out)
            }
        }
    }
}

topsSING3 <- head(SING3[order(SING3$LDAsvd, decreasing =T),],500)
s3 <- sort(table(topsSING3$Genes))
s4 <- s3[which(s3>1)]
s4n <- names(s4)
s5 <- topsSING3[which(topsSING3$Genes %in% s4n),]
s5[order(s5$Genes),]

islarger <- SING3[which(SING3$direction == "RespIsLargerInMut" & SING3$Sing == "VHL"),]
il <- sort(table(islarger$Genes))
ilg <- data.frame(islarger %>% group_by(Genes) %>% summarize(Yo=sum(LDAsvd)))
ilgs <- ilg[order(ilg$Yo),]



#Some possible FAVs
favsup_TP53 <- c("DDB2","IDNK","FBXO22","NT5C","RABGEF1")
favsup_KRAS <- c(?)
favsup_PIK3CA <- c() 


#Questions, How should we layer on Cancer type? 







#ONCOPROTEINS 
#TP53 HITS 
4462            SNF8 16.700997474
1174           DCAF8 16.721945252
5219           VPS25 16.776185863
4586            SSH3 16.813448142
3654          PRKAB1 16.925512642
2230          IMPACT 16.926428450
4912          TOLLIP 16.945970666
910             CMBL 17.028219298
2052            HID1 17.143211674
296            APPL2 17.226734992
3425          PHYHD1 17.242830225
5214            VMAC 17.286149839
1635             FAS 17.291442219
1778            GAMT 17.292832676
2954            NANS 17.342294228
526           BORCS7 17.410748068
3426          PHYKPL 17.418717047
3176          NUDT12 17.519369853
3832          RAB27B 17.566280606
2763             MPI 17.605190445
1747           FUCA1 17.730474440
3169           NUCB2 17.741445580
4139           RRM2B 17.802215169
2359          KIF13B 17.909091505
1843           GLOD4 17.937865774
1614         FAM172A 18.015313731
2187            IDUA 18.127047115
3023          NDUFS4 18.409066904
4271        SELENBP1 18.446641105
998            COPS4 18.476621722
3421          PHLDA3 18.496670715
4807           THSD4 18.543337976
5229           VPS36 18.548117710
4355            SHPK 18.748429058
473             BAG1 18.829287694
3371           PEBP1 19.663869015
3003          NDUFA7 19.772338093
3776          PTPN23 19.903330369
3025          NDUFS6 20.122871223
1755           FYCO1 20.312078139
4286         SEPSECS 20.437141865
2932           MYO5C 20.762170701
3423           PHPT1 20.777330966
3860         RABGEF1 21.104087691
3856          RABEP1 21.267151804
3156            NT5C 21.607588740
5018           TTC19 22.167995449
1654          FBXO22 23.066819770
1197            DDB2 30.019094207



#KRAS HITS 
1            PDPR 5.6838405159
1391        EPM2AIP1 5.6993745323
3438          PRRC2A 5.7052076620
1553           FKBP8 5.7338277631
1881          HIRIP3 5.7526324599
688            CDK12 5.7581431175
1261           EBAG9 5.7629617251
4575          TOMM34 5.7703711280
4107         SMARCB1 5.7910040462
400           ATXN10 5.8076012974
3584            RAE1 5.8142024901
727            CFDP1 5.8142473288
4277          STEEP1 5.8226467301
4672            TSR2 5.8492640646
3556           RAB29 5.8526339916
3938         SDR39U1 5.8628653759
3434           PRPS1 5.8758708270
680            CDC5L 5.9278243591
7               AAR2 5.9556826937
1382           EPB41 6.0741760308
2869          NMRAL1 6.0919782220
4811           USP13 6.1048173226
3155            PEX3 6.1293384693
1801          GTF3C3 6.1294865960
1196          DNAJC8 6.1311122657
124              AHR 6.1861230839
4335          SUPT6H 6.2418888274
3423          PRPF19 6.3837475540
3685            RFC3 6.3959687608
2310            LSG1 6.4359267464
4920           WDR70 6.4437135360
242         APOBEC3F 6.5104315207
1346            ELOA 6.5147973650
1281           EEF2K 6.6457114885
2194            KRI1 6.6595807610
636             CCNK 6.7003895176
1145           DHX38 6.7883642400
3684            RFC2 6.8684482979
3686            RFC4 6.9025924220
2414            MAU2 6.9094297527
2939          NUDCD2 7.0153496286
821            CMTR1 7.9072700737
1270            EDC4 8.3987942563
3059           PALS2 8.5678880992

#PIK3CA
4364         SELENOH 7.754477e+00
4684           SSBP1 7.773867e+00
2815          MRPL10 7.788609e+00
2823          MRPL19 7.823581e+00
2856           MRPL9 7.824266e+00
1003           COQ8A 7.932347e+00
2081          HIBADH 7.933945e+00
4020           RCC1L 7.987218e+00
2870          MRPS25 7.998698e+00
2814           MRPL1 8.004825e+00
1308          DNAJA3 8.024822e+00
2058           HDHD5 8.049492e+00
5050            TPX2 8.052484e+00
2843          MRPL41 8.071565e+00
2817          MRPL12 8.080321e+00
3155           NOL11 8.128905e+00
2851          MRPL50 8.148827e+00
2853          MRPL53 8.191063e+00
2921           MTIF2 8.215201e+00
5293           UTP25 8.215239e+00
1788      GADD45GIP1 8.217574e+00
4880           TFB2M 8.237088e+00
2682            MDH2 8.264846e+00
125              AGK 8.272315e+00
2845          MRPL44 8.310380e+00
1656         FASTKD2 8.322308e+00
4075            RIDA 8.354990e+00
5384           WDR12 8.387816e+00
4233           RRP7A 8.388959e+00
1168           DARS2 8.395635e+00
1271           DHX30 8.467470e+00
2649           MARS2 8.474837e+00
3058         NDUFAF2 8.479735e+00
1216           DDX28 8.520076e+00
1967           GRWD1 8.566957e+00
3056         NDUFAB1 8.627898e+00
4879           TFB1M 8.687046e+00
1001            COQ5 8.718753e+00
2000         GTPBP10 8.806684e+00
2879           MRPS6 8.860501e+00
2842          MRPL40 8.892301e+00
3149            NOA1 8.914043e+00
3999           RBM28 8.993845e+00
3219           NUBPL 9.106210e+00
2880           MRPS7 9.137583e+00
2459           LARS2 9.168775e+00
1164            DAP3 9.384888e+00
3599           PNPT1 9.408863e+00
1872             GLS 9.769087e+00
2197      HSPE1-MOB4 9.888349e+00
1322         DNAJC19 9.936926e+00
2811            MRM1 1.015438e+01
1548           ERAL1 1.086252e+01



#PTEN 
4587            SSB  8.447609462
1078        CTNNBL1  8.495469196
4146          RRP7A  8.518715461
4994        TRMT61A  8.564771114
12             AATF  8.648064012
1481           EMG1  8.659678484
2633           MED1  8.681961990
2310           KAT8  8.684076475
5177         UTP14A  8.718605302
507            BOP1  8.741530900
3096           NOL6  8.803929755
3688          PRPF6  8.825894859
5281          WDR46  9.007146321
1151         DCAF13  9.014410457
2826           MSH6  9.051354493
1202          DDX52  9.090443909
1627            FBL  9.101955962
5183           UTP3  9.145300430
517            BRD4  9.157690025
3097           NOL7  9.175375160
763           CEBPZ  9.268450954
1712          FOXK1  9.306519077
3109          NOP58  9.323839725
5289          WDR75  9.389293606
3110           NOP9  9.424787073
3331         PDCD11  9.484263107
3094          NOL11  9.541549941
2671           MFN1  9.657781389
2209           IMP4  9.687781051
891           CMSS1  9.711832537
2011         HEATR1  9.866616777
2737      MPHOSPH10  9.904924584
2974            NCL  9.980818279
1180          DDX10  9.980944882
1539           ESF1  9.984963347
3908          RBM19 10.059813974
5279          WDR43 10.060260220
3967          REXO4 10.106538974
5184           UTP4 10.148786475
4993          TRMT6 10.727461444
39             ABT1 10.972067167
3108          NOP56 11.121290423
1189          DDX27 11.218619083
3093          NOL10 11.687752877
3752           PTEN 12.457783018


#EGFR
378             BAG2 4.9487496911
2318           MTFMT 4.9653517862
1792          IGFBP4 4.9781785573
1743         HSP90B1 4.9914464760
3227          RNASE3 5.0388322198
3056             PVR 5.0529825605
1056          DNAJB4 5.0603235502
3433            SDK1 5.0680943862
1143          EEF1B2 5.0822826685
634            CEMIP 5.1178763237
420              BPI 5.1260032452
1860           ITGB1 5.1712741401
3446          SEC24D 5.1776495246
3099           RAB35 5.1960626496
3442         SEC23IP 5.1989403656
3314           RPS12 5.2033453852
4280           WDR43 5.2090552276
1333             FBL 5.2171948337
73              ADA2 5.2777324493
483               C9 5.3581342591
801            COPG1 5.3912646691
979            DDX21 5.4117029862
2161            MESD 5.4460011925
799            COPB2 5.4488673556
2332           MTPAP 5.4992654841
1372             FGR 5.5738251740
1373          FHIP2A 5.5957559190
2086            MANF 5.6245799216
4156            UFM1 5.6613776544
2055            LYAR 5.6789691231
1794          IGFBP6 5.7597916898
1407           FOCAD 5.9583301977
4151            UCK2 5.9706030176
3783           SULF1 5.9774434987
4216           VAT1L 5.9886060466
864         CRISPLD2 6.0132939032
1444           GALK2 6.0333936429
4057           TRUB1 6.0498854314
150          ALDH1L2 6.0698505669
2321          MTHFD2 6.0738663203
2390            NANS 6.0971612446
4378          ZNF622 6.0985170483
897             CTSB 6.1026378294
2760            PFN2 6.1029752922
2899            PPIB 6.1598945323
1920             KIT 6.2091938940
1623           HACL1 6.2795542644
2508            NNMT 6.3560673110
1233            EOGT 6.4075445587
263           ARID3A 6.4580541755
4102          TXNDC5 6.4904606208
364             AVEN 6.5972291793
581            CD209 6.6559192002
1180           EIF3H 6.6980070162
1556            GPC6 7.2682395562
1569            GRB2 7.3704953618


#VHL 
189           ARPC2 3.171344846
3015          VSIG4 3.171621675
589          CTHRC1 3.175207693
1781          NUP85 3.183244438
400            CD99 3.187777571
2153           RALY 3.189448985
711         DNAJB11 3.197914941
938           FKBP9 3.210769912
3006         VPS26B 3.218061738
2209           RCN3 3.237214631
1989          POTEJ 3.250687243
2229          RGS19 3.254147355
1610          MXRA5 3.267202976
298            C1QA 3.282196396
1614         MYCBP2 3.290427593
1652           NAV1 3.295137899
2994           VCAN 3.297366248
1023           GLUL 3.304869682
505           CNPY3 3.311981469
1194         HSPA13 3.336871439
478           CKAP4 3.340040901
1515          MED18 3.347064675
2164          RASA3 3.387762775
386           CD163 3.390060434
192           ARPC5 3.417752376
1119            HCK 3.421640105
300            C1QC 3.443703495
2543         SMCHD1 3.458275685
1564           MNDA 3.519087593
2815         TMSB4X 3.540966565
2823      TNFAIP8L2 3.569163382
1542         MICAL1 3.604116985
1407           LRP1 3.616472124
394            CD44 3.650040293
1041           GNG5 3.672481147
2754          TGFBI 3.696050676
2125          RAB31 3.723908618
476           CILK1 3.809514789
1282          ITGAM 3.825353286
1479         MAPRE1 3.887213008
743            DPYD 3.891900216
2321          RPS10 3.893920420
1469           MAP4 3.922721129




issmaller <- SING3[which(SING3$direction == "RespIsSmallerInMut" & SING3$Sing == "EGFR"),]
is <- sort(table(issmaller$Genes))
isg <- data.frame(issmaller %>% group_by(Genes) %>% summarize(Yo=sum(LDAsvd)))
isgs <- isg[order(isg$Yo),]
FAVS = c("PKM","NOL11",""


#TP53
5411            XPO1 20.695250204
3429            PES1 20.733924349
2107           HMCES 20.751756229
3654           PPM1G 20.762660736
5266           USP39 20.809904411
1177          DCAF13 20.834188059
40              ABT1 20.964335831
867            CIP2A 20.982555801
3010           NCAPG 21.118469903
1424           EEF1G 21.391134217
784            CEBPZ 21.501457604
2748           MKI67 21.555843368
2402           KPNA2 21.621487122
3008          NCAPD2 21.625960769
4225          RSL1D1 21.678194398
1883            GMPS 21.818584419
2669            MDC1 22.079558501
1748           FOXK1 22.095654480
1218           DDX27 22.106084793
391           ATAD3A 22.136411677
5391           WDR75 22.138299103
3383          PDCD11 22.151708744
5008            TOP1 22.202838023
2378           KIF11 22.211695207
764             CDK1 22.280347003
1211           DDX18 22.387376028
4489        SLC4A1AP 22.541204613
2905          MTHFD2 22.784416692
3243           NUP93 23.025690010
2384           KIF2C 23.564323043
5178           UBE2C 23.821946944
5380           WDR43 23.959240102
2056          HEATR1 24.443968668
5073          TRIP13 24.539773914
1081           CSE1L 24.698237356
5413            XPO5 24.859585227 
5281            UTP4 25.236289186 ?
5035            TPX2 26.542790047 ?
3140           NOL11 26.836648449 #Nucleolar integrity during interphase supports faithful Cdk1 activation and mitotic entry, "Depletion of NOL11 delayed entry into the mitotic phase owing to increased inhibitory phosphorylation of cyclin-dependent kinase 1 (Cdk1) and aberrant accumulation of Wee1, a kinase that phosphorylates and inhibits Cdk1."
5009           TOP2A 27.819347774 #HER2 association all over 
4039           REXO4 32.898225212 #A novel nine gene signature integrates stemness characteristics associated with prognosis


#KRAS
3508         RNASE6 5.8129391559
836          COL1A2 5.8167820679
3491           RIN1 5.8594533259
4469          UBE2N 5.8689728861
628           CD209 5.8861633170
1102          DHX32 5.9370672687
1436           FAT1 5.9467999244
4657          XIRP2 5.9609351124
2011          ITGA2 5.9636597353
3664          S100P 5.9719594190
318          ARPC1B 5.9878054813
3407          RAPH1 5.9995487945
4278          TMOD3 6.0136466102
125           AIFM2 6.0485372962
1675           GOT1 6.0642606024
3385          RABL6 6.0699748280
4648          WIPI1 6.0928790042
4501           UGP2 6.1984843797
2823          OLFM4 6.1997062553
3061        PLEKHG2 6.2068728372
4389         TSPAN8 6.2358037309
481            BZW1 6.2427836342
798           CMPK1 6.2598153770
2033       IVNS1ABP 6.2689316026
3654        S100A11 6.3978349632
2104          KRT18 6.4691808607
217           AP2A1 6.4763844166
1141        DNAJC10 6.5120278041
287         ARHGAP4 6.5247144811
1484           FHL2 6.5528729142
1754           H1-5 6.7734176150
4063          STRN4 6.8931660662
2139           LCN2 6.9139686074
1162          DOCK5 7.0060224319
2110           KRT8 7.0374776716 #wrong direction. 
2186         LPGAT1 7.3062749692 #Fifteen mRNAs (ENDOU, MFN2, FASLG, SHOC2, VEGFA, ZFPM2, HOXC6, KLK10, DDIT4, LPGAT1, BEX4, DENND5B, PHF20L1, HSP90B1, and PSPC1) were identified as prognostic biomarkers for CRC by multivariate Cox regression. 
73            ACTR3 7.4739988472
4428           TWF1 7.7531332057
2233            LXN 8.0221689447 #Latexin (Lxn) is a negative regulator of stem cell proliferation and we investigate the effects of Lxn on CD133+ pancreatic cancer stem-like cells.

#PIK3CA
13             ATL3 7.5481943168
1917           GNG12 7.5786549124
237            ANXA1 7.6206640615
533           BORCS7 7.6328729217
1744            FLII 7.7537871555
3854          PTGFRN 7.7732624889
4473           SHTN1 7.7845757875
3138          NIBAN2 7.7864481527
3124           NFKB1 7.8559374864
5184            TWF1 7.8630277877
4458          SH3GL1 7.8838778345
521          BLOC1S4 7.8983077963
296            APPL2 7.9187809484
513             BIN3 7.9217825122
1960             GPI 8.0398157592
216           ANKFY1 8.0630892298
3581           PLIN3 8.1474302768
3532          PIK3R1 8.1577129673
4105           RIPK1 8.2191862122
644            CAPN1 8.2559703564
649           CAPNS1 8.3231735792
3452           PDZD8 8.3262307628
1849          GDPGP1 8.3317239327
1811          GALNT1 8.3466758990
5191         TXNDC17 8.3629637519
5018           TMOD3 8.4371247045
2448           KRT17 8.4871510817
239            ANXA2 8.4925524715
3485            PGK2 8.5858924979
5447         XPNPEP1 8.7729330977
4414        SERPINB6 8.9204081665
2510           LIMA1 8.9858308417
5105           TRIM5 8.9868436225
3702        PPP1R13L 9.1459130670
5265            UGDH 9.2133146155
3677            PPIC 9.3887377566
5518          ZFYVE1 9.6777663969
3549             PKM 9.7076787415 #Wrong direction? 


#PTEN
879           CLINT1 7.669102e+00
4559            ST14 7.671146e+00
5084           UEVLD 7.702216e+00
1048             CRK 7.704589e+00
4920         TRAPPC8 7.734885e+00
3690           PSMC2 7.790275e+00
1678          FHIP2B 7.805198e+00
3906           RCHY1 7.812647e+00
1758             GAK 7.852770e+00
5190          VPS37B 7.916676e+00
976            COPB2 7.919552e+00
2122         HSD17B4 7.950249e+00
398            ATG12 7.963127e+00
3422          PITPNA 7.986451e+00
236            ANXA2 7.994022e+00
5195           VPS4B 8.014121e+00
4722           TCF25 8.047443e+00
3360            PFKP 8.065985e+00
1626             FAS 8.150132e+00
5341         ZFYVE16 8.160232e+00
2902           MYO5C 8.190355e+00
3965           RIOK3 8.466665e+00
3431         PLA2G4A 8.473318e+00
493            BCL10 8.562791e+00
854            CILK1 8.732067e+00
1764          GALNT1 8.739294e+00
4198            SCP2 8.760525e+00
2286         ITPRID2 8.820630e+00
4462            SNX4 8.845307e+00
253            AP3M1 8.863504e+00
2343           KIF1C 8.866745e+00
4252         SEPSECS 8.892823e+00
657            CASP6 8.937647e+00
4915        TRAPPC11 8.940681e+00
3550            PPCS 9.037558e+00
217           ANKMY2 9.095355e+00
213           ANKFY1 9.111313e+00
3428             PKM 9.234985e+00 #There it is again
2340          KIF13B 9.296457e+00
4226          SEC24B 9.314006e+00
4465            SNX7 9.414421e+00
251            AP3B1 9.538609e+00
4322           SHMT1 9.840229e+00
2028          HECTD3 9.931909e+00
1975            GYS1 1.032752e+01 #A number of observations have suggested that the turnover of glycogen is altered in tumor cells. The levels of glycogen were demonstrated to be particularly high in breast, kidney, uterus, bladder, ovary, skin, and brain cancer cell lines. Glycogen content in these cells was inversely correlated with proliferation rate (Rousset et al., 1981).


#VHL 
1683         PLEKHA6 2.9405325694
979          HEATR5B 2.9607098753
2503           ZMYM1 2.9612088916
1215            MAP7 2.9747029575
529             CRYZ 2.9762571122
1118           KNTC1 3.0373675396
274            BPNT1 3.0396814814
1892           RDH10 3.0644983304
617            DIP2C 3.0930830201
1175          LRATD2 3.1046425849
1922           RMND1 3.1122389301
1619           PEBP1 3.1206284022
441            CNDP2 3.1308504016
2353             TXN 3.1744327241
551             CUL9 3.1827450120
96            AKR1C3 3.1911625036
584            DDAH1 3.2072840039
1236            MDH1 3.2119957543
1161           LIN7A 3.2348390237
1225           MAT2B 3.2568738064
905           GNPDA1 3.3656186567
1015           HOOK1 3.4209643444
325             CBR1 3.5461961587
1224           MAT2A 3.6739426712
207            ATOX1 3.7777442766
2087           SMYD3 3.7812792606
197           ASRGL1 3.8368461851
184            ARMT1 3.8635266223
743           EPS8L2 3.9563073652

#EGFR 
1500           HSPB1 5.3518756834
1347           GSPT1 5.3869537402
3659        VKORC1L1 5.3877191618
1720          LRATD2 5.3953869438
1254            GET3 5.4758126697
1331          GPR89A 5.5795186111
3713            WWC1 5.5967024515
318             ATRX 5.6101178081
2692          RABEP2 5.6547094267
2162           NSUN5 5.6588610596
1650           KMT2C 5.7437994627
1265            GHDC 5.8191695636
1458           HMOX2 5.8931003398
1120         FAM120C 5.9082479390
2126          NMRAL1 5.9532357393
1713           LMTK2 5.9821553538
1160            FCSK 5.9960713398
359            BLVRA 6.0388878804
1107            FAAH 6.0481143660
2798           RINT1 6.0504540346
2881          RSPRY1 6.0936421975
2110        NIPSNAP2 6.2137054093
315            ATP9A 6.4164679800
2016            MYO6 6.7896184531
2373          PHLDA3 7.0914392240
1425            HIP1 7.2735518598
1671          LANCL2 7.3711983624
999             EGFR 7.5125454567


#Good paper maybe 


favs = c("ILK","COPB2","FOXK1","HSPA5")
c("DDB2","IDNK","FBXO22","NT5C","RABGEF1")

g = "PKM"
mut = "PTEN"

sing = unique(maf[which(maf$Hugo_Symbol == mut),]$Tumor_Sample_Barcode)
sing3 <- meta[which(meta$Sample %in% sing),]$Proteome_Sample_ID # " 

df <- prot[which(prot$external_gene_name == g),]
df$external_gene_name <- NULL
mdf <- melt(df)
mdf$Tissue = ifelse(mdf$variable %in% unique(meta$Proteome_Normal_Sample_ID),"N","T")
mdf$Mut = ifelse(mdf$variable %in% sing3,"MT","WT")
this_triple <- lda(formula=Tissue~value,data=mdf)
this_mut <- lda(formula=Mut~value,data=mdf)
plot(this_triple)
plot(this_mut)



########################################################KRAS#############################################################
#Step 1: Find the samples with somatic mutations in TP53 
#NOTE: I should/can change this to only include Collins/Eduards/Amelia's predictions
allsamps <- unique(maf$Tumor_Sample_Barcode)
tp53samps <- unique(maf[which(maf$Hugo_Symbol == "KRAS"),]$Tumor_Sample_Barcode)
xtp53samps <- setdiff(allsamps,tp53samps)

#Step 2: Split the Proteomics data up into groups and look for the most differentiated
#NOTE: This can be done a couple different ways. I'm going to start with LDA 
#NOTE: Also try to build some cross validation sets in this. 70 20 10 on this front. 
#Make sure that the column label match between proteins and mutations..
#All _ to . and all - to _ 

#TODO: I need to change this code to pull in the Meta data idenfiers instead of guessing. 
allsamps3 <- meta[which(meta$Sample %in% allsamps),]$Proteome_Sample_ID #Replace -
tp53samps3 <- meta[which(meta$Sample %in% tp53samps),]$Proteome_Sample_ID # " 
xtp53samps3 <- meta[which(meta$Sample %in% xtp53samps),]$Proteome_Sample_ID # " 

#Change the column names in prot 
genes <- prot$external_gene_name

#This just goes through and compares TP53 and TP53wt 
g = "NDUFS6"
look <- NULL
for(g in genes){
    df <- prot[which(prot$external_gene_name == g),]
    df$external_gene_name <- NULL
    mdf <- melt(df)
    mdf$Grouped = ifelse(mdf$variable %in% tp53samps3, "M", NA)#This will represent the mutated group (single genes) or cooc. v wt
    mdf$Grouped = ifelse(mdf$variable %in% xtp53samps3, "W", mdf$Grouped) #This will represent the more wt group 
    clean_mdf <- mdf[which(!is.na(mdf$Grouped)),]
    this_triple <- lda(formula=Grouped~value,data=clean_mdf)
    direction <- ifelse(this_triple$means[2] < this_triple$means[1],"RepIsLarger","RespIsSmaller")
    out = data.frame(Genes=g,LDAsvd=this_triple$svd,direction=direction)
    look = rbind(look,out)
}

#This not takes the top hit from TP53 and looks at which looks best from T/N comparisons
topsTP53 <- head(look[order(look$LDAsvd, decreasing =T),],100)
tumor <- na.omit(unique(meta$Proteome_Sample_ID))
normal <- na.omit(unique(meta$Proteome_Normal_Sample_ID))

lookTN <- NULL
for(g in topsTP53$Genes){
    df <- prot[which(prot$external_gene_name == g),]
    tG <- topsTP53[which(topsTP53$Genes == g),]
    df$external_gene_name <- NULL
    mdf <- melt(df)
    mdf$Tissue = ifelse(mdf$variable %in% normal,"N","T")
    this_triple <- lda(formula=Tissue~value,data=mdf)
    direction <- ifelse(this_triple$means[2] < this_triple$means[1],"RespIsSmall","RespIsLarger")
    out = data.frame(Genes=g,LDAgene=tG$LDAsvd,LDAsvd=this_triple$svd,dirMUT=tG$direction,dirNAT=direction)
    lookTN = rbind(lookTN,out)
}

lookTN$LDAsum = lookTN$LDAgene+lookTN$LDAsvd

#Does this makes sense with a plot: 
#Think throught he co-occuring events 
ordLDA <- head(lookTN[order(lookTN$LDAsum, decreasing =T),],100)

#ID by Cohort
tcode <- na.omit(data.frame(meta %>% select(Proteome_Sample_ID,cohort)))
ncode <- na.omit(data.frame(meta %>% select(Proteome_Normal_Sample_ID,cohort)))
ncode <- ncode[which(ncode$Proteome_Normal_Sample_ID != ""),]
colnames(tcode) <- c("ID","cohort")
colnames(ncode) <- c("ID","cohort")
codes <- rbind(tcode,ncode)


favs3 <- c("DHRS11") #KRAS hits 
#Cell proliferation and met in TNBC NUP93 (2020)
#Promotes migration (CSE1L) 
#DNA damage: Topoisomerases (TOP1)

g = "DHRS11"
g = "SYNMA"
g = "PRSS1"
g = "RPL5"
g = "TMEM214"

#for(g in favs){
    df <- prot[which(prot$external_gene_name == g),]
    df$external_gene_name <- NULL
    mdf <- melt(df)
    mdf$TN = ifelse(mdf$variable %in% normal,"N","T")
    mdf$MW = ifelse(mdf$variable %in% tp53samps3,"M","W")
    mdf$Group = paste(mdf$TN,mdf$MW,sep="_")
    mdf2 <- merge(mdf,codes,by.x="variable",by.y="ID")
#}

p <- ggplot(mdf2,aes(x=Group,y=value,fill=cohort))
p <- p + geom_boxplot()
p <- p + theme_classic()
p <- p + ggtitle(g)
p




#LETs look at this a bit differently  For what LDA was meant to do. 
yo <- prot
myo <- melt(prot)
dyo <- dcast(myo,variable~external_gene_name)
#dyo$TN <- ifelse(dyo$variable %in% normal,"N","T")
dyo$MW <- ifelse(dyo$variable %in% tp53samps3,"M","W")
dyoc <- merge(dyo,codes,by.x="variable",by.y="ID")
rownames(dyoc) <- dyoc$variable
dyoc$variable <- NULL 

bigLDA <- lda(MW ~ ., data = dyoc)



#Now can I think about a way to do the cisTrans Wilcox T

















 


