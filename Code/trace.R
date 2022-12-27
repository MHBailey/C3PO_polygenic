library(data.table)
library(stringr)
library(reshape2)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyr)

'%!in%' <- function(x,y)!('%in%'(x,y))
args = commandArgs(trailingOnly=TRUE)

#args = c("Data/Hallmarks/nanostring.gl.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/NESTs/TheNEST.csv","Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Somatic_cnv/Broad_pipeline_wxs_formatted/GISTIC_all_thresholded.by_genes.tsv","Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_raw_gene_level_v1.3.tsv","Processed_data/BRCA.PolyRisk.dna.v2.txt","Processed_data/BRCA.PolyRisk.amp.v2.txt","Processed_data/BRCA.PolyRisk.del.v2.txt","Processed_data/BRCA.cnv.p.effect.mannu.amplification.txt","Processed_data/BRCA.cnv.p.effect.mannu.deletion.txt","Processed_data/BRCA.dna.p.effect.mannu.txt","X01BR043","Adaptive_immunity")

args = c("Data/Hallmarks/nanostring.gl.txt","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","Data/NESTs/TheNEST.csv","Data/Somatic_mutation_wxs/Broad_WashU_union/PanCan_Union_Maf_Broad_WashU_v1.1.maf","Data/Somatic_cnv/Broad_pipeline_wxs_formatted/GISTIC_all_thresholded.by_genes.tsv","Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_raw_gene_level_v1.3.tsv","Processed_data/PDAC.PolyRisk.dna.v2.txt","Processed_data/PDAC.PolyRisk.amp.v2.txt","Processed_data/PDAC.PolyRisk.del.v2.txt","Processed_data/PDAC.cnv.p.effect.mannu.amplification.txt","Processed_data/PDAC.cnv.p.effect.mannu.deletion.txt","Processed_data/PDAC.dna.p.effect.mannu.txt","X01BR043","Adaptive_immunity","KRAS","GZMK","PDAC")


#Get all of the genes in the hallmark
halls = fread(args[1])

samp = args[13]
myhall = c("Genes",args[14])

myh = halls %>% select(all_of(myhall))
genesInHall = myh[which(myh[,2] == "Y"),]$Genes


#Get all of the sample ids from meta
meta = fread(args[2])
myids = meta[which(meta$Proteome_Sample_ID == samp),]
mutid = myids$WXS
cnvid = myids$CASE_ID
protid = samp


#Get all of the genes in a nest
nest = fread(args[3],sep=",",header=F)

#Get the maf
keep_muts <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins")
maf <- fread(args[4])
maf_cnt <- maf[which(maf$Variant_Classification %in% keep_muts),]
sort(table(maf_cnt$Tumor_Sample_Barcode))
submaf <- maf[which(maf$Tumor_Sample_Barcode == mutid),]
smafm <- submaf[which(submaf$Variant_Classification %in% keep_muts),] 

#Get the CNV
cnv <- fread(args[5])
cnv_keep = c("Gene Symbol",cnvid)
subcnv <- cnv %>% select(all_of(cnv_keep))

#Get Proteins in the list
prot <- fread(args[6])
protgenes <- prot$external_gene_name

#Intersect hallmarks with Proteins
all_genes <- intersect(protgenes,genesInHall)

#What are the best nests for all_genes DNA
dna_eff <- fread(args[12]) 
dna_seff <- dna_eff[which(dna_eff$Protein %in% all_genes),]
dna_good_nest <- dna_seff[which(dna_seff$TPVALUE < 0.1),]

#What are the best nests for cnva
amp_eff <- fread(args[10])
amp_seff <- amp_eff[which(amp_eff$Protein %in% all_genes),]
amp_good_nest <- amp_seff[which(amp_seff$TPVALUE < 0.1),] 

#What are the best nests for cnvd
del_eff <- fread(args[11])
del_seff <- del_eff[which(del_eff$Protein %in% all_genes),]
del_good_nest <- del_seff[which(del_seff$TPVALUE < 0.1),] 

#So now we are down to 3000+ Nests for 26 proteins
#dna_good_nest[order(dna_good_nest$COHEN_D),]
#amp_good_nest[order(amp_good_nest$COHEN_D),]
#del_good_nest[order(del_good_nest$COHEN_D),]



#Skip a couple steps and go strait to scores 
to_keep = c("Proteome_Sample_ID",all_genes)
prsdna <- fread(args[7])
prsdna_sub <- prsdna %>% select(all_of(to_keep))
prsdna_sub[which(prsdna_sub$Proteome_Sample_ID == samp),]

prsamp <- fread(args[8])
prsamp_sub <- prsamp %>% select(all_of(to_keep))
prsamp_sub[which(prsamp_sub$Proteome_Sample_ID == samp),]

prsdel <- fread(args[9])
prsdel_sub <- prsdel %>% select(all_of(to_keep))
prsdel_sub[which(prsdel_sub$Proteome_Sample_ID == samp),]


#Which mutations are contributing Protein score
DNA_contribs = NULL
PROTEIN_DNA_contribs = NULL
for(g in all_genes){
    yo = dna_good_nest[which(dna_good_nest$Protein == g),]
    for(n in yo$NESTv1){
        yoGenes <- str_split(nest[which(nest$V1 == n),]$V3," ")[[1]]
        yo2 = yo[which(yo$NESTv1 == n),]
        sampmuts <- unique(smafm[which(smafm$Hugo_Symbol %in% yoGenes),]$Hugo_Symbol)
        DNA_contribs = c(sampmuts,DNA_contribs)
        if(length(sampmuts) == 0){
            hit_collapse = NA 
        }else{
            hit_collapse = paste(unique(sampmuts),collapse=" ")
        } 
        out = data.frame(Sample=samp,PROTEIN=g,NEST=n,HITs=hit_collapse, EFFECT=yo2$COHEN_D, PVALUE=yo2$TPVALUE,yo2$nMUT, yo2$nWT, nGENES = length(yoGenes) )
        PROTEIN_DNA_contribs = rbind(out,PROTEIN_DNA_contribs)
    }
}



#Now I'm going to take a slightly different approach. 
#Count the shared nest mutations for all immune hot subtypes 

immune_hot = c("X01BR008", "X01BR017", "X01BR018", "X01BR026", "X01BR030", "X01BR040","X01BR043", "X05BR001", "X11BR003", "X11BR004", "X11BR011", "X11BR016","X11BR017", "X11BR018", "X11BR020", "X11BR023", "X11BR024", "X11BR030","X11BR038", "X11BR043", "X11BR073", "X11BR080", "X18BR003", "X18BR009","X21BR001")

immune_cooler = c("CPT000814", "CPT001846", "X01BR001", "X01BR010", "X01BR015", "X01BR023", "X01BR025", "X01BR027", "X01BR031", "X01BR032", "X03BR002", "X03BR004", "X03BR005", "X03BR006", "X03BR010", "X05BR003", "X05BR004", "X05BR005", "X05BR016", "X05BR038", "X05BR044", "X05BR045", "X06BR014", "X09BR001", "X09BR004", "X09BR005", "X11BR006", "X11BR012", "X11BR013", "X11BR015", "X11BR019", "X11BR022", "X11BR027", "X11BR028", "X11BR036", "X11BR040", "X11BR042", "X11BR047", "X11BR051", "X11BR053", "X11BR054", "X11BR055", "X11BR058", "X11BR060", "X11BR074", "X11BR075", "X13BR009", "X14BR005", "X14BR008", "X14BR014", "X16BR012", "X18BR004", "X18BR006", "X18BR010", "X18BR016", "X18BR019", "X20BR001", "X20BR002", "X20BR005", "X20BR006", "X20BR008", "X21BR010")

myhots = meta[which(meta$Proteome_Sample_ID %in% immune_hot),]
muthots = myhots$WXS
cnvhots = myhots$CASE_ID
prothots = immune_hot

mycooler = meta[which(meta$Proteome_Sample_ID %in% immune_cooler),]
mutcooler = mycooler$WXS
cnvcooler = mycooler$CASE_ID
protcooler = immune_cooler


samples = maf_cnt %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  distinct() %>%
  rename(Genes = Hugo_Symbol) %>%
  print()

colnames(nest) <- c("NEST", "Description", "Genes")
nests = nest %>%
  select(NEST, Genes) %>%
  separate_rows(Genes, sep = " ") %>% 
  print()

matchMatrix = full_join(samples, nests, by="Genes") %>%
  ungroup() %>%
  #print()
  drop_na(NEST) %>%
  drop_na(Tumor_Sample_Barcode) %>%
  select(-Genes) %>%
  distinct() %>%
  mutate(Value = 1) %>%
  pivot_wider(names_from = NEST, values_from = Value, names_sort = TRUE, values_fill = 0) %>%
  print()

sub_mm <- matchMatrix[which(matchMatrix$Tumor_Sample_Barcode %in% muthots),]
sort(colSums(sub_mm[-1]))





#TOP HOT (21 of 25 have a mutation in this nest) 

1  NEST:91    BRCA1
2  NEST:91    BUB1B
3  NEST:91    CCNA2
4  NEST:91    CCNB1
5  NEST:91    CCND1
6  NEST:91    CCNE1
7  NEST:91     CCNH
8  NEST:91     CDC6
9  NEST:91     CDK1
10 NEST:91     CDK2
11 NEST:91   CDKN1A
12 NEST:91   CDKN1B
13 NEST:91   CDKN2A
14 NEST:91     CDT1
15 NEST:91    CKS1B
16 NEST:91    COPS5
17 NEST:91   CREBBP
18 NEST:91   CTNNB1
19 NEST:91     CUL1
20 NEST:91      DTL
21 NEST:91     E2F1
22 NEST:91    EP300
23 NEST:91     FEN1
24 NEST:91     FZR1
25 NEST:91    HDAC1
26 NEST:91 HIST1H3A
27 NEST:91    HLA-A
28 NEST:91     HUS1
29 NEST:91    KAT2B
30 NEST:91     MCM4
31 NEST:91    MYBL2
32 NEST:91      MYC
33 NEST:91     PCNA
34 NEST:91    RAD17
35 NEST:91      RB1
36 NEST:91     RBL1
37 NEST:91     RRM2
38 NEST:91     SKP2
39 NEST:91    TOP2A
40 NEST:91     TP53


#Now I just want to get some plots 
#Go after pairs
qgene = args[15]
qprot = args[16]
qtype = args[17]

canmaf <- maf_cnt[which(maf_cnt$COHORT == qtype),]
allsamps <- unique(canmaf$Tumor_Sample_Barcode)
mutsamps <- unique(canmaf[which(canmaf$Hugo_Symbol == qgene),]$Tumor_Sample_Barcode)
wtsamps <- setdiff(allsamps,mutsamps)

#Now I need to convert the dna to protids 
protmuts <- meta[which(meta$WXS %in% mutsamps),]$Proteome_Sample_ID
protwts <- meta[which(meta$WXS %in% wtsamps),]$Proteome_Sample_ID
ctmeta <- meta[which(meta$cohort == qtype),]
normals <- ctmeta$Proteome_Normal_Sample_ID[ctmeta$Proteome_Normal_Sample_ID != ""]

#Now I need to 
prott <- prot
GN <- prot$external_gene_name
prott$external_gene_name = NULL
prottt <- data.frame(t(prott))
colnames(prottt) <- GN
prottt$PSAMP <- rownames(prottt)
prottt$Status <- ifelse(rownames(prottt) %in% normals,"Normal",NA)
prottt$Status <- ifelse(rownames(prottt) %in% protmuts,"NEST_mut",prottt$Status)
prottt$Status <- ifelse(rownames(prottt) %in% protwts,"NEST_wt",prottt$Status)

#Grab all GZM genes 
gzm <- c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM")
plotp <- prottt %>% select(all_of(qprot),Status,all_of(gzm))
plotp$SAMP = rownames(plotp)
mp <- melt(plotp)

p <- ggplot(mp,aes(x=Status,y=value))
p <- p + geom_boxplot()
p <- p + ggtitle(paste(qgene,qprot,sep="  "))
p <- p + facet_wrap(.~variable)
p

