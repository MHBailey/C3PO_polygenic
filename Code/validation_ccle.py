def makeGeneList(genespace):
    genes = genespace.split(" ")
    return genes


def getNestDict(fname):
    NEST = {}
    with open(fname,"r") as f:
        for line in f:
            nest,description,genes = line[:-1].split(",")
            Genes = makeGeneList(genes)
            NEST[nest] = Genes
    f.close()
    return NEST


def getMAFadj(fname,isDriver):
    #make a pandas data frame of the MAF
    #return a high impact WT/Mut 0,1 table for all genes
    #This will call all drivers as 
    #With high Impact, SIFT, Polyphen mutations (missense)
    #All frameshift and nonsense
    #NOTE: Bring in VAF differences
    maf = pandas.read_csv(fname,low_memory=False)
    if(isDriver):
        maf_sub = maf
    else:
        smaf = maf[(maf['Variant_Classification'] == "Frame_Shift_Del") | \
        (maf['Variant_Classification'] == "Frame_Shift_Ins") | \
        (maf['Variant_Classification'] == "Missense_Mutation") | \
        (maf['Variant_Classification'] == "Nonsense_Mutation")] #all of these conditions 
        maf_sub = smaf[["Hugo_Symbol","Chromosome","Start_position","Reference_Allele","DepMap_ID","Annotation_Transcript","Codon_Change","Variant_annotation"]]
    return(maf_sub)

def makeDELMatrix(cnv,prots,nests,can):
    if can != "PANCAN":
        ccnv = cnv[cnv["cohort"] == can]
    else:
        ccnv = cnv
    allsamps = ccnv["DepMap_ID"].unique().tolist()
    SAMPxNEST = pandas.DataFrame()
    for i in nests:
        sgenes = nests[i] #This gets all genes in a nest 
        cnv2 = ccnv[ccnv['Gene'].isin(sgenes)]
        nAMP = len(cnv2[cnv2['value'] > 0])
        amp = cnv2[cnv2['value'] > 0]
        uamp = amp['DepMap_ID'].unique().tolist() #Samples with amp
        nDEL = len(cnv2[cnv2['value'] < 0])
        dell = cnv2[cnv2['value'] < 0]
        udel = dell['DepMap_ID'].unique().tolist()
        nWT = len(cnv2[cnv2['value'] == 0])
        cnv2test = cnv2[cnv2['DepMap_ID'].isin(udel)] #BUILD AND LOOP TO DO THIS FOR DELETIONS AND AMP Seperately
        cnvsamps = cnv2test["DepMap_ID"].unique().tolist()
        prots2 = prots[prots.index.isin(allsamps)]
        mutated = prots2[prots2.index.isin(cnvsamps)]
        wt = prots2[~prots2.index.isin(cnvsamps)]
        yowt = wt.index.tolist()
        wt0 = [0]*len(yowt)
        nestid_wt = [i]*len(yowt)
        yomut = mutated.index.tolist()
        mut1 = [1]*len(yomut)
        nestid_mut = [i]*len(yomut)
        wt_tuples = list(zip(yowt, nestid_wt, wt0))
        mut_tuples = list(zip(yomut, nestid_mut, mut1))
        df_wt = pandas.DataFrame(wt_tuples, columns = ['DepMap_ID', 'NEST', "Mut_WT"])
        df_mut = pandas.DataFrame(mut_tuples, columns = ['DepMap_ID', 'NEST', "Mut_WT"])
        mydf = pandas.concat([df_wt,df_mut])
        yuck = mydf.pivot_table(index='DepMap_ID',columns='NEST',values='Mut_WT')
        if SAMPxNEST.empty:
             SAMPxNEST = yuck
        else:
             SAMPxNEST = pandas.merge(SAMPxNEST,yuck,left_index=True, right_index=True, how='outer').fillna(0)
    return(SAMPxNEST)


def makeAMPMatrix(cnv,prots,nests,can):
    if can != "PANCAN":
        ccnv = cnv[cnv["cohort"] == can]
    else:
        ccnv = cnv
    allsamps = ccnv["DepMap_ID"].unique().tolist()
    SAMPxNEST = pandas.DataFrame()
    for i in nests:
        sgenes = nests[i] #This gets all genes in a nest 
        cnv2 = ccnv[ccnv['Gene'].isin(sgenes)]
        nAMP = len(cnv2[cnv2['value'] > 0])
        amp = cnv2[cnv2['value'] > 0]
        uamp = amp['DepMap_ID'].unique().tolist() #Samples with amp
        nDEL = len(cnv2[cnv2['value'] < 0])
        dell = cnv2[cnv2['value'] < 0]
        udel = dell['DepMap_ID'].unique().tolist()
        nWT = len(cnv2[cnv2['value'] == 0])
        cnv2test = cnv2[cnv2['DepMap_ID'].isin(uamp)] #BUILD AND LOOP TO DO THIS FOR DELETIONS AND AMP Seperately
        cnvsamps = cnv2test["DepMap_ID"].unique().tolist()
        prots2 = prots[prots.index.isin(allsamps)]
        mutated = prots2[prots2.index.isin(cnvsamps)]
        wt = prots2[~prots2.index.isin(cnvsamps)]
        yowt = wt.index.tolist()
        wt0 = [0]*len(yowt)
        nestid_wt = [i]*len(yowt)
        yomut = mutated.index.tolist()
        mut1 = [1]*len(yomut)
        nestid_mut = [i]*len(yomut)
        wt_tuples = list(zip(yowt, nestid_wt, wt0))
        mut_tuples = list(zip(yomut, nestid_mut, mut1))
        df_wt = pandas.DataFrame(wt_tuples, columns = ['DepMap_ID', 'NEST', "Mut_WT"])
        df_mut = pandas.DataFrame(mut_tuples, columns = ['DepMap_ID', 'NEST', "Mut_WT"])
        mydf = pandas.concat([df_wt,df_mut])
        yuck = mydf.pivot_table(index='DepMap_ID',columns='NEST',values='Mut_WT')
        if SAMPxNEST.empty:
             SAMPxNEST = yuck
        else:
             SAMPxNEST = pandas.merge(SAMPxNEST,yuck,left_index=True, right_index=True, how='outer').fillna(0)
    return(SAMPxNEST)


def makeDNAMatrix(muts,prots,nests,can):
    if can == "PANCAN":
        cmuts = muts
    elif can == "COAD":
        cmuts = muts[(muts["CODE"] == "COAD") | (muts["CODE"] == "READ")]
    else:
        cmuts = muts[muts["CODE"] == can]
    allsamps = cmuts["DepMap_ID"].unique().tolist() #can specific
    SAMPxNEST = pandas.DataFrame()
    for i in nests:
        sgenes = nests[i] #This gets all genes in a nest 
        muts2 = cmuts[cmuts['Hugo_Symbol'].isin(sgenes)]
        mutsamps = muts2["DepMap_ID"].unique().tolist() #Samples with mut in any gene in nest
        #Subset PROTS
        prots2 = prots[prots.index.isin(allsamps)]
        mutated = prots2[prots2.index.isin(mutsamps)]
        wt = prots2[~prots2.index.isin(mutsamps)]
        columns = prots2.columns.tolist()
        yowt = wt.index.tolist()
        wt0 = [0]*len(yowt)
        nestid_wt = [i]*len(yowt)
        yomut = mutated.index.tolist()
        mut1 = [1]*len(yomut)
        nestid_mut = [i]*len(yomut)
        wt_tuples = list(zip(yowt, nestid_wt, wt0))
        mut_tuples = list(zip(yomut, nestid_mut, mut1))
        df_wt = pandas.DataFrame(wt_tuples, columns = ['DepMap_ID', 'NEST', "Mut_WT"])
        df_mut = pandas.DataFrame(mut_tuples, columns = ['DepMap_ID', 'NEST', "Mut_WT"])
        mydf = pandas.concat([df_wt,df_mut])
        yuck = mydf.pivot_table(index='DepMap_ID',columns='NEST',values='Mut_WT')
        if SAMPxNEST.empty:
             SAMPxNEST = yuck
        else:
             SAMPxNEST = pandas.merge(SAMPxNEST,yuck,left_index=True, right_index=True,how='outer').fillna(0)
    return(SAMPxNEST)



if __name__ == "__main__":
    import sys
    import numpy 
    import pandas
    from scipy import stats

    wdna = pandas.read_csv(sys.argv[1],sep="\t")
    wcnva = pandas.read_csv(sys.argv[2],sep="\t")
    wcnva = pandas.read_csv(sys.argv[2],sep="\t")
    wcnvd = pandas.read_csv(sys.argv[3],sep="\t")
    pandna = pandas.read_csv(sys.argv[4],sep="\t")
    pancnva = pandas.read_csv(sys.argv[5],sep="\t")
    pancnvd = pandas.read_csv(sys.argv[6],sep="\t")
    ccled = getMAFadj(sys.argv[7],False)
    cclec = pandas.read_csv(sys.argv[8],sep="\t") 
    cclep = pandas.read_csv(sys.argv[9])
    nest = getNestDict(sys.argv[10])
    odna = sys.argv[11] 
    oamp = sys.argv[12]
    odel = sys.argv[13]
    cancer = sys.argv[14]
    pdna = sys.argv[15]
    pcnva = sys.argv[16]
    pcnvd = sys.argv[17]
    sampinfo = pandas.read_csv(sys.argv[18])
    oprotts = sys.argv[19]

    #Now I want to subset the cell lines 
    subset_can = []
    if cancer == "OV":
        subset_can = sampinfo[sampinfo['primary_disease'] == 'Ovarian Cancer']
    elif cancer == "BRCA":
        subset_can = sampinfo[sampinfo['primary_disease'] == 'Breast Cancer']
    elif cancer == "COAD":
        subset_can = sampinfo[sampinfo['primary_disease'] == 'Colon/Colorectal Cancer']
    else:
        subset_can = sampinfo 


    #CLEAN up Protein 
    todrop = ["Protein_Id","Description","Group_ID","Uniprot","Uniprot_Acc","TenPx01_Peptides","TenPx02_Peptides","TenPx03_Peptides","TenPx04_Peptides","TenPx06_Peptides","TenPx07_Peptides","TenPx08_Peptides","TenPx09_Peptides","TenPx10_Peptides","TenPx11_Peptides","TenPx12_Peptides","TenPx13_Peptides","TenPx16_Peptides","TenPx17_Peptides","TenPx19_Peptides","TenPx20_Peptides","TenPx21_Peptides","TenPx22_Peptides","TenPx27_Peptides","TenPx28_Peptides","TenPx29_Peptides","TenPx33_Peptides","TenPx34_Peptides","TenPx35_Peptides","TenPx36_Peptides","TenPx37_Peptides","TenPx38_Peptides","TenPx39_Peptides","TenPx40_Peptides","TenPx42_Peptides","TenPx23_Peptides","TenPx05_Peptides","TenPx30_Peptides","TenPx31_Peptides","TenPx32_Peptides","TenPx14_Peptides","TenPx15_Peptides","TenPx41_Peptides","TenPx26_Peptides","TenPx25_Peptides","TenPx18_Peptides","TenPx24_Peptides"]

    cclep2 = cclep.drop(todrop,axis=1)
    yo = cclep2.columns.to_list()
    yo2 = [i.rsplit("_",1)[0] for i in yo ]
    yo3 = [value for value in yo2 if value in sampinfo['CCLE_Name'].to_list()]
    dmid = dict(zip(sampinfo.CCLE_Name,sampinfo.DepMap_ID))
    depmapid = ['Gene'] + [dmid[i] for i in yo2 if i in dmid.keys()]
    cclep2.set_axis(depmapid, axis=1, inplace=True)
    prott = cclep2.set_index('Gene').T.rename_axis('Sample')
    protts = prott[prott.index.isin(subset_can.DepMap_ID)]
    protts.to_csv(oprotts,sep="\t")

    #Clean up CNV
    mcnv = pandas.melt(cclec.reset_index(), id_vars='index')
    mcnv.set_axis(['DepMap_ID','Gene','value'], axis=1, inplace=True)
    mcnv_sub = mcnv[mcnv.DepMap_ID.isin(protts.index.to_list())]
    mcnv_sub['cohort'] = cancer

    #Clean up DNA (subset =to only data that we have protein on
    ccled2 = ccled[ccled.DepMap_ID.isin(protts.index.to_list())]
    ccled2['CODE'] = cancer

    #Make DNA Matrix
    DNA_MAT = makeDNAMatrix(ccled2,protts,nest,cancer)
    #print(DNA_MAT)

    #MAKE CNVa MATRIX
    AMP_MAT = makeAMPMatrix(mcnv_sub,protts,nest,cancer)
    #print(AMP_MAT)

    #DEL MAT
    DEL_MAT = makeDELMatrix(mcnv_sub,protts,nest,cancer)
    #print(DEL_MAT)

    #NOTE: I may be able to squeak these data a bit tighter for higher 
    #R2 values if I raise this value SET A PVALUE Threshold
    thresh = 0.1

    #THIS IS FOR USING THE CIS.CANCER WEIGHTS
    #DNA weights
    sub_wdna = wdna[wdna['TPVALUE'] < thresh]
    pwdna = sub_wdna.pivot_table(index='Protein',columns='NESTv1',values='COHEN_D').transpose()
    DNA_MATclean = DNA_MAT.reindex(columns=pwdna.index.tolist())
    dna_result = DNA_MATclean.dot(pwdna.fillna(0))

    #AMP weights
    sub_wcnva = wcnva[wcnva['TPVALUE'] < thresh]
    pwcnva = sub_wcnva.pivot_table(index='Protein',columns='NESTv1',values='COHEN_D').transpose()
    AMP_MATclean = AMP_MAT.reindex(columns=pwcnva.index.tolist())
    amp_result = AMP_MATclean.dot(pwcnva.fillna(0))

    #DEL weights
    sub_wcnvd = wcnvd[wcnvd['TPVALUE'] < thresh]
    pwcnvd = sub_wcnvd.pivot_table(index='Protein',columns='NESTv1',values='COHEN_D').transpose()
    DEL_MATclean = DEL_MAT.reindex(columns=pwcnvd.index.tolist())
    del_result = DEL_MATclean.dot(pwcnvd.fillna(0))

    #THIS IF FOR USING THE PANCAN WEIGHTS
    #DNA weights
    sub_pandna = pandna[pandna['TPVALUE'] < thresh]
    ppandna = sub_pandna.pivot_table(index='Protein',columns='NESTv1',values='COHEN_D').transpose()
    DNA_MATclean = DNA_MAT.reindex(columns=ppandna.index.tolist())
    pan_dna_result = DNA_MATclean.dot(ppandna.fillna(0))

    #AMP weights
    sub_pancnva = pancnva[pancnva['TPVALUE'] < thresh]
    ppancnva = sub_pancnva.pivot_table(index='Protein',columns='NESTv1',values='COHEN_D').transpose()
    AMP_MATclean = AMP_MAT.reindex(columns=ppancnva.index.tolist())
    pan_amp_result = AMP_MATclean.dot(ppancnva.fillna(0))

    #DEL weights
    sub_pancnvd = pancnvd[pancnvd['TPVALUE'] < thresh]
    ppancnvd = sub_pancnvd.pivot_table(index='Protein',columns='NESTv1',values='COHEN_D').transpose()
    DEL_MATclean = DEL_MAT.reindex(columns=ppancnvd.index.tolist())
    pan_del_result = DEL_MATclean.dot(ppancnvd.fillna(0)) 
   

    dna_result.to_csv(odna, sep="\t")
    amp_result.to_csv(oamp, sep="\t")
    del_result.to_csv(odel, sep="\t")
    pan_dna_result.to_csv(pdna, sep="\t")
    pan_amp_result.to_csv(pcnva, sep="\t")
    pan_del_result.to_csv(pcnvd, sep="\t")

#This is for testing purposes: 
#argv = ["CODE","Processed_data/COAD.dna.p.effect.mannu.txt", "Processed_data/COAD.cnv.p.effect.mannu.amplification.txt", "Processed_data/COAD.cnv.p.effect.mannu.deletion.txt", "Processed_data/PANCAN.dna.p.effect.mannu.txt", "Processed_data/PANCAN.cnv.p.effect.mannu.amplification.txt", "Processed_data/PANCAN.cnv.p.effect.mannu.deletion.txt", "Data/Validation_CCLE/CCLE_mutations.csv", "Validation/CCLE.thresholded.cnv.tsv", "Data/Validation_CCLE/protein_quant_current_normalized.csv", "Data/NESTs/TheNEST.csv", "Validation/CCLE.COAD.cis.PolyRisk.dna.v2.txt", "Validation/CCLE.COAD.cis.PolyRisk.amp.v2.txt", "Validation/CCLE.COAD.cis.PolyRisk.del.v2.txt", "COAD", "Validation/CCLE.COAD.pancan.PolyRisk.dna.v2.txt", "Validation/CCLE.COAD.pancan.PolyRisk.amp.v2.txt", "Validation/CCLE.COAD.pancan.PolyRisk.del.v2.txt", "Data/Validation_CCLE/sample_info.csv"]

#wdna = pandas.read_csv(argv[1],sep="\t")
#wcnva = pandas.read_csv(argv[2],sep="\t")
#wcnva = pandas.read_csv(argv[2],sep="\t")
#wcnvd = pandas.read_csv(argv[3],sep="\t")
#pandna = pandas.read_csv(argv[4],sep="\t")
#pancnva = pandas.read_csv(argv[5],sep="\t")
#pancnvd = pandas.read_csv(argv[6],sep="\t")
#ccled = getMAFadj(argv[7],False)
#cclec = pandas.read_csv(argv[8],sep="\t") 
#cclep = pandas.read_csv(argv[9])
#nest = getNestDict(argv[10])
#odna = argv[11] 
#oamp = argv[12]
#odel = argv[13]
#cancer = argv[14]
#pdna = argv[15]
#pcnva = argv[16]
#pcnvd = argv[17]
#sampinfo = pandas.read_csv(argv[18])
#
##Validate finding in only these three cancer types. 
#my_dict = {i:sampinfo.primary_disease.to_list().count(i) for i in sampinfo.primary_disease.to_list()}
#
#
#todrop = ["Protein_Id","Description","Group_ID","Uniprot","Uniprot_Acc","TenPx01_Peptides","TenPx02_Peptides","TenPx03_Peptides","TenPx04_Peptides","TenPx06_Peptides","TenPx07_Peptides","TenPx08_Peptides","TenPx09_Peptides","TenPx10_Peptides","TenPx11_Peptides","TenPx12_Peptides","TenPx13_Peptides","TenPx16_Peptides","TenPx17_Peptides","TenPx19_Peptides","TenPx20_Peptides","TenPx21_Peptides","TenPx22_Peptides","TenPx27_Peptides","TenPx28_Peptides","TenPx29_Peptides","TenPx33_Peptides","TenPx34_Peptides","TenPx35_Peptides","TenPx36_Peptides","TenPx37_Peptides","TenPx38_Peptides","TenPx39_Peptides","TenPx40_Peptides","TenPx42_Peptides","TenPx23_Peptides","TenPx05_Peptides","TenPx30_Peptides","TenPx31_Peptides","TenPx32_Peptides","TenPx14_Peptides","TenPx15_Peptides","TenPx41_Peptides","TenPx26_Peptides","TenPx25_Peptides","TenPx18_Peptides","TenPx24_Peptides"]
#
#cclep2 = cclep.drop(todrop,axis=1)
#
##Notes: Colon = 83, Ovarian = 74, Breast = 83 
##Those are fine numbers to work with 
##In need to rework the header for the proteomics data to get to DEPMAP or CCLE files 
#yo = cclep.columns.to_list()
#yo2 = [i.rsplit("_",1)[0] for i in yo ]
#yo3 = [value for value in yo2 if value in sampinfo['CCLE_Name'].to_list()]
#dmid = dict(zip(sampinfo.CCLE_Name,sampinfo.DepMap_ID))
#depmapid = [dmid[i] for i in yo2 if i in dmid.keys()]
#todrop = ["Protein_Id","Gene_Symbol","Description","Group_ID","Uniprot","Uniprot_Acc","TenPx01_Peptides","TenPx02_Peptides","TenPx03_Peptides","TenPx04_Peptides","TenPx06_Peptides","TenPx07_Peptides","TenPx08_Peptides","TenPx09_Peptides","TenPx10_Peptides","TenPx11_Peptides","TenPx12_Peptides","TenPx13_Peptides","TenPx16_Peptides","TenPx17_Peptides","TenPx19_Peptides","TenPx20_Peptides","TenPx21_Peptides","TenPx22_Peptides","TenPx27_Peptides","TenPx28_Peptides","TenPx29_Peptides","TenPx33_Peptides","TenPx34_Peptides","TenPx35_Peptides","TenPx36_Peptides","TenPx37_Peptides","TenPx38_Peptides","TenPx39_Peptides","TenPx40_Peptides","TenPx42_Peptides","TenPx23_Peptides","TenPx05_Peptides","TenPx30_Peptides","TenPx31_Peptides","TenPx32_Peptides","TenPx14_Peptides","TenPx15_Peptides","TenPx41_Peptides","TenPx26_Peptides","TenPx25_Peptides","TenPx18_Peptides","TenPx24_Peptides"]
#
#cclep2 = cclep.drop(todrop)
#
#
#if cancer == "OV":
#    subset_can = sampinfo[sampinfo['primary_disease'] == 'Ovarian Cancer']
#elif cancer == "BRCA":
#    subset_can = sampinfo[sampinfo['primary_disease'] == 'Breast Cancer']
#elif cancer == "COAD":
#    subset_can = sampinfo[sampinfo['primary_disease'] == 'Colon/Colorectal Cancer']
#else:
#    subset_can = sampinfo 
#
