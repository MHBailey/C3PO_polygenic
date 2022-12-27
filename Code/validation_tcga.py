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
    maf = pandas.read_csv(fname,sep="\t",low_memory=False)
    if(isDriver):
        maf_sub = maf
    else:
        smaf = maf[(maf['Variant_Classification'] == "Frame_Shift_Del") | \
        (maf['Variant_Classification'] == "Frame_Shift_Ins") | \
        (maf['Variant_Classification'] == "Missense_Mutation") | \
        (maf['Variant_Classification'] == "Nonsense_Mutation")]#all of these conditions 
        maf_sub = smaf[["Hugo_Symbol","Chromosome","Start_Position","Reference_Allele","Tumor_Sample_Barcode","Transcript_ID","HGVSp_Short","CODE"]]
    return(maf_sub)

def makeDELMatrix(cnv,prots,nests,can):
    if can != "PANCAN":
        ccnv = cnv[cnv["cohort"] == can]
    else:
        ccnv = cnv
    allsamps = ccnv["variable"].unique().tolist()
    SAMPxNEST = pandas.DataFrame()
    for i in nests:
        sgenes = nests[i] #This gets all genes in a nest 
        cnv2 = ccnv[ccnv['attrib_name'].isin(sgenes)]
        nAMP = len(cnv2[cnv2['value'] > 0])
        amp = cnv2[cnv2['value'] > 0]
        uamp = amp['variable'].unique().tolist() #Samples with amp
        nDEL = len(cnv2[cnv2['value'] < 0])
        dell = cnv2[cnv2['value'] < 0]
        udel = dell['variable'].unique().tolist()
        nWT = len(cnv2[cnv2['value'] == 0])
        cnv2test = cnv2[cnv2['variable'].isin(udel)] #BUILD AND LOOP TO DO THIS FOR DELETIONS AND AMP Seperately
        cnvsamps = cnv2test["variable"].unique().tolist()
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
        df_wt = pandas.DataFrame(wt_tuples, columns = ['variable', 'NEST', "Mut_WT"])
        df_mut = pandas.DataFrame(mut_tuples, columns = ['variable', 'NEST', "Mut_WT"])
        mydf = pandas.concat([df_wt,df_mut])
        yuck = mydf.pivot_table(index='variable',columns='NEST',values='Mut_WT')
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
    allsamps = ccnv["variable"].unique().tolist()
    SAMPxNEST = pandas.DataFrame()
    for i in nests:
        sgenes = nests[i] #This gets all genes in a nest 
        cnv2 = ccnv[ccnv['attrib_name'].isin(sgenes)]
        nAMP = len(cnv2[cnv2['value'] > 0])
        amp = cnv2[cnv2['value'] > 0]
        uamp = amp['variable'].unique().tolist() #Samples with amp
        nDEL = len(cnv2[cnv2['value'] < 0])
        dell = cnv2[cnv2['value'] < 0]
        udel = dell['variable'].unique().tolist()
        nWT = len(cnv2[cnv2['value'] == 0])
        cnv2test = cnv2[cnv2['variable'].isin(uamp)] #BUILD AND LOOP TO DO THIS FOR DELETIONS AND AMP Seperately
        cnvsamps = cnv2test["variable"].unique().tolist()
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
        df_wt = pandas.DataFrame(wt_tuples, columns = ['variable', 'NEST', "Mut_WT"])
        df_mut = pandas.DataFrame(mut_tuples, columns = ['variable', 'NEST', "Mut_WT"])
        mydf = pandas.concat([df_wt,df_mut])
        yuck = mydf.pivot_table(index='variable',columns='NEST',values='Mut_WT')
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

    allsamps = cmuts["char12"].unique().tolist() #can specific
    SAMPxNEST = pandas.DataFrame()
    for i in nests:
        sgenes = nests[i] #This gets all genes in a nest 
        muts2 = cmuts[cmuts['Hugo_Symbol'].isin(sgenes)]
        mutsamps = muts2["char12"].unique().tolist() #Samples with mut in any gene in nest
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
        df_wt = pandas.DataFrame(wt_tuples, columns = ['char12', 'NEST', "Mut_WT"])
        df_mut = pandas.DataFrame(mut_tuples, columns = ['char12', 'NEST', "Mut_WT"])
        mydf = pandas.concat([df_wt,df_mut])
        yuck = mydf.pivot_table(index='char12',columns='NEST',values='Mut_WT')
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
    tcgad = getMAFadj(sys.argv[7],False)
    tcgac = pandas.read_csv(sys.argv[8],sep="\t") 
    tcgap = pandas.read_csv(sys.argv[9],sep="\t")
    nest = getNestDict(sys.argv[10])
    odna = sys.argv[11] 
    oamp = sys.argv[12]
    odel = sys.argv[13]
    cancer = sys.argv[14]
    pdna = sys.argv[15]
    pcnva = sys.argv[16]
    pcnvd = sys.argv[17]

    #Clean up DNA
    tcgad['char12'] = tcgad.Tumor_Sample_Barcode.str.slice(0,12).str.replace("-",".").tolist()
    #CLEAN up Protein 
    prott = tcgap.set_index('attrib_name').T.rename_axis('Sample')
    #Clean up CNV
    mcnv = tcgac.melt(id_vars=['attrib_name'])
    mcnv['cohort'] = cancer
    #Make DNA Matrix
    DNA_MAT = makeDNAMatrix(tcgad,prott,nest,cancer)
    #print(DNA_MAT)

    #MAKE CNVa MATRIX
    AMP_MAT = makeAMPMatrix(mcnv,prott,nest,cancer)
    #print(AMP_MAT)

    #DEL MAT
    DEL_MAT = makeDELMatrix(mcnv,prott,nest,cancer)
    #print(DEL_MAT)

    #SET A PVALUE Threshold
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
    #argv = ["CODE",'Processed_data/BRCA.dna.p.effect.mannu.txt','Processed_data/BRCA.cnv.p.effect.mannu.amplification.txt','Processed_data/BRCA.cnv.p.effect.mannu.deletion.txt','Processed_data/PANCAN.dna.p.effect.mannu.txt','Processed_data/PANCAN.cnv.p.effect.mannu.amplification.txt','Processed_data/PANCAN.cnv.p.effect.mannu.deletion.txt','Data/Validation_TCGA/BRCA.maf','Data/Validation_TCGA/BRCA_GISTIC2_threshold.LinkedOmics.20220321.cgt','Data/Validation_TCGA/BRCA.LinkedOmics.20220321.tsv','Data/NESTs/TheNEST.csv','Validation/BRCA.PolyRisk.dna.v2.txt','Validation/BRCA.PolyRisk.amp.v2.txt','Validation/BRCA.PolyRisk.del.v2.txt','BRCA','Validation/BRCA.pancan.PolyRisk.dna.v2.txt','Validation/BRCA.pancan.PolyRisk.amp.v2.txt','Validation/BRCA.pancan.PolyRisk.del.v2.txt']

#wdna = pandas.read_csv(argv[1],sep="\t")
#wcnva = pandas.read_csv(argv[2],sep="\t")
#wcnva = pandas.read_csv(argv[2],sep="\t")
#wcnvd = pandas.read_csv(argv[3],sep="\t")
#pandna = pandas.read_csv(argv[4],sep="\t")
#pancnva = pandas.read_csv(argv[5],sep="\t")
#pancnvd = pandas.read_csv(argv[6],sep="\t")
#tcgad = getMAFadj(argv[7],False)
#tcgac = pandas.read_csv(argv[8],sep="\t") 
#tcgap = pandas.read_csv(argv[9],sep="\t")
#nest = getNestDict(argv[10])
#odna = argv[11] 
#oamp = argv[12]
#odel = argv[13]
#cancer = argv[14]
#pdna = argv[15]
#pcnva = argv[16]
#pcnvd = argv[17]
