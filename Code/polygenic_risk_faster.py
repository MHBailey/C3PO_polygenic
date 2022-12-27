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
        maf_sub = smaf[["Hugo_Symbol","Chromosome","Start_Position","Reference_Allele","Tumor_Sample_Barcode","Annotation_Transcript","Protein_Change","COHORT"]]
    return(maf_sub)

def makeDELMatrix(cnv,prots,nests,can):
    if can != "PANCAN":
        ccnv = cnv[cnv["cohort"] == can]
    else:
        ccnv = cnv

    allsamps = ccnv["Proteome_Sample_ID"].unique().tolist()

    SAMPxNEST = pandas.DataFrame()
    for i in nests:
        sgenes = nests[i] #This gets all genes in a nest 
        cnv2 = ccnv[ccnv['Gene Symbol'].isin(sgenes)]
        nAMP = len(cnv2[cnv2['value'] > 0])
        amp = cnv2[cnv2['value'] > 0]
        uamp = amp['variable'].unique().tolist() #Samples with amp
        nDEL = len(cnv2[cnv2['value'] < 0])
        dell = cnv2[cnv2['value'] < 0]
        udel = dell['variable'].unique().tolist()
        nWT = len(cnv2[cnv2['value'] == 0])
        cnv2test = cnv2[cnv2['variable'].isin(udel)] #BUILD AND LOOP TO DO THIS FOR DELETIONS AND AMP Seperately
        cnvsamps = cnv2test["Proteome_Sample_ID"].unique().tolist()
        
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
        
        df_wt = pandas.DataFrame(wt_tuples, columns = ['Proteome_Sample_ID', 'NEST', "Mut_WT"]) 
        df_mut = pandas.DataFrame(mut_tuples, columns = ['Proteome_Sample_ID', 'NEST', "Mut_WT"]) 
        mydf = pandas.concat([df_wt,df_mut]) 
        yuck = mydf.pivot_table(index='Proteome_Sample_ID',columns='NEST',values='Mut_WT')
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

    allsamps = ccnv["Proteome_Sample_ID"].unique().tolist()

    SAMPxNEST = pandas.DataFrame()
    for i in nests:
        sgenes = nests[i] #This gets all genes in a nest 
        cnv2 = ccnv[ccnv['Gene Symbol'].isin(sgenes)]
        nAMP = len(cnv2[cnv2['value'] > 0])
        amp = cnv2[cnv2['value'] > 0]
        uamp = amp['variable'].unique().tolist() #Samples with amp
        nDEL = len(cnv2[cnv2['value'] < 0])
        dell = cnv2[cnv2['value'] < 0]
        udel = dell['variable'].unique().tolist()
        nWT = len(cnv2[cnv2['value'] == 0])
        cnv2test = cnv2[cnv2['variable'].isin(uamp)] #BUILD AND LOOP TO DO THIS FOR DELETIONS AND AMP Seperately
        cnvsamps = cnv2test["Proteome_Sample_ID"].unique().tolist()

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

        df_wt = pandas.DataFrame(wt_tuples, columns = ['Proteome_Sample_ID', 'NEST', "Mut_WT"])
        df_mut = pandas.DataFrame(mut_tuples, columns = ['Proteome_Sample_ID', 'NEST', "Mut_WT"])
        mydf = pandas.concat([df_wt,df_mut])
        yuck = mydf.pivot_table(index='Proteome_Sample_ID',columns='NEST',values='Mut_WT')
        if SAMPxNEST.empty:
             SAMPxNEST = yuck
        else:
             SAMPxNEST = pandas.merge(SAMPxNEST,yuck,left_index=True, right_index=True, how='outer').fillna(0)
    return(SAMPxNEST)

def makeDNAMatrix(muts,prots,nests,can):
    if can != "PANCAN":
        cmuts = muts[muts["COHORT"] == can]
    else:
        cmuts = muts

    allsamps = cmuts["Proteome_Sample_ID"].unique().tolist() #can specific

    SAMPxNEST = pandas.DataFrame()
    for i in nests:
        sgenes = nests[i] #This gets all genes in a nest 
        muts2 = cmuts[cmuts['Hugo_Symbol'].isin(sgenes)]
        mutsamps = muts2["Proteome_Sample_ID"].unique().tolist() #Samples with mut in any gene in nest
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

        df_wt = pandas.DataFrame(wt_tuples, columns = ['Proteome_Sample_ID', 'NEST', "Mut_WT"])
        df_mut = pandas.DataFrame(mut_tuples, columns = ['Proteome_Sample_ID', 'NEST', "Mut_WT"])
        mydf = pandas.concat([df_wt,df_mut])
        yuck = mydf.pivot_table(index='Proteome_Sample_ID',columns='NEST',values='Mut_WT')
        if SAMPxNEST.empty:
             SAMPxNEST = yuck
        else:
             SAMPxNEST = pandas.merge(SAMPxNEST,yuck,left_index=True, right_index=True,how='outer').fillna(0)

    return(SAMPxNEST)


#def poly_genic_risk(wdna,muts,prots,nests,can,subs):



if __name__ == "__main__":
    import sys 
    import pandas
    import numpy
    from scipy import stats
    
    wdna = pandas.read_csv(sys.argv[1],sep="\t")
    wcnva = pandas.read_csv(sys.argv[2],sep="\t")
    wcnvd = pandas.read_csv(sys.argv[3],sep="\t")
    dna = getMAFadj(sys.argv[4],False)
    cnv = pandas.read_csv(sys.argv[5],sep="\t")
    nest = getNestDict(sys.argv[6])
    meta = pandas.read_csv(sys.argv[7],sep="\t")
    prot = pandas.read_csv(sys.argv[8],sep="\t")
    cancer = sys.argv[9]
    odna = sys.argv[10]
    oamp = sys.argv[11]
    odel = sys.argv[12]
    #
    dnam = dna.merge(meta,left_on='Tumor_Sample_Barcode',right_on='WXS',how='left')
    prott = prot.set_index('external_gene_name').T.rename_axis('Sample') # #Transpose the table

    #substrates for the data structure and so that I can add them later.
    substrates = ['dna','cnva','cnvd'] #Eventually add methylation H2a H2b H3 and H4. 


    mcnv = cnv.melt(id_vars = ["Gene Symbol","Gene ID","Cytoband"])

    #ADD PROT IDS to MAF
    mm = mcnv.merge(meta,left_on='variable',right_on='CASE_ID',how='left')

    #Make DNA Matrix
    DNA_MAT = makeDNAMatrix(dnam,prott,nest,cancer)
    print(DNA_MAT)

    #MAKE CNVa MATRIX
    AMP_MAT = makeAMPMatrix(mm,prott,nest,cancer)
    print(AMP_MAT)

    #DEL MAT
    DEL_MAT = makeDELMatrix(mm,prott,nest,cancer)
    print(DEL_MAT)

    #SET A PVALUE Threshold
    thresh = 0.1
    
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

    #TAKE A LOOK AT results
    print(dna_result)
    print(amp_result)
    print(del_result)

    
    dna_result.to_csv(odna, sep="\t")
    amp_result.to_csv(oamp, sep="\t")
    del_result.to_csv(odel, sep="\t")

    #Now i just need to is multiply the weights by the mutations, Sum the nests, and correlate with the proteins
    #Multiply dataframes
    
     
