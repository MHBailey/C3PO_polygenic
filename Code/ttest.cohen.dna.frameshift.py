def readlist(fname):
    mylist = []
    with open(fname,"r") as f:
        for line in f:
            l = line[:-1]
            mylist.append(l)
    f.close()
    return mylist

def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (numpy.mean(x) - numpy.mean(y)) / numpy.sqrt(((nx-1)*numpy.std(x, ddof=1) ** 2 + (ny-1)*numpy.std(y, ddof=1) ** 2) / dof)


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
        (maf['Variant_Classification'] == "Nonsense_Mutation")]#all of these conditions 
        maf_sub = smaf[["Hugo_Symbol","Chromosome","Start_Position","Reference_Allele","Tumor_Sample_Barcode","Annotation_Transcript","Protein_Change","COHORT"]]
    return(maf_sub)

 
def ttest_cohen(muts,prots,nests,can,ofname):
    o = open(ofname,"w")
    if can != "PANCAN":
        cmuts = muts[muts["COHORT"] == can]
    else:
        cmuts = muts

    allsamps = cmuts["Proteome_Sample_ID"].unique().tolist() #can specific

    #NOTE keep track of these missing samples and (driver specific genes)  
    header = ["Protein","NESTv1","COHORT","SUBSTRATE","TSTAT","TPVALUE","COHEN_D","nMUT","nWT"]
    o.write("\t".join(header))
    o.write("\n")

    for i in nests:
        sgenes = nests[i] #This gets all genes in a nest 
        muts2 = cmuts[cmuts['Hugo_Symbol'].isin(sgenes)]
        mutsamps = muts2["Proteome_Sample_ID"].unique().tolist() #Samples with mut in any gene in nest
        #Subset PROTS
        prots2 = prots[prots.index.isin(allsamps)]

        mutated = prots2[prots2.index.isin(mutsamps)]
        wt = prots2[~prots2.index.isin(mutsamps)] 
        columns = prots2.columns.tolist()
        
        if len(mutated.index) > 5 and len(wt.index) > 5:
            for g in columns:
                a = mutated[g].dropna()
                b = wt[g].dropna()
                #Make sure there isn't missing data
                t_value,p_value=stats.mannwhitneyu(a,b,nan_policy='omit')
                cohensd = cohen_d(a,b)
                out = [str(g),str(i),str(can),"DNA",str(round(t_value,4)),str(p_value),str(round(cohensd,4)),str(len(a)),str(len(b))]
                o.write("\t".join(out))
                o.write("\n")
    o.close()

if __name__ == "__main__":
    import sys
    import numpy
    import pandas
    from scipy import stats
    nest = getNestDict(sys.argv[1]) 
    adjmaf = getMAFadj(sys.argv[2],False) #labeled with 299 Drivers or all genes NOTE turn this into a param
    prot = pandas.read_csv(sys.argv[3],sep="\t")
    meta = pandas.read_csv(sys.argv[4],sep="\t")
    #ADD PROT IDS to MAF
    mm = adjmaf.merge(meta,left_on='Tumor_Sample_Barcode',right_on='WXS',how='left')
    #REMOVE Normals from PROT 
    prott = prot.set_index('external_gene_name').T.rename_axis('Sample') #Transpose the table
    np = prott[prott.index.isin(meta['Proteome_Sample_ID'].tolist())]
    can = sys.argv[5]
    ofname = sys.argv[6]
    #TTEST and COHENS
    sys.setrecursionlimit(2500)
    ttest_cohen(mm,np,nest,can,ofname)
     
    
    
 
