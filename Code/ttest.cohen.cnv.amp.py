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

def ttest_cohen(cnv,prots,nests,can,ofname):
    #Think about CNV AMP, DEL directionality in NESTs>? 
    #NOTE: code this up for driver specific stuff too.
    #cnv = mm
    #prots = np
    #nests = nest

    o = open(ofname,"w")
    if can != "PANCAN":
        ccnv = cnv[cnv["cohort"] == can]
    else:
        ccnv = cnv

    allsamps = ccnv["Proteome_Sample_ID"].unique().tolist() #can specific list

    header = ["Protein","NESTv1","COHORT","SUBSTRATE","TSTAT","TPVALUE","COHEN_D","nMUT","nWT","nAMP","nDEL"]
    o.write("\t".join(header))
    o.write("\n")
    for i in nests:
        #i = "NEST:333"
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

        #Subset to can PROTS
        prots2 = prots[prots.index.isin(allsamps)]

        mutated = prots2[prots2.index.isin(cnvsamps)]
        wt = prots2[~prots2.index.isin(cnvsamps)] 
        columns = prots2.columns.tolist()
        if len(mutated.index) > 5 and len(wt.index) > 5:
            for g in columns:
                a = mutated[g].dropna()
                b = wt[g].dropna()
                #Make sure there isn't missing data
                if len(a) > 5 and len(b) > 5:
                    t_value,p_value=stats.mannwhitneyu(a,b,nan_policy='omit')
                    cohensd = cohen_d(a,b)
                    out = [str(g),str(i),can,"CNV",str(round(t_value,4)),str(p_value),str(round(cohensd,4)),str(len(a)),str(len(b)),str(nAMP),str(nDEL)]
                    o.write("\t".join(out))
                    o.write("\n")
    o.close()

if __name__ == "__main__":
    import sys
    import numpy
    import pandas
    from scipy import stats
    sys.setrecursionlimit(2500)

    nest = getNestDict(sys.argv[1]) 
    #nest = getNestDict(argv[1])
    cnv = pandas.read_csv(sys.argv[2],sep="\t")
    #cnv = pandas.read_csv(argv[2],sep="\t")
    prot = pandas.read_csv(sys.argv[3],sep="\t")
    #prot = pandas.read_csv(argv[3],sep="\t")

    meta = pandas.read_csv(sys.argv[4],sep="\t")
    #meta = pandas.read_csv(argv[4],sep="\t")
    mcnv = cnv.melt(id_vars = ["Gene Symbol","Gene ID","Cytoband"])

    #ADD PROT IDS to MAF
    mm = mcnv.merge(meta,left_on='variable',right_on='CASE_ID',how='left')
    #REMOVE Normals from PROT 
    prott = prot.set_index('external_gene_name').T.rename_axis('Sample') #Transpose the table
    np = prott[prott.index.isin(meta['Proteome_Sample_ID'].tolist())]

    ######## Bring in Cancers ####### 
    can = sys.argv[5]
    #can = "PANCAN"
    ofname = sys.argv[6]
    #ofname = "testcnv"

    #TTEST and COHENS
    ttest_cohen(mm,np,nest,can,ofname)
     
    
    
#argv = ["CODE","Data/NESTs/TheNEST.csv","Data/Somatic_cnv/Broad_pipeline_wxs_formatted/GISTIC_all_thresholded.by_genes.tsv","Data/Proteome/Broad_updated_tumor_NAT/Proteome_Broad_updated_tumor_NAT_imputed_gene_level.tsv","Data/Meta_table/CPTAC-Pancan-Data_metatable.txt","PANCAN","Processed_data/PANCAN.cnv.p.effect.mannu.txt"]
 
