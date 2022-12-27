def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (numpy.mean(x) - numpy.mean(y)) / numpy.sqrt(((nx-1)*numpy.std(x, ddof=1) ** 2 + (ny-1)*numpy.std(y, ddof=1) ** 2) / dof)



def delta_pairs(prots,meta,deltaout,can):
    o = open(deltaout,"w")

    if can != "PANCAN":
        cmeta = meta[meta["cohort"] == can]
    else:
        cmeta = meta

    tumor = cmeta['Proteome_Sample_ID'].unique().tolist()
    normal = cmeta['Proteome_Normal_Sample_ID'].unique().tolist()
    
    tumordf = prots[prots.index.isin(tumor)]
    normaldf = prots[prots.index.isin(normal)]

    if len(normaldf) > 5: #Some cancer types don't have normal
        PAIRS = cmeta[['Proteome_Sample_ID','Proteome_Normal_Sample_ID']].dropna() 
        PAIRS.reset_index(inplace=True, drop=True)

        header = ["Protein","TumorID","NormalID","DeltaTN"]
        o.write("\t".join(header))
        o.write("\n")

        columns = prots.columns.tolist()
        for i in range(len(PAIRS)):        
            for g in columns:
                tid = PAIRS.at[i,"Proteome_Sample_ID"]
                nid = PAIRS.at[i,"Proteome_Normal_Sample_ID"]
                tprot = tumordf.at[tid,g]
                nprot = normaldf.at[nid,g]
                deltap = tprot - nprot 
                out = [str(g),str(tid),str(nid),str(deltap)]
                o.write("\t".join(out))
                o.write("\n")
    o.close()



def tumor_normal(prots,meta,rawout,can):
    o = open(rawout,"w")

    if can != "PANCAN":
        cmeta = meta[meta["cohort"] == can]
    else:
        cmeta = meta

    tumor = cmeta['Proteome_Sample_ID'].unique().tolist()
    normal = cmeta['Proteome_Normal_Sample_ID'].unique().tolist()
    
    tumordf = prots[prots.index.isin(tumor)]
    normaldf = prots[prots.index.isin(normal)]


    header = ["Protein","NESTv1","COHORT","SUBSTRATE","TSTAT","TPVALUE","COHEN_D","nTumor","nNormal","meanTumor","meanNormal"]
    o.write("\t".join(header))
    o.write("\n")


    columns = prots.columns.tolist()
    for g in columns:
        a = tumordf[g].dropna()
        b = normaldf[g].dropna()
        if len(a) > 5 and len(b) > 5:
            t_value,p_value=stats.mannwhitneyu(a,b,nan_policy='omit')
            cohensd = cohen_d(a,b)
            meana = sum(a)/len(a)
            meanb = sum(b)/len(b)
            out = [str(g),"NA",can,"TN",str(round(t_value,4)),str(p_value),str(round(cohensd,4)),str(len(a)),str(len(b)),str(meana),str(meanb)]
            o.write("\t".join(out))
            o.write("\n")

    o.close()    




if __name__ == "__main__":
    import sys 
    import numpy
    import pandas 
    from scipy import stats 
    
    prot=pandas.read_csv(sys.argv[1],sep="\t")
    prott = prot.set_index('external_gene_name').T.rename_axis('Sample')
    meta = pandas.read_csv(sys.argv[2],sep="\t")
    can = sys.argv[3]
    rawout = sys.argv[4]
    deltaout = sys.argv[5]

    tumor_normal(prott,meta,rawout,can)
    delta_pairs(prott,meta,deltaout,can)
