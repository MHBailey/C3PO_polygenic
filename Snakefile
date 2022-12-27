configfile: 'config.yaml'

from pathlib import Path
import datetime as dt

DATE = dt.datetime.today().strftime("%Y_%m_%d")
DAY = dt.datetime.today().strftime("%d")

CANS=['BRCA','CCRCC','COAD','GBM','HNSCC','LSCC','LUAD','OV','PDAC','UCEC']
PANCAN=['BRCA','CCRCC','COAD','GBM','HNSCC','LSCC','LUAD','OV','PDAC','UCEC',"PANCAN"]
PAIRS=['CCRCC','COAD','HNSCC','LUAD','LSCC','OV','PDAC','UCEC','PANCAN']
TCGA=['COAD','BRCA','OV']
CCLE=['']


rule gen_ttests_dna:
    input:
        code='Code/ttest.cohen.dna.py',
        nest=config['NEST'],
        maf=config['FULLMAF'],
        proteomics=config['PROT'],
        meta=config['META']
    output:
        outdna='Processed_data/{can}.dna.p.effect.mannu.txt'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        python {input.code} {input.nest} {input.maf} {input.proteomics} {input.meta} {params.cancer} {output.outdna}
        '''

#rule plot_NESTs_dna:
#    input:
#        code='Code/plot_nests.dna.R'
#        nest=config['NEST'],
#        maf=config['FULLMAF'],
#        proteomics=config['PROT'],
#        meta=config['META']
        #NOTE: Be sure to finish this off, just an output file


rule gen_ttests_dna_miss:
    input:
        code='Code/ttest.cohen.dna.missense.py',
        nest=config['NEST'],
        maf=config['FULLMAF'],
        proteomics=config['PROT'],
        meta=config['META']
    output:
        outdna='Processed_data/{can}.dnam.p.effect.mannu.txt'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        python {input.code} {input.nest} {input.maf} {input.proteomics} {input.meta} {params.cancer} {output.outdna}
        '''


rule gen_ttests_dna_frame:
    input:
        code='Code/ttest.cohen.dna.frameshift.py',
        nest=config['NEST'],
        maf=config['FULLMAF'],
        proteomics=config['PROT'],
        meta=config['META']
    output:
        outdna='Processed_data/{can}.dnaf.p.effect.mannu.txt'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        python {input.code} {input.nest} {input.maf} {input.proteomics} {input.meta} {params.cancer} {output.outdna}
        '''



rule gen_ttest_cnv_AMP:
    input:
        code='Code/ttest.cohen.cnv.amp.py',
        nest=config['NEST'],
        cnv=config['FULLCNV'],
        proteomics=config['PROT'],
        meta=config['META'] 
    output:
        outdna='Processed_data/{can}.cnv.p.effect.mannu.amplification.txt'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        python {input.code} {input.nest} {input.cnv} {input.proteomics} {input.meta} {params.cancer} {output.outdna}
        '''

rule gen_ttest_cnv_DEL:
    input:
        code='Code/ttest.cohen.cnv.del.py',
        nest=config['NEST'],
        cnv=config['FULLCNV'],
        proteomics=config['PROT'],
        meta=config['META']
    output:
        outdna='Processed_data/{can}.cnv.p.effect.mannu.deletion.txt'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        python {input.code} {input.nest} {input.cnv} {input.proteomics} {input.meta} {params.cancer} {output.outdna}
        '''

rule protein_tn: 
    input:
        code='Code/tumor_normal_compare.py',
        proteomics=config['PROT'],
        meta=config['META'] 
    output:
        raw='Processed_data/{can}.tn.prots.txt',
        delta_pairs='Processed_data/{can}.pairs.prots.txt'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        python {input.code} {input.proteomics} {input.meta} {params.cancer} {output.raw} {output.delta_pairs}
        '''

rule plot_qc:
    input:
        qc='Code/qc_scores.R',
        tn='Processed_data/{can}.tn.prots.txt',
        pairs='Processed_data/{can}.pairs.prots.txt',
        dna='Processed_data/{can}.dna.p.effect.mannu.txt',
        proteomics=config['PROT'],
        g299='Data/GeneLists/299.genes.txt',
        gCPTAC='Data/GeneLists/CPTAC-Pancan-summary.txt',
        cnva='Processed_data/{can}.cnv.p.effect.mannu.amplification.txt',
        cnvd='Processed_data/{can}.cnv.p.effect.mannu.deletion.txt',
    output:
        bigfig='Figures/{can}.bigFig.V1.pdf',
        bigdata='Processed_data/{can}.bigFig.data.V1.txt',
        tndiff='Figures/{can}.tn.diff.V1.pdf',
        pairdiff='Figures/{can}.pairs.diff.V1.pdf'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.qc} {input.tn} {input.pairs} {input.dna} {input.proteomics} {input.g299} {input.gCPTAC} {input.cnva} {input.cnvd} {params.cancer} {output.bigfig} {output.bigdata} {output.tndiff} {output.pairdiff}
        '''

rule sample_weigts:
    input:
        code='Code/polygenic_risk_faster.py',
        wdna='Processed_data/{can}.dna.p.effect.mannu.txt',
        wcnva='Processed_data/{can}.cnv.p.effect.mannu.amplification.txt',
        wcnvd='Processed_data/{can}.cnv.p.effect.mannu.deletion.txt',
        maf=config['FULLMAF'],
        cnv=config['FULLCNV'],
        nest=config['NEST'],
        meta=config['META'],
        proteomics=config['PROT'],
    output:
        odna='Processed_data/{can}.PolyRisk.dna.v2.txt',
        oamp='Processed_data/{can}.PolyRisk.amp.v2.txt',
        odel='Processed_data/{can}.PolyRisk.del.v2.txt',
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        python {input.code} {input.wdna} {input.wcnva} {input.wcnvd} {input.maf} {input.cnv} {input.nest} {input.meta} {input.proteomics} {params.cancer} {output.odna} {output.oamp} {output.odel}
        '''

rule correlation_prs:
    input:
        code='Code/prs.plots.R',
        proteomics=config['PROT'],
        prsdna='Processed_data/{can}.PolyRisk.dna.v2.txt',
        prsamp='Processed_data/{can}.PolyRisk.amp.v2.txt',
        prsdel='Processed_data/{can}.PolyRisk.del.v2.txt',
    output:
        dna_corr='Processed_data/{can}.SampleCorrelations.dna.v1.txt',
        amp_corr='Processed_data/{can}.SampleCorrelations.amp.v1.txt',
        del_corr='Processed_data/{can}.SampleCorrelations.del.v1.txt',
        combined='Processed_data/{can}.SampleCorrelations.combined.v1.txt'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.proteomics} {input.prsdna} {input.prsamp} {input.prsdel} {params.cancer} {output.dna_corr} {output.amp_corr} {output.del_corr} {output.combined} 
        '''

rule combineC3PO:
    input:
        code='Code/combineC3PO.R',
        prsdna='Processed_data/{can}.PolyRisk.dna.v2.txt',
        prsamp='Processed_data/{can}.PolyRisk.amp.v2.txt',
        prsdel='Processed_data/{can}.PolyRisk.del.v2.txt',
    output:
        combined='Processed_data/{can}.C3PO.combined.v1.txt'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.prsdna} {input.prsamp} {input.prsdel} {params.cancer} {output.combined} 
        '''


rule tcga_validation:
    input: 
        code='Code/validation_tcga.py',
        wdna='Processed_data/{can}.dna.p.effect.mannu.txt',
        wcnva='Processed_data/{can}.cnv.p.effect.mannu.amplification.txt',
        wcnvd='Processed_data/{can}.cnv.p.effect.mannu.deletion.txt',
        pandna='Processed_data/PANCAN.dna.p.effect.mannu.txt',
        pancnva='Processed_data/PANCAN.cnv.p.effect.mannu.amplification.txt',
        pancnvd='Processed_data/PANCAN.cnv.p.effect.mannu.deletion.txt',
        tcgad='Data/Validation_TCGA/{can}.maf',
        tcgac='Data/Validation_TCGA/{can}_GISTIC2_threshold.LinkedOmics.20220321.cgt',
        tcgap='Data/Validation_TCGA/{can}.LinkedOmics.20220321.tsv',
        nest=config['NEST'],
    output:
        odna='Validation/{can}.cis.PolyRisk.dna.v2.txt',
        oamp='Validation/{can}.cis.PolyRisk.amp.v2.txt',
        odel='Validation/{can}.cis.PolyRisk.del.v2.txt',
        pdna='Validation/{can}.pancan.PolyRisk.dna.v2.txt',
        pamp='Validation/{can}.pancan.PolyRisk.amp.v2.txt',
        pdel='Validation/{can}.pancan.PolyRisk.del.v2.txt'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        python {input.code} {input.wdna} {input.wcnva} {input.wcnvd} {input.pandna} {input.pancnva} {input.pancnvd} {input.tcgad} {input.tcgac} {input.tcgap} {input.nest} {output.odna} {output.oamp} {output.odel} {params.cancer} {output.pdna} {output.pamp} {output.pdel}
        '''

rule tcga_prs_corrs:
    input:
        code='Code/validate_prs_tcga.R',
        tcgap='Data/Validation_TCGA/{can}.LinkedOmics.20220321.tsv',
        tdna='Validation/{can}.cis.PolyRisk.dna.v2.txt',
        tamp='Validation/{can}.cis.PolyRisk.amp.v2.txt',
        tdel='Validation/{can}.cis.PolyRisk.del.v2.txt',
        tpdna='Validation/{can}.pancan.PolyRisk.dna.v2.txt',
        tpamp='Validation/{can}.pancan.PolyRisk.amp.v2.txt',
        tpdel='Validation/{can}.pancan.PolyRisk.del.v2.txt',
    output:
        dna_corr='Validation/{can}.SampleCorrelations.dna.v1.txt',
        amp_corr='Validation/{can}.SampleCorrelations.amp.v1.txt',
        del_corr='Validation/{can}.SampleCorrelations.del.v1.txt',
        combined='Validation/{can}.SampleCorrelations.combined.v1.txt',
        pandna_corr='Validation/{can}.pancan.SampleCorrelations.dna.v1.txt',
        panamp_corr='Validation/{can}.pancan.SampleCorrelations.amp.v1.txt',
        pandel_corr='Validation/{can}.pancan.SampleCorrelations.del.v1.txt',
        pancombined='Validation/{can}.pancan.SampleCorrelations.combined.v1.txt',
        figdna='Figures/validation.{can}.dna.qqplot.pdf',
        figcnva='Figures/validation.{can}.cnva.qqplot.pdf',
        figcnvd='Figures/validation.{can}.cnvd.qqplot.pdf',
        figcomb='Figures/validation.{can}.combined.qqplot.pdf',
        figdnaP='Figures/validation.{can}.pdna.qqplot.pdf',
        figcnvaP='Figures/validation.{can}.pcnva.qqplot.pdf',
        figcnvdP='Figures/validation.{can}.pcnvd.qqplot.pdf',
        figcombP='Figures/validation.{can}.pcombined.qqplot.pdf',
    params:
        cancer=lambda wildcards: wildcards.can
    shell: 
        '''
        Rscript --quiet --vanilla {input.code} {input.tcgap} {input.tdna} {input.tamp} {input.tdel} {params.cancer} {input.tpdna} {input.tpamp} {input.tpdel} {output.dna_corr} {output.amp_corr} {output.del_corr} {output.combined} {output.pandna_corr} {output.panamp_corr} {output.pandel_corr} {output.pancombined} {output.figdna} {output.figcnva} {output.figcnvd} {output.figcomb} {output.figdnaP} {output.figcnvaP} {output.figcnvdP} {output.figcombP}
        '''

rule prep_CCLE_cnv:
    input:
        code='Code/clean_cnv_ccle.R',
        cnv='Data/Validation_CCLE/CCLE_gene_cn.csv',
    output: 
        ampthresh='Figures/CCLE.validation.cnv.MYC.pdf',
        delthresh='Figures/CCLE.validation.cnv.PTEN.pdf',
        cnvthresh='Validation/CCLE.thresholded.cnv.tsv'
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.cnv} {output.ampthresh} {output.delthresh} {output.cnvthresh}
        '''

rule CCLE_validation:
    input:
        code='Code/validation_ccle.py',
        wdna='Processed_data/{can}.dna.p.effect.mannu.txt',
        wcnva='Processed_data/{can}.cnv.p.effect.mannu.amplification.txt',
        wcnvd='Processed_data/{can}.cnv.p.effect.mannu.deletion.txt',
        pandna='Processed_data/PANCAN.dna.p.effect.mannu.txt',
        pancnva='Processed_data/PANCAN.cnv.p.effect.mannu.amplification.txt',
        pancnvd='Processed_data/PANCAN.cnv.p.effect.mannu.deletion.txt',
        ccled='Data/Validation_CCLE/CCLE_mutations.csv',
        cclec='Validation/CCLE.thresholded.cnv.tsv',
        cclep='Data/Validation_CCLE/protein_quant_current_normalized.csv',
        nest=config['NEST'],
        sampinfo='Data/Validation_CCLE/sample_info.csv'
    output:
        odna='Validation/CCLE.{can}.cis.PolyRisk.dna.v2.txt',
        oamp='Validation/CCLE.{can}.cis.PolyRisk.amp.v2.txt',
        odel='Validation/CCLE.{can}.cis.PolyRisk.del.v2.txt',
        pdna='Validation/CCLE.{can}.pancan.PolyRisk.dna.v2.txt',
        pamp='Validation/CCLE.{can}.pancan.PolyRisk.amp.v2.txt',
        pdel='Validation/CCLE.{can}.pancan.PolyRisk.del.v2.txt',
        newtable='Validation/CCLE.{can}.updated.protein.table.tsv'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        python {input.code} {input.wdna} {input.wcnva} {input.wcnvd} {input.pandna} {input.pancnva} {input.pancnvd} {input.ccled} {input.cclec} {input.cclep} {input.nest} {output.odna} {output.oamp} {output.odel} {params.cancer} {output.pdna} {output.pamp} {output.pdel} {input.sampinfo} {output.newtable}
        '''

rule ccle_prs_corrs:
    input:
        code='Code/validate_prs_ccle.R',
        tcgap='Validation/CCLE.{can}.updated.protein.table.tsv',
        tdna='Validation/CCLE.{can}.cis.PolyRisk.dna.v2.txt',
        tamp='Validation/CCLE.{can}.cis.PolyRisk.amp.v2.txt',
        tdel='Validation/CCLE.{can}.cis.PolyRisk.del.v2.txt',
        tpdna='Validation/CCLE.{can}.pancan.PolyRisk.dna.v2.txt',
        tpamp='Validation/CCLE.{can}.pancan.PolyRisk.amp.v2.txt',
        tpdel='Validation/CCLE.{can}.pancan.PolyRisk.del.v2.txt',
    output:
        dna_corr='Validation/CCLE.{can}.SampleCorrelations.dna.v1.txt',
        amp_corr='Validation/CCLE.{can}.SampleCorrelations.amp.v1.txt',
        del_corr='Validation/CCLE.{can}.SampleCorrelations.del.v1.txt',
        combined='Validation/CCLE.{can}.SampleCorrelations.combined.v1.txt',
        pandna_corr='Validation/CCLE.{can}.pancan.SampleCorrelations.dna.v1.txt',
        panamp_corr='Validation/CCLE.{can}.pancan.SampleCorrelations.amp.v1.txt',
        pandel_corr='Validation/CCLE.{can}.pancan.SampleCorrelations.del.v1.txt',
        pancombined='Validation/CCLE.{can}.pancan.SampleCorrelations.combined.v1.txt',
        figdna='Figures/CCLE.validation.{can}.dna.qqplot.pdf',
        figcnva='Figures/CCLE.validation.{can}.cnva.qqplot.pdf',
        figcnvd='Figures/CCLE.validation.{can}.cnvd.qqplot.pdf',
        figcomb='Figures/CCLE.validation.{can}.combined.qqplot.pdf',
        figdnaP='Figures/CCLE.validation.{can}.pdna.qqplot.pdf',
        figcnvaP='Figures/CCLE.validation.{can}.pcnva.qqplot.pdf',
        figcnvdP='Figures/CCLE.validation.{can}.pcnvd.qqplot.pdf',
        figcombP='Figures/CCLE.validation.{can}.pcombined.qqplot.pdf',
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.tcgap} {input.tdna} {input.tamp} {input.tdel} {params.cancer} {input.tpdna} {input.tpamp} {input.tpdel} {output.dna_corr} {output.amp_corr} {output.del_corr} {output.combined} {output.pandna_corr} {output.panamp_corr} {output.pandel_corr} {output.pancombined} {output.figdna} {output.figcnva} {output.figcnvd} {output.figcomb} {output.figdnaP} {output.figcnvaP} {output.figcnvdP} {output.figcombP}
        '''


rule cptac_hallmarks:
    input:
        code='Code/hallmakers.v1.R',
        prsdna='Processed_data/{can}.PolyRisk.dna.v2.txt',
        prsamp='Processed_data/{can}.PolyRisk.amp.v2.txt',
        prsdel='Processed_data/{can}.PolyRisk.del.v2.txt',
        halls='Data/Hallmarks/nanostring.gl.txt',
        meta=config['META'],
        hallyaml='./hallmarks.v2.yaml'
    output:
        figs='Figures/CPTAC.{can}.hallmark.boxplots.pdf',
        samphall='Processed_data/CPTAC.{can}.hallmark.scores.txt',
        circos='Figures/Circos.CPTAC.{can}.hallmarks.pdf'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.prsdna} {input.prsamp} {input.prsdel} {params.cancer} {input.halls} {output.figs} {output.samphall} {output.circos} {input.meta} {input.hallyaml}
        '''

rule cptac_all:
#This code may need to change based on the location of cptac_hallmarks output.
    input:
        code='Code/all_circos.R',
        hallyaml='./hallmarks.v2.yaml',
        meta=config['META'],
    output:
        bigfig='Figures/top3_all_circos.pdf'
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.hallyaml} {output.all_top_circos}
        '''

rule cptac_r2:
    input:
        code='Code/percent_explained.R',
        prsdna='Processed_data/{can}.PolyRisk.dna.v2.txt',
        prsamp='Processed_data/{can}.PolyRisk.amp.v2.txt',
        prsdel='Processed_data/{can}.PolyRisk.del.v2.txt',
        halls='Data/Hallmarks/nanostring.gl.txt',
        meta=config['META'],
    output:
        figs='Figures/CPTAC.{can}.hallmark.boxplots.pdf',
        samphall='Processed_data/CPTAC.{can}.hallmark.scores.txt',
        circos='Figures/Circos.CPTAC.{can}.hallmarks.pdf'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.prsdna} {input.prsamp} {input.prsdel} {params.cancer} {input.halls} {output.figs} {output.samphall} {output.circos} {input.meta} {input.hallyaml}
        '''
 
rule tcga_hallmarks:
    input:
        code='Code/hallmakers.tcga.R',
        tdna='Validation/{can}.cis.PolyRisk.dna.v2.txt',
        tamp='Validation/{can}.cis.PolyRisk.amp.v2.txt',
        tdel='Validation/{can}.cis.PolyRisk.del.v2.txt',
        halls='Data/Hallmarks/nanostring.gl.txt',
        tpdna='Validation/{can}.pancan.PolyRisk.dna.v2.txt',
        tpamp='Validation/{can}.pancan.PolyRisk.amp.v2.txt',
        tpdel='Validation/{can}.pancan.PolyRisk.del.v2.txt',
    output:
        figs='Figures/TCGA.{can}.hallmark.boxplots.pdf',
        figsp='Figures/TCGA.pancan.{can}.hallmark.boxplots.pdf'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.tdna} {input.tamp} {input.tdel} {params.cancer} {input.halls} {output.figs} {input.tpdna} {input.tpamp} {input.tpdel} {output.figsp}
        '''

rule ccle_hallmarks:
    input:
        code='Code/hallmakers.ccle.R',
        tdna='Validation/CCLE.{can}.cis.PolyRisk.dna.v2.txt',
        tamp='Validation/CCLE.{can}.cis.PolyRisk.amp.v2.txt',
        tdel='Validation/CCLE.{can}.cis.PolyRisk.del.v2.txt',
        halls='Data/Hallmarks/nanostring.gl.txt',
        tpdna='Validation/CCLE.{can}.pancan.PolyRisk.dna.v2.txt',
        tpamp='Validation/CCLE.{can}.pancan.PolyRisk.amp.v2.txt',
        tpdel='Validation/CCLE.{can}.pancan.PolyRisk.del.v2.txt',
    output:
        figs='Figures/CCLE.{can}.hallmark.boxplots.pdf',
        figsp='Figures/CCLE.pancan.{can}.hallmark.boxplots.pdf'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.tdna} {input.tamp} {input.tdel} {params.cancer} {input.halls} {output.figs} {input.tpdna} {input.tpamp} {input.tpdel} {output.figsp}
        '''


rule clinical:
    input: 
        code='Code/clinical_relevance.R',
        samphall='Processed_data/CPTAC.{can}.hallmark.scores.txt',
        clinical=config['CLUST'],
    output:
        story1='Processed_data/immune_correlations.{can}.txt',
        s1fig='Figures/C3POxImmune.{can}.pdf',
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.samphall} {input.clinical} {params.cancer} {output.story1} {output.s1fig} 
        '''

rule cptac_hallmarks_protein:
    input:
        code='Code/hallmakers.v2.protein.R',
        proteomics=config['PROT'],
        halls='Data/Hallmarks/nanostring.gl.txt',
        meta=config['META'],
    output:
        figs='Figures/ProteinOnly.CPTAC.{can}.hallmark.boxplots.pdf',
        samphall='Processed_data/ProteinOnly.CPTAC.{can}.hallmark.scores.txt',
        circos='Figures/ProteinOnly.Circos.CPTAC.{can}.hallmarks.pdf'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.proteomics} {params.cancer} {input.halls} {output.figs} {output.samphall} {output.circos} {input.meta}
        '''


rule hallmarks_plus:
    input:
        code='Code/hallmakers.v3.R',
        prsdna='Processed_data/{can}.PolyRisk.dna.v2.txt',
        prsamp='Processed_data/{can}.PolyRisk.amp.v2.txt',
        prsdel='Processed_data/{can}.PolyRisk.del.v2.txt',
        halls='Data/Hallmarks/nanostring.gl.txt',
        meta=config['META'],
        hallyaml='hallmarks.v2.yaml',
        proteomics=config['PROT']
    output:
        figs='Figures/ProteinOnly.plusProtein.CPTAC.{can}.hallmark.boxplots.pdf',
        samphall='Processed_data/ProteinOnly.CPTAC.plusProtein.{can}.hallmark.scores.txt',
        circos='Figures/Circos.plusProtein.CPTAC.{can}.hallmarks.pdf'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.prsdna} {input.prsamp} {input.prsdel} {params.cancer} {input.halls} {output.figs} {output.samphall} {output.circos} {input.meta} {input.hallyaml} {input.proteomics}
        '''
#NEST paper: https://www.science.org/doi/10.1126/science.abf3067
#Cuckoos and Cowbirds
rule extra:
    input:
        code="Code/trace.R",
        halls='Data/Hallmarks/nanostring.gl.txt',
        meta=config['META'],
        nest=config['NEST'],
        maf=config['FULLMAF'],
        cnv=config['FULLCNV'],
        prot=config['PROT'],
        prsdna='Processed_data/UCEC.PolyRisk.dna.v2.txt',
        prsamp='Processed_data/UCEC.PolyRisk.amp.v2.txt',
        prsdel='Processed_data/UCEC.PolyRisk.del.v2.txt',
        cnva='Processed_data/UCEC.cnv.p.effect.mannu.amplification.txt',
        cnvd='Processed_data/UCEC.cnv.p.effect.mannu.deletion.txt',
        dnaw='Processed_data/UCEC.dna.p.effect.mannu.txt',
    output:
        contfigs='Figures/ONEOFF.UCEC.Immune.pdf', 
        contribs='Processed_data/ONEOFF.UCEC.Immune.pdf', 
    params:
        sample="X01BR043",
        hallmark="Adaptive_immunity"
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.halls} {input.nest} {input.maf} {input.cnv} {input.prsdna} {input.prsamp} {input.prsdel}
        '''

rule percent_r2:
#This is going to read through all of the Correlations tables and build some overall plots
    input:
        code='Code/percent_r2.R',
        meta=config['META'],
    output:
        figR2='Figures/percent_captured.pdf'
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.meta} {output.figR2}  
        '''

rule sample_hallmarks:
    input:
        code='Code/hall_heatmap.v2.R',
        halls="Processed_data/CPTAC.{can}.hallmark.scores.txt",
    output:
        fig='Figures/CPTAC.{can}.hall_heatmap.v2.pdf'
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.halls} {params.cancer} {output.fig}  
        '''

rule ssGSEA:
#source activate heatmap
    input:
        code="ssGSEA_2_heatmap.R",
        dna="Data/ssGSEA/Pancan_driver_DNA_ssGSEA-combined.gct",
        rna="Data/ssGSEA/Pancan_driver_RNA_ssGSEA-combined.gct",
        protein="Data/ssGSEA/Pancan_driver_protein_ssGSEA-combined.gct",
        meta=config['META']
    output:
        dna_out="Processed_data/BRCA.ssGSEA.sample.matrix.RNA.txt",
        dna_fig="Figures/BRCA.ssGSEA.heatmap.RNA.pdf",
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.dna} {input.rna} {input.protein}{input.meta} {output.dna_out} {output.out_fig}
        ''' 

rule ssGSEA_scores:
    input:
        code="Code/score.ssGSEQ.R",
        dna_out="Processed_data/{can}.ssGSEA.sample.matrix.DNA.txt",
        rna_out="Processed_data/{can}.ssGSEA.sample.matrix.RNA.txt",
        protein="Processed_data/{can}.ssGSEA.sample.matrix.Protein.txt",
    output:
        rankS="Processed_data/ssGSEA.hallmark.wilcox.signed.SampleRank.{can}.txt",
        rankSfig="Figures/ssGSEA.pvalue.distrib.{can}.SampleRank.pdf",
        rankH="Processed_data/ssGSEA.hallmark.wilcox.signed.HallmarkRank.{can}.txt",
        rankHfig="Figures/ssGSEA.pvalue.distrib.{can}.HallmarkRank.pdf",
    params:
        cancer=lambda wildcards: wildcards.can    
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.dna_out} {input.rna_out} {input.protein} {params.cancer} {output.rankS} {output.rankSfig} {output.rankH} {output.rankHfig}
        '''

rule add_features_C3PO:
#source activate heatmap
    input: 
        code='Code/add_features_C3PO.R',
        c3po='Processed_data/{can}.C3PO.combined.v1.txt',
        hall='Processed_data/CPTAC.{can}.hallmark.scores.txt',
        meta=config['META'],
        smoking="Data/Meta_table/count.Neoantigen.Mutation.AllCancerType.v2.smoking.tsv",
    output: 
        fig="Figures/C3PO.{can}.Smoking.features.heatmap.pdf"
    params:
        cancer=lambda wildcards: wildcards.can
    shell:
        '''
        Rscript --quiet --vanilla {input.code} {input.c3po} {input.hall} {input.meta} {input.smoking} {params.cancer} {output.fig}
        '''

rule TMB:
    input:
        fig=expand("Figures/C3PO.{can}.Smoking.features.heatmap.pdf",can=CANS)


rule ssGSEA_wilcox:
    input:
        wilcox=expand("Figures/ssGSEA.pvalue.distrib.{can}.SampleRank.pdf",can=CANS)


rule ssGSEA_plots:
    input:
        dna_o="Processed_data/ssGSEA.dna.txt",

rule basic_heat_hall:
    input:
        fig=expand('Figures/CPTAC.{can}.hall_heatmap.v2.pdf',can=CANS)



rule hall_protein_plus:
    input:
        figs=expand('Figures/Circos.plusProtein.CPTAC.{can}.hallmarks.pdf',can=CANS)

rule hall_protein:
    input:
        figs=expand('Figures/ProteinOnly.CPTAC.{can}.hallmark.boxplots.pdf',can=CANS)

rule hall_cptac:
    input:
        figs=expand('Figures/CPTAC.{can}.hallmark.boxplots.pdf',can=CANS)


rule hall_tcga:
    input:
        figs=expand('Figures/TCGA.{can}.hallmark.boxplots.pdf',can=TCGA)

rule hall_ccle:
    input:
        figs=expand('Figures/CCLE.{can}.hallmark.boxplots.pdf',can=TCGA)

rule cclev:
    input:
        valid_ccle=expand('Validation/CCLE.{can}.SampleCorrelations.dna.v1.txt',can=TCGA)

rule tcgav:
    input:
        valid_tcga=expand('Validation/{can}.SampleCorrelations.dna.v1.txt',can=TCGA)

rule ttests_dna:
    input:
        dna=expand('Processed_data/{can}.dna.p.effect.mannu.txt',can=CANS)

rule ttests_cnv_amp:
    input:
        dna=expand('Processed_data/{can}.cnv.p.effect.mannu.amplification.txt',can=CANS)

rule ttests_cnv_del:
    input:
        dna=expand('Processed_data/{can}.cnv.p.effect.mannu.deletion.txt',can=CANS)

rule tumor_normal:
    input:
        tn=expand('Processed_data/{can}.tn.prots.txt',can=CANS)

rule qc:
    input:
        qc=expand('Figures/{can}.bigFig.V1.pdf',can=PAIRS)

rule next3:
    input:
        weights=expand('Processed_data/{can}.PolyRisk.dna.v2.txt',can=CANS)

rule corrs:
    input:
        corrs=expand('Processed_data/{can}.SampleCorrelations.dna.v1.txt',can=CANS)

rule combC3PO:
    input:
        combined=expand('Processed_data/{can}.C3PO.combined.v1.txt',can=CANS)

