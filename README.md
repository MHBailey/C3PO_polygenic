# C3PO_polygenic
Basic linear polygenic risk score to predict protein levels in CPTAC


## Background: Polygenic risk scores 
Simply put, this project uses very simple linear algebra to calculalte a Combined Polygenic Protein Prediction in Onocolgy (C3PO). I hope to develop this more over the next few years to incoporporate Mendeliean randomization, colocalization, and many more covariates. 

Like with all polygenic scores the foundation of the score centers on two principles. 1) Genetic variation and 2) the effect of that variation on a phenotype. If you do this over enough samples then you estimate the effect of any given change on a protein. For large germline studies this can be accomplished using Germline mutations. We looked to explore some biologically guided steps to apply this approach to somatic mutations. This is apparent in the code, but we leveraged [NeST_v1.0](https://idekerspaper.com) to collapse mutations into these networks. Then using a basic union, looked at the effect of DNA varition in these nests, specificially DNA mutations and copy number events. 

Training and test data was performed on [CPTCA](https://cptac-data-portal.georgetown.edu/cptac/) data, that is discussed in more detail here: [CPTAC explained](https://proteomics.cancer.gov/programs/cptac). 


## Installation tips: 
We use conda to library and package management. Basic scripts (written in R and Python) are organized and executed using snakemake [SNAKEMAKE](https://snakemake.readthedocs.io/en/stable/). I have found that the conda enviornmental replictation doesn't always get snakemake right. If the conda doesn't install snakemake properly, please use their code to get snakemake installed.

## Installation setting upt the environment:
STEP 1) Install conda according to [documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
  [Windows](https://conda.io/projects/conda/en/latest/user-guide/install/windows.html)
  [macOS](https://conda.io/projects/conda/en/latest/user-guide/install/macos.html)
  [Linux](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html)

STEP 2) Mimic the conda enviornment that accompanies this code. 
~~~
   conda env create -f C3PO.environment.yml 
~~~
NOTE: If a "Pip suprocess error" is shown. To my knowledge, this can be ignored.

STEP 3) Enter the conda environment 
~~~
   source activate C3PO
~~~  
or
~~~
   conda activate C3PO
~~~
Then add a coupld more libraries that didn't seem to come with the yaml. 
   conda install -c conda-forge pandas


STEP 4) Check snakemake functionality 
~~~
   which snakemake
~~~
NOTE: If this command returns an empty string or comes back empty, please install snakemake using the following commands (while in the conda environment).
~~~
   conda install -c bioconda snakemake
~~~  

STEP 5) This next step probably breaks many conda rules against builing nested environments, but it was "a" way to get ComplexHeatmap working with the more up-to-date R libraries. This will be called before generated the heatmap figures as separate snakemake runs. I will specify below when to enter this alternative env. While in the C3PO environment create a nested environment with the following code: 
~~~
   conda env create -f ComplexHeatmap.environment.yml
~~~  
NOTE: you may have to repeat step 4 within the ComplexHeatmap conda environment. 

STEP 6) Reconfigure the config.yaml file to match the data files. 
NOTE: These are not available at this time but the full paths are shown, and Data are not uploaded to GitHub. I will provide these links when these data become widely available. 

STEP  7) Rules that I can and should run to reproduce the figure


STEP 7.1) Generate the circos plots.
~~~
   #Dry run
   snakemake -np hall_protein_plus 

   #Actual run
   snakemake -p hall_protein_plus -c1
~~~ 

STEP 7.2) Generate the C3PO heatmaps 
~~~
   source activate ComplexHeatmap
   snakemake -p TMB -c1 
   conda deactivate 
~~~ 

STEP 7.3) Validation with TCGA samples (CCLE-TCGA samples)
~~~
   snakemake -p tcgav -c1 
~~~

STEP 7.4) Validation with CCLE samples
~~~
   snakemake -p cclev -c1 
~~~

STEP 7.5) Generate the entropy plots
~~~
   snakemake -p percent_r2 -c1 
~~~


Trouble shooting) 
1) While I was building this from scratch, I rand into the follwing error. https://github.com/snakemake/snakemake/issues/1899. You may need to downgrade to tabulate=0.8.10. Trying that now. To fix use the following: 
    conda install -c conda-forge mamba
    mamba install 'tabulate=0.8.10'

Please recognize that this is a preliminary attempt at applying polygenic risk scores to better understand the genomic impact of somatic mutations on protein abudances. We plan, over the course of the next few years, to optimize this algorithm and overcome the many simple assumptions we are making with these data. Additionally, we plan to add a methylation and germline component to our model in the near future as the data and strategies become available. 

Thank you for using the tool and we hope you find it useful in your own research. 


All my best,
MHBailey 




##RULES TO KEEP
gen_ttest_cnv_AMP
gen_ttest_cnv_DEL
gen_ttests_dna   
hallmarks_plus   
sample_weigts 

end rules: hall_protein_plus
