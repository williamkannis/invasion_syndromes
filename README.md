# Source code for *Syndromes of multifaceted beta diversity change in invaded metacommunities*

# Anonymous peer review version

Code and data for reproducing the analyses presented in:

> Author et al. (YEAR). Syndromes of multidimensional beta diversity change in
invaded metacommunities.[Journal, DOI]

This repository contains the R code and data products required to
reproduce the analyses, figures, and tables presented in the manuscript.

The general null-model workflow used to generate the beta-diversity
metrics is available in the companion repository:

[Beta diversity null-model workflow](https://anonymous.4open.science/r/beta_null_hpc-6490/README.md)

For questions about this data set or analysis, please contact:

```bash
Name: William K. Annis

Email: wannis@fsu.edu, williamkannis@gmail.com

OrcID: 0009-0003-3541-8503
```

If you use this code or data, please cite:

>BLANK. Syndromes of multidimensional beta diversity change in 
invaded metacommunities. in review

>BLANK. Data and code for Syndromes of multidimensional beta diversity change in 
invaded metacommunities.[citation / DOI]

and if you adapt the null model workflow:

>BLANK. Beta diversity change null modelling workflow for Slurm-based high 
>performance computer (HPC) clusters. [Workflow citation / DOI]


# Overview
This repository is separated into two separate analysis routes:

1. Calculation of null model standardized effect sizes

2. Conducting of manuscript analyses

Replicating the analyses from route 1, requires the use of a high performance 
computer cluster (HPC) to replicate the effect size calculations. We provide
data, R scripts, and shell scripts to accomplish this task. For those wanting
to replicate the analyses without high performance computing, we also
provide intermediate effect size diversity data require to replicate analyses
in route 2.

<ins>NOTE:</ins> We are unable to directly provide the raw community or trait 
data used in our analyses. As such, this repository does not permit 
reconstruction of the observed community data from which the diversity metrics 
were originally calculated. Instead, we provide the observed diversity values 
and complete null-model iterations generated from those data. These allow 
independent reproduction of the effect-size calculations and all downstream 
manuscript analyses. Additionally we provide code for community and trait data
preparation [here](#diversity-input-data), and the beta diversity null modelling 
workflow used in this manuscript can be found 
[here](https://anonymous.4open.science/r/beta_null_hpc-6490/README.md).

# Data

The data below can be found at the manuscript's 
[Zenodo repository](https://doi.org/10.5281/zenodo.16949657):

### Observed diversity data
Observed alpha and beta diversity data can be downloaded from ```obs_out.zip``` 
at the [Zenodo repository](https://doi.org/10.5281/zenodo.16949657) and unzipped into ```HPC_data/``` directory. 
These data are used to estimate null model empirical effect sizes, used in 
analyses, and used for summary statistics. 

**Required** for both route 1 and route 2.

### Null diversity iterations
The raw null model iteration data for beta and alpha diversity can be downloaded 
from ```null_out.zip``` at the [Zenodo repository](https://doi.org/10.5281/zenodo.16949657) and unzipped into 
the ```HPC_data/``` directory. These data are the result of the randomization of 
traits and phylogenies (beta) or randomization of communities (alpha). These 
data will be compiled and used to create null distributions used to create null 
model standardized diversity values.

**Required** only for route 1.

### Summarized null model outputs
The summarized null model analyses results for all diversity metrics can be 
downloaded from ```ses_out.zip``` at the [Zenodo repository](https://doi.org/10.5281/zenodo.16949657) 
and unzipped into the ```HPC_data/``` directory. These data contain effect sizes,
can be used to evaluate the properties of the null distributions to choose 
between standardized effect  sizes and empirical effect sizes, and can be 
merged together for plotting and to create the final data.frames for the use in 
analyses.

**Required** only for route 2.

### Delta lcbd effects sizes
The combined effect size data for delta LCBD can be downloaded from ```delta_lcbd.rds``` 
at the [Zenodo repository](https://doi.org/10.5281/zenodo.16949657) and unzipped 
into the ```Analysis_data/``` directory. These data can be used to replicate
clustering analysis in ```02_cluster_analysis.R```.

Intermediate data, **NOT required** for replicating data using the full
workflows of route 1 or route 2.

### Community invadedness
We summarized the native and nonnative species richness across communities for each
syndrome. We used the raw community data to estimate the total species 
richness, native species richness, and the richness of nonnative species. 
These data can be downloaded from ```com_invaded.rds``` at the 
[Zenodo repository](https://doi.org/10.5281/zenodo.16949657) and unzipped 
into ```Analysis_data/```.

**Required** for both route 1 and route 2 to recreate figure 6.


## File directories
Download the entire repository. Then download required data from 
[Zenodo Repository](https://doi.org/10.5281/zenodo.16949657), and data sources listed 
[here](#publically-available-data), and unzip into the following 
file structure:

#UPDATE
```bash
├── Analysis_data
│   │── delta_lcbd*
│   └── com_invaded.rds*
│
├── Figures
├── HPC
│   │── null_out*
│   │── obs_out*
│   │── ses_inputs
│   └── ses_outputs*
│
└── Scripts


(*) directories or files downloaded from data sources other than current repository
```

<ins>NOTE:</ins> When downloading data from the [Zenodo repository](https://doi.org/10.5281/zenodo.16949657), the .zip 
files will create duplicate folders. It is recommended to download all data 
before creating new file structures.

# Workflows

## Required Software

**R version**: 4.5.0
#UPDATE

R packages

* ```'DescTools'``` version: 0.99.60
* ```'dplyr'``` version: 1.1.4
* ```'ggplot2'``` version: 4.0.2
* ```'matrixStats'``` version: 1.5.0
* ```'parallel'``` version: 4.5.0
* ```'purrr'``` version: 1.1.0
* ```'sf'``` version: 1.0.20
* ```'StreamCatTools'``` version: 0.10.0
* ```'tibble'``` version: 3.2.1
* ```'tidyr'``` version: 1.3.1
* ```'vegan'``` version: 2.6.10
* ```'VGAM'```  version: 1.1.14


## Route 1: Calculate effect sizes using HPC

R and shell scripts can be found in ```HPC``` directory and are listed in order of 
workflow. Users starting here do not need to download ```ses_out.zip```.

<ins>NOTE:</ins> Completing this part of workflow requires HPC cluster
with Slurm interface and downloading large files (```null_out.zip```). If this 
is not possible, skip to [Route 2](#route-2-conduct-analayses-using-local-machine).

### 1. Upload entire ```HPC``` directory to high performance cluster storage.
For all steps utilizing HPC clusters, users will run shell scripts that will
run the respective R script using specified HPC resources.

### 2. Prepare effect size input data
**Scripts:** 

* ```01_tax_beta_null_model_prep.sh``` - ```01_tax_beta_null_model_prep.R```
* ```02_beta_null_model_prep.sh``` - ```02_beta_null_model_prep.R```

<ins>Purpose:</ins> Calculates difference in LCBD between contemporary and 
native pools (delta) for the observed values, and for each null iteration. 
Consolidates outputs into single files, with separate delta, native, and 
contemporary species pool values.

<ins>Outputs:</ins> Intermediate files used for effect size calculations


### 3. Calculate effect sizes
**Scripts:** ```03_batch_ses.sh``` - ```03_batch_ses.R```

<ins>Purpose:</ins> Estimates standardize effect sizes (SES) of each single
metric using the custom function found in ```00_null_model_effect_size_function.R```.

This function summarizes null distributions and calculates 
standardize effect sizes in the traditional z score method (SES), empirical 
p-values, and p-value based effect sizes (ES), and reports optional diagnostic 
metrics used to select between the two effect size methods. Asymmetrical null 
distributions should be assessed using empirical p-value based effect sizes 
rather than z-score based SES. See 
[Botta-Dukát (2018)](https://doi.org/10.1556/168.2018.19.1.8) for more 
information on selecting SES or p-value based ES.

<ins>Outputs:</ins> 

* [Summarized null model outputs](#summarized-null-model-outputs)


### 3. Download the entire ```HPC``` directory to local machine.



## Route 2: Conduct analayses using local machine

R scripts can be found in ```Scripts``` directory and are listed in order of 
workflow. Users starting here do not need to download ```null_out.zip``` from
Zenodo repository


### 1. Summarize null model results
**Script:** ```01_ses_comp.R```

<ins>Purpose:</ins> Compiles and formats the resulting SES, ES, and diagnostic 
stats across files. Also visualizes normality diagnostics to allow users to 
decide between SES or ES values for further analyses, and creates plots
summarizing raw and effect size beta diversity values.

<ins>Outputs:</ins> 

* [```delta_lcbd.rds```](#delta-lcbd-effects-sizes), 
* Summary statistics
* Plots for figure 4


### 2. Perform clustering analysis to identify syndromes
**Script:** ```02_cluster_analysis.R```

<ins>Purpose:</ins> Performs model based clustering to identify invasion 
syndromes based on shared changes in multidimensional LCBD effect sizes. Filters
sites so only those with at least one dimension with a significant change are
included in cluster analyses. Creates summary tables and plots.

<ins>Outputs:</ins> 

* Clustering summary statistics
* Plots for figure 5
* Summary statistics for figure 6
* Plots for figure 6
* Plots for appendix 6
* Tables for appendix 7 - *Requires unprovided raw community data*


# Figures
Contains raw R plots and .csv tables used to create the figures and tables in 
the main manuscript and appendices. Most figures were edited in Adobe Illustrator
for aesthetic purposes, and we included Illustrator files as well.

# Diversity Input Data
We are not able to publicly provide all the data necessary to replicate the 
multidimensional diversity data estimated using the 
[HPC null model workflow](https://anonymous.4open.science/r/beta_null_hpc-6490/README.md). 
We do, however, provide the code used to prepare said data and can offer
data upon completed data requests. See below for more information.

The fish occurrence data were obtained through data sharing agreements 
with United States governmental agencies. While these raw data are not directly 
available from the authors for redistribution due to data sharing agreements, 
they can be accessed through formal requests to the agencies listed in Appendix 1. 
Individuals with completed data requests may contact the corresponding 
author for harmonized versions of the data. 

Trait data were obtained from public and private data sets and harmonized data 
cannot be shared without permission. Trait data can be found publicly available at
[Frimpong & Angermeier (2009)](https://www.sciencebase.gov/catalog/item/5a7c6e8ce4b00f54eb2318c0),
and upon request from the authors of 
[Giam & Olden (2016)](https://doi.org/10.1111/geb.12475).

Despite not being able to share all data, we provide R scripts used 
to harmonize the community ```community_data_prep.R``` and 
trait data ```trait_data_prep.R```, as well as the list of species of the 
entire data set ```full_species_list.csv```, and species list for the 
analysis ```filtered_species_list.csv```.
These scripts and data can be found in the ```Diversity Input Data``` directory.

The phylogenetic tree used for multidimensional diversity metrics, and 
information on its methodology can be found [here](https://anonymous.4open.science/r/fishscales_super_tree-01FC/README.md).
