# SpeciationClock

This repository gathers the data and the scripts for data analyses corresponding to the article ***Speciation accelerates in extreme inbreeders*** by Gizaw, Gustafsson et al. It contains to main folder: one for the treatment of GBS data and one for statistical data analyses

## gbs
- **data**
	- ***contigs***: list of gbs contigs information for each species (including contig length)
	- ***vcf***: vcf files for each species
- **scripts**: to treat and analyse vcf files
- **output**: popgen statistics and SNP information for each species


## data_analyses
- **data**: initial datasets
	- ***experiment_data***: data from crosses experiment
	- ***sfs***: site frequency spectra from GBS data
	- other data files
- **scripts**
	- ***dadi***: python scripts for running ∂a∂i
	- R scripts for statistical analyses
- **output**: intermediate files
	- ***dadi_results***: results from ∂a∂i models
	- ***gbs_results***: results from gbs data analyses
	- ***synthetized data***: cleaned experimental data and combinaton with genetic data
- **results**: final results files
- **figures**: final figures