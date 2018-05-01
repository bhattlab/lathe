# Metagenomics Workflows
Workflows for metagenomic sequence data processing and analysis.  Further documentation found in each workflow folder.

## bin_label_and_evaluate

Snakemake workflow for aligning, binning, classifying and evaluating a
metagenomic assembly.

## assembly_comparison_circos
Snakemake workflow for visualizing assemblies of a particular genome across conditions and time points.  Calls out pre-identified sequences, highlights selected contigs.



# Installing workflows

First, install miniconda3. This is an environment management system that should keep everything organized.

Once installed, clone this github directory to some location where it can be stored permanently.

    git clone https://github.com/elimoss/metagenomics_workflows.git

Then enter the newly downlaoded mustache directory

    cd metagenomics_workflows

And create the new conda environment with

	conda env create -f envs/environment.yml

Future updates can be done with 
	conda env update -f envs/environment.yaml
	
Now activate the environment with

    conda activate mgwf

(This is a step that must be repeated whenever using any metagenomics workflow.)


This repository is forked from elimoss/metagenomics_workflows.  To sync updates, 
	git fetch upstream
	git checkout master
	git merge upstream/master
