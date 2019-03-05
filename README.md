# Metagenomics Workflows
Workflows for metagenomic sequence data processing and analysis.  Further documentation found in each workflow folder.

# Installing workflows

First, install [miniconda3](https://conda.io/en/latest/miniconda.html)

Then install snakemake.  This can be done with the following.

```
conda install snakemake
snakemake --version #please ensure this is >=5.4.3
```

Next, clone this github directory to some location where it can be stored permanently.  Remember to keep it updated with `git pull`.

```
git clone https://github.com/elimoss/metagenomics_workflows.git
```

Snakemake does not have native support for SLURM. Instructions to enable Snakemake to schedule cluster jobs with SLURM can be found at https://github.com/bhattlab/slurm


# long_read_assembly: Nanopore long read basecalling, assembly and post-processing workflow

### Note: it is highly recommended to run this workflow on a cluster.  

# Inputs
### Alter config.yaml to provide the following:
 * **sample_name**: Name of sample and output directory

 * **fast5_dirs_list**: text file containing a list of absolute paths to run/fast5/* subfolders containing .fast5 files.  A good way to generate this is with `find -maxdepth 2 -mindepth 2 fast5_parent_dir > fodn.txt`

 * **flowcell**: flowcell code, e.g. FLO-MIN106, passed to basecaller

 * **kit**: kit code, e.g. SQK-LSK109, passed to basecaller

 * **genome_size**: Estimated genome size, e.g. 50m, passed to Canu.

 * **singularity**: location (including on the internet) of a singularity image to be used for the workflow.  Don't change this.

 * **short_reads**: location of short reads to be used for Pilon polishing, or empty quotes for long-read polishing.

 * **use_grid**: should Canu execute in distributed mode on a cluster?

 * **grid_options**: Extra options for execution on a cluster

 * **canu_args**: Extra options for Canu

 * **skip_circularization**: Should circularization be omitted from the workflow?


For cluster Canu execution, please note: if set to True, you will need to install Canu in your environment, e.g. `conda install -c conda-forge -c bioconda Canu=1.8` as well as provide any additional required parameters for your job scheduler in the config.yaml file.  When executing on a cluster, Canu will appear to Snakemake to fail, as the first process does not produce an assembly and instead spawns subsequent jobs on the cluster.  Don't worry, just re-run Snakemake when the assembly completes.

To execute this workflow, please run the following.  Please note, you must substitute a parent directory containing all of your data and working directories for `/labs/`.

```
snakemake --use-singularity --singularity-args '--bind /labs/' -s path/to/metagenomics_workflows/long_read_assembly/Snakefile --configfile path/to/modified_config.yaml --restart-times 0 --keep-going
```


# bin_label_and_evaluate

Snakemake workflow for aligning, binning, classifying and evaluating a
metagenomic assembly.

## Before running this workflow, please do the following:

	source activate mgwf #activate the environment
	cd <checkm data directory of your choice>
	wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz #download checkm databases
	tar -zxf checkm_data_2015_01_16.tar.gz
	checkm data setRoot #set the location for checkm data and wait for it to initialize

## Inputs
### Alter config.yaml to provide the following:
 * **Assembly**: Sequence to bin. Fasta format.

 * **Sample**: names the output directory.

 * **Reads 1, Reads 2**: forward and reverse reads in fastq or fastq.gz format.

 * **Krakendb**: Kraken database with which to classify asssembly contigs.

 * **Read length**: read length.

Known problems: occasionally fails after binning step. Just re-run snakemake.  This is a problem with dynamic job scheduling, and will hopefully be fixed in a future snakemake update.


# assembly_comparison_circos
Snakemake workflow for visualizing assemblies of a particular genome across conditions and time points.  Calls out pre-identified sequences, highlights selected contigs.
