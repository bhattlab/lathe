# Nanopore long read basecalling, assembly and post-processing workflow

## To Run

To execute this workflow, please run the following.  Please note, you must substitute a parent directory containing all of your data for `/labs/`.  Also, for cluster Canu execution, see note for `use_grid` below.

```
snakemake --use-singularity --singularity-args '--bind /labs/' -s path/to/metagenomics_workflows/long_read_assembly/Snakefile --configfile path/to/config.yaml
```

### Note: it is highly recommended to run this workflow on a cluster.  

# Inputs
### Alter config.yaml to provide the following:
 * **sample_name**: Name of sample and output directory

 * **fast5_dirs_list**: textual list of absolute paths to run/fast5/* subfolders containing .fast5 files

 * **fast5_parent_dir**: parent to the folders in the above file

 * **flowcell**: flowcell code, e.g. FLO-MIN106, passed to basecaller

 * **kit**: kit code, e.g. SQK-LSK109, passed to basecaller

 * **genome_size**: Estimated genome size, e.g. 50m, passed to Canu.

 * **singularity**: location (including on the internet) of a singularity image to be used for the workflow.

 * **short_reads**: location of short reads to be used for Pilon polishing, or empty quotes for long-read polishing.

 * **use_grid**: should Canu execute in distributed mode on a cluster?  Note: if True, you will need to install Canu manually, e.g. `conda install -c conda-forge -c bioconda Canu` as well as provide any additional required parameters for your job scheduler below.  When executing on a cluster, Canu will appear to Snakemake to fail, as the first process does not produce an assembly, but instead spawns a subsequent job on the cluster.  Don't worry, just re-run Snakemake when the assembly eventually completes.  You may need to add --cleanup-metadata <assembly> before Snakemake will continue.

 * **grid_options**: Extra options for execution on a cluster

 * **canu_args**: Extra options for Canu
