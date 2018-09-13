# Nanopore long read basecalling, assembly and post-processing workflow

### Before running this workflow, please do the following:

	source activate longread_assembly #activate the environment (every time)

# Inputs
### Alter config.yaml to provide the following:
 * **Sample Name**: Name of sample and output directory

 * **fast5_dirs_list**: list of absolute paths to run/fast5/* subfolders containing .fast5 files

 * **Flowcell**: flowcell code, e.g. FLO-MIN106, passed to basecaller

 * **Kit**: kit code, e.g. SQK-LSK109, passed to basecaller

 * **Genome size**: Estimated genome size, e.g. 50m, passed to canu.

 * **Grid options**: Extra options for execution on the local grid, passed to canu.

Known problems: occasionally fails after calculating DAG. Just re-run snakemake.  This is a problem with dynamic job scheduling, and will hopefully be fixed in a future snakemake update.
