# Sequence binning and bin annotation and evaluation workflow

### Before running this workflow, please do the following:

	source activate mgwf #activate the environment
	cd <checkm data directory of your choice>
	wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz #download checkm databases
	tar -zxf checkm_data_2015_01_16.tar.gz
	checkm data setRoot #set the location for checkm data and wait for it to initialize


# Inputs
### Alter config.yaml to provide the following:
 * **Assembly**: Sequence to bin. Fasta format.

 * **Sample**: names the output directory.

 * **Reads 1, Reads 2**: forward and reverse reads in fastq or fastq.gz format.

 * **Krakendb**: Kraken database with which to classify asssembly contigs.
 
 * **Read length**: read length.

Known problems: occasionally fails after binning step. Just re-run snakemake.  This is a problem with dynamic job scheduling, and will hopefully be fixed in a future snakemake update.
