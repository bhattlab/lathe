# Sequence binning and bin annotation and evaluation workflow
assembly:assembly_links/6753_foo/contigs.fasta,
sample:sr_6753_foo,
reads1:read_links/6753_foo_1.fq.gz,
reads2:read_links/6753_foo_2.fq.gz,
krakendb:~/bhattlab/data/program_indices/kraken/kraken_custom/

###Before running this workflow, please do the following:

	source activate mgwf #activate the environment
	cd <checkm data directory of your choice>
	wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz #download checkm databases
	tar -zxf checkm_data_2015_01_16.tar.gz
	checkm data setRoot #set the location for checkm data and wait for it to initialize


### Inputs
# Alter config.yaml to provide the following:
 * **Assembly**: Sequence to bin. Fasta format.

 * **Sample**: names the output directory.

 * **Reads 1, Reads 2**: forward and reverse reads in fastq or fastq.gz format.

 * **Krakendb**: Kraken database with which to classify asssembly contigs.
