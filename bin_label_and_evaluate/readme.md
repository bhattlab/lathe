# Sequence binning and bin annotation and evaluation workflow
assembly:assembly_links/6753_foo/contigs.fasta,
sample:sr_6753_foo,
reads1:read_links/6753_foo_1.fq.gz,
reads2:read_links/6753_foo_2.fq.gz,
krakendb:~/bhattlab/data/program_indices/kraken/kraken_custom/

### Inputs
# Alter config.yaml to provide the following:
 * **Assembly**: Sequence to bin. Fasta format.

 * **Sample**: names the output directory.

 * **Reads 1, Reads 2**: forward and reverse reads in fastq or fastq.gz format.

 * **Krakendb**: Kraken database with which to classify asssembly contigs.
