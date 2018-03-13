#!/usr/bin/env python

localrules: bwa_index

rule targets:
	input:
		config['reference'],
		config['light_contigs'],
		config['highlight_sequences'],
		config['assemblies']

rule bwa_index_setup:
	input:
		config['reference']
	output:
		"{input}/idx/{input}.fa"
	resources:
		mem=1,
		time=1
	threads: 1
	shell:
		"cp {input} {output}"


rule bwa_index:
	input:
		rules.bwa_index_setup.output
	output:
		"{input}.amb",
		"{input}.ann",
		"{input}.bwt",
		"{input}.pac",
		"{input}.sa",
	log:
		"log/bwa_index.log".format(samp=config['sample'])
	resources:
		mem=8,
		time=24
	threads: 1
	params:
		prefix="{rules.bwa_index_setup.output}",
		algorithm="bwtsw"
	wrapper:
		"0.22.0/bio/bwa/index"

#16S (highlight sequence) alignments
rule highlight_alignments:
	input:
		rules.bwa_index.output,
		config['highlight_sequences']
	output:
		"highlight_alignments/highlights


#contig alignments
rule assembly_alignments:
#nucmer ../../$REF ../../$TRUSEQ
#delta-filter -1 out.delta > $(echo $ORGANISM)_truseq.delta.filt
#rm out.delta
#show-coords $(echo $ORGANISM)_truseq.delta.filt | tr '|' ' ' | tr -s ' ' '\t' > $(echo $ORGANISM)_truseq.coords
#/srv/gsfs0/projects/bhatt/moss/projects/10x/tenex/coord_filter $(echo $ORGANISM)_truseq.coords 5000 | cut -f1,2,8,9 | tr ' ' '\t' | awk '{print $3, $1, $2, $4}' | tr ' ' '\t' > $(echo $ORGANISM)_truseq.bed


rule circos:
	input:
		config['light_contigs'],
		rules.
	script:
		scripts/circos.R
