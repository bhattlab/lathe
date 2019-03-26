'''
Long read assembly and post-processing workflow for assembling genomes from metagenomes.

Author: Eli Moss
'''

import os
import glob

localrules: pilon_ranges, pilon_aggregate_vcf, assemble_final, faidx, extract_bigtigs, circularize_final, polish_final

sample = config['sample_name']
fast5_dirs = [l.strip() for l in open(config['fast5_dirs_list'], 'r').readlines()]
fast5_abspath_run_subfolder = {}
fast5_run_subfolder_abspath = {}

singularity_image = config['singularity']

for d in fast5_dirs:
	run_subfolder = "/".join([d.split("/")[::-1][2], d.split("/")[::-1][0]])
	fast5_abspath_run_subfolder[run_subfolder] = d
	fast5_run_subfolder_abspath[d] = run_subfolder

rule all:
	input:
		#"{sample}/2.polish/{sample}_polished.fa".format(sample = sample),
		"{sample}/5.final/{sample}_final.fa".format(sample = sample), #this must be commented out until after the workflow has reached the first output
		'{sample}/0.basecall/nanoplots/Weighted_LogTransformed_HistogramReadlength.png'.format(sample = sample)

#Basecalling and assembly
#########################################################

rule basecall:
	input: lambda wildcards: fast5_abspath_run_subfolder["/".join([wildcards.run, wildcards.subfolder])]
	output: '{sample}/0.basecall/raw_calls/{run}/{subfolder}/sequencing_summary.txt'
	threads: 4
	resources:
		time=2,
		mem=16
	singularity: singularity_image
	shell:
		"guppy_basecaller --cpu_threads_per_caller {threads} -i {input} -s {sample}/0.basecall/raw_calls/{wildcards.run}/{wildcards.subfolder}/ " +
		"--flowcell {fc} --kit {k}" .format(fc=config['flowcell'], k = config['kit'])

rule basecall_final:
	input: expand('{{sample}}/0.basecall/raw_calls/{foo}/sequencing_summary.txt', foo = fast5_abspath_run_subfolder.keys())
	output: '{sample}/0.basecall/{sample}.fq'
	shell:
		"find {sample}/0.basecall/raw_calls/*/*/*.fastq | xargs cat > {output}"

rule nanoplot:
	input: rules.basecall_final.output
	output: '{sample}/0.basecall/nanoplots/Weighted_LogTransformed_HistogramReadlength.png'
	resources:
		time=4
	threads: 12
	singularity: "docker://quay.io/biocontainers/nanoplot:1.20.0--py35_0"
	shell: "NanoPlot --fastq {input} -t {threads} --color green  -o " + "{s}/0.basecall/nanoplots".format(s = sample)

rule assemble:
	input: rules.basecall_final.output
	output:
		'{sample}/1.assemble/assemble_{genome_size}/{sample}_{genome_size}.contigs.fasta',
		'{sample}/1.assemble/assemble_{genome_size}/{sample}_{genome_size}.correctedReads.fasta.gz'
	threads: 1
	resources:
		mem=100,
		time=80
	singularity: singularity_image if config['usegrid'] != 'True' and config['usegrid'] != True else '' #switch between image and local environment for canu depending on whether cluster execution is required
	shell:
		"canu -p {sample}_{wildcards.genome_size} -d {sample}/1.assemble/assemble_{wildcards.genome_size}/ -nanopore-raw {input} " +
		config['canu_args'] +
		" stopOnReadQuality=false genomeSize={wildcards.genome_size} " +
		"useGrid={grid} gridOptions='{opts}'".format(
			grid = config['usegrid'],
			opts = config['grid_options']
		)


rule misassemblies_detect:
	input:
		"{sample}/{sequence}.fa{sta}", #rules.assemble.output[0],
		"{sample}/{sequence}.fa{sta}.fai", #rules.assemble.output[0] + '.fai',
		"{sample}/{sequence}.fa{sta}.bam",
		"{sample}/{sequence}.fa{sta}.bam.bai"
	output: "{sample}/{sequence}.fa{sta}.misassemblies.tsv" #"{sample}/1.assemble/assemble_{genome_size}/misassemblies/misassemblies.tsv"
	params:
		window_width = 2000,
		min_tig_size = 50000
	resources:
		mem=24,
		time=6
	singularity: "shub://elimoss/lathe:htsbox"
	shell:
		"""
	    bedtools makewindows -g {input[1]} -w {params.window_width} | join - {input[1]} | tr ' ' '\t' | \
	    cut -f1-4 | awk '{{if ($2 > {params.window_width} && $3 < $4 - {params.window_width} && $4 > {params.min_tig_size}) print $0}}' | \
	    xargs -P 16 -l bash -c '
	     htsbox samview {input[2]} $0:$1-$1 -p | \
	     cut -f8,9 | awk "{{if (\$1 < $1 - ({params.window_width}/2) && \$2 > $1 + ({params.window_width}/2)) print \$0}}" | wc -l | \
	     paste <(echo $0) <(echo $1) - ' | awk '{{if ($3 < 2) print $0}}
	    ' > {output}
		"""

rule misassemblies_correct:
	input:
		"{sample}/{sequence}.fa{sta}.misassemblies.tsv", # rules.misassemblies_detect_subasm.output,
		"{sample}/{sequence}.fa{sta}.fai",
		"{sample}/{sequence}.fa{sta}" #rules.assemble.output[0]
	output:
		"{sample}/{sequence}.corrected.fa{sta}"
	shell:
		"""
		cat {input[0]} | grep -v ^# | sort -k1,1 -k2,2g | join - <(cat {input[1]} | sort -k1,1 -k2,2g) | \
		awk '{{
		if ($1 == prev_tig){{
		 print($1,prev_coord,$2)
		 }}
		else{{
			if (prev_len > 0){{
				print(prev_tig,prev_coord,prev_len)
			}}
			print($1,"1",$2)
			}}
		prev_tig = $1
		prev_coord = $2
		prev_len = $4
		}}
		END {{ print(prev_tig,prev_coord,prev_len) }}' | sed "s/\(.*\) \(.*\)\ \(.*\)/\\1:\\2-\\3/g" |  xargs samtools faidx {input[2]} \
		| cut -f1 -d ':' | awk '(/^>/ && s[$0]++){{$0=$0\"_\"s[$0]}}1;' > {output[0]}

		cut -f1 {input[0]} > {sample}/{wildcards.sequence}.tigs.toremove
		grep -vf {sample}/{wildcards.sequence}.tigs.toremove {input[1]} | cut -f1 | xargs samtools faidx {input[2]} >> {output[0]}
		rm {sample}/{wildcards.sequence}.tigs.toremove
		"""

rule merge:
	input: expand("{{sample}}/1.assemble/assemble_{g}/{{sample}}_{g}.contigs.corrected.fasta", g = config['genome_size'].split(","))
	output: "{sample}/1.assemble/{sample}_merged.fasta"
	resources:
		time=6,
		mem=24
	singularity: "shub://elimoss/lathe:quickmerge"
	shell:
		"merge_wrapper.py {input} -ml 10000 -c 5 -hco 10; mv merged_out.fasta {output}"

rule contig_size_filter:
	input:
		rules.merge.output[0],
		rules.merge.output[0] + '.fai'
	output: '{sample}/1.assemble/{sample}_merged_mincontig_{contig_cutoff}.fasta',
	shell: "sort -k2,2gr {input[1]} | awk '{{if ($2 > {wildcards.contig_cutoff}) print $1}}' | xargs samtools faidx {input[0]} > {output}"

def choose_contig_cutoff(wildcards):
	if 'min_contig_size' in config and int(config['min_contig_size'] > 0):
		return(expand("{{sample}}/1.assemble/assemble_{g}/{{sample}}_{g}.contigs.mincontig_{contig_cutoff}.fasta",
			g = config['genome_size'].split(","),
			contig_cutoff = config['min_contig_size']
			))
	else:
		return(rules.merge.output)

rule assemble_final:
	input: choose_contig_cutoff
	output: "{sample}/1.assemble/{sample}_raw_assembly.fa"
	shell:
		"cat {input} | cut -f1 -d '_' | fold -w 120 | awk '(/^>/ && s[$0]++){{$0=$0\"_\"s[$0]}}1;' > {output}"


		#"cp {input[0]} {output}" #request all specified assemblies, but use the first genome size listed in the config
		#"cat {input} | cut -f1 -d '_' | fold -w 120 | awk '(/^>/ && s[$0]++){{$0=$0\"_\"s[$0]}}1;' > {output}" #sed 's/\.tig/\\ttig/g' | tr '_' '\\t' | " +
		#"cut -f1,9,17 | tr '\\t' '_'

#Polishing
#########################################################

def choose_polish(wildcards):
	#skip_polish_or_not(wildcards)
	if 'skip_polishing' in config and (config['skip_polishing'] == True or config['skip_polishing'] == 'True'):
		return(rules.assemble_final.output)
	elif config['short_reads'] != '':
		return(rules.pilon_consensus.output)
	else:
		return(rules.medaka.output)

def choose_pilon_input():
	if 'polish_both' in config and config['polish_both'] == True:
		return(rules.medaka.output)
	else:
		return(rules.assemble_final.output)

def get_racon_input(wildcards):
	#this method will choose a previous iteration of racon or the original assembly, depending on the value of the iteration wildcard
	if int(wildcards.iteration) == 1:
		return(rules.assemble_final.output[0] + '.paf', rules.assemble_final.output[0])
	else:
		result = "{sample}/2.polish/racon/{sample}_racon_{iter}.fa".format(sample = wildcards.sample, iter = str(int(wildcards.iteration) - 1))
		return(result + '.paf', result)

rule racon:
	input:
		rules.basecall_final.output,
		get_racon_input
	output:
		'{sample}/2.polish/racon/{sample}_racon_{iteration}.fa'
	threads: 16
	singularity: singularity_image
	resources:
		mem=48,
		time=8
	shell:
		"""
		racon -m 8 -x -6 -g -8 -w 500 -t {threads} {input} > {output}
		"""

rule medaka:
	input:
		rules.basecall_final.output,
		'{sample}/2.polish/racon/{sample}_racon_4.fa' # request four iterations of racon, as specified by medaka docs
	output: '{sample}/2.polish/medaka/{sample}_medaka.fa'
	threads: 16
	resources:
		mem=32,
		time=100
	singularity: singularity_image
	shell:
		"""
		medaka_consensus -i {input[0]} -d {input[1]} -o {sample}/2.polish/medaka -t {threads} -m r941_flip213
		cut -f1 -d ':' {sample}/2.polish/medaka/consensus.fasta > {output}
		"""

rule align_short_reads:
	input:
		choose_pilon_input(),
		config['short_reads'].split(',')
	output: "{sample}/2.polish/pilon/short_reads.bam"
	threads: 16
	params:
		reads = config['short_reads'].split(',')
	resources:
		mem=100,
		time=24
	singularity: singularity_image
	shell:
		"bwa index {input[0]}; bwa mem -t {threads} {input} | samtools sort --threads {threads} > {output}"

checkpoint pilon_ranges:
	input:
		choose_pilon_input(),
		choose_pilon_input()[0] + '.fai'
	output: directory('{sample}/2.polish/pilon/ranges')
	singularity: singularity_image
	shell:
		"""
		mkdir {output}
		bedtools makewindows -w 100000 -g {input[1]} | awk '{{print $1,\":\", $2+ 1, \"-\", $3}}'  | tr -d ' ' |
		xargs -n 1 -I foo touch {sample}/2.polish/pilon/ranges/foo
		"""

rule pilon_subsetrun:
	input:
		choose_pilon_input(),
		rules.align_short_reads.output,
		rules.align_short_reads.output[0] + '.bai',
		'{sample}/2.polish/pilon/ranges/{range}', #exclude hidden files like .snakemake_timestamp
	output:
		'{sample}/2.polish/pilon/sub_runs/{range}/{sample}_{range}.vcf.gz'
	resources:
		time=4,
		mem=32
	singularity: singularity_image
	params:
		java_mem = 30,
		bam = "{sample}/2.polish/pilon/sub_runs/{range}/{sample}_{range}.bam",
		fa = "{sample}/2.polish/pilon/sub_runs/{range}/{sample}_{range}.fa",
		subrun_folder = "{sample}/2.polish/pilon/sub_runs/{range}",
		target_coverage = 50
	shell:
		"""
		# set env var $i to be the smallest read subset decimal (in increments of 0.1, with a couple very low values thrown in, too)
    	# sufficient to generate at least 40x coverage depth of the target sequence, or
		# 1 if 40x coverage cannot be achieved with the available read data

		for i in 0.01 0.05 $(seq 0.1 0.1 1);
		do
		   cov=$(samtools view {input[1]} -s $i -h {wildcards.range} | samtools depth - | cut -f3 | awk '{{sum+=$1}}END{{print sum/(NR+1)}}')
		   if [ $(echo $cov'>'{params.target_coverage}|bc) -eq 1 ]
		   then
		       break
		   fi
		done
		echo Using $i x of total coverage;

		samtools view -h -O BAM -s $i {input[1]} {wildcards.range} > {params.bam}
		samtools index {params.bam}
		samtools faidx {input[0]} $(echo {wildcards.range}| cut -f1 -d ':') | cut -f1 -d ':' > {params.fa}
		java -Xmx{params.java_mem}G -jar $(which pilon | sed 's/\/pilon//g')/../share/pilon*/pilon*.jar \
			--genome {params.fa} \
			--unpaired {params.bam} --output {sample}_{wildcards.range} --outdir {params.subrun_folder} \
			--vcf --nostrays --mindepth 3
		bgzip {params.subrun_folder}/{sample}_{wildcards.range}.vcf
		tabix -fp vcf {params.subrun_folder}/{sample}_{wildcards.range}.vcf.gz
		"""

def aggregate_pilon_subsetruns(wildcards):
	checkpoint_output = checkpoints.pilon_ranges.get(**wildcards).output[0]
	result = expand(rules.pilon_subsetrun.output,
		sample=wildcards.sample,
		range=glob_wildcards(os.path.join(checkpoint_output, '{range,tig.+}')).range)
	return(result)

rule pilon_aggregate_vcf:
	input:
		aggregate_pilon_subsetruns
	output:
		'{sample}/2.polish/pilon/corrections.vcf.gz'
	resources:
		time=4,
		mem=8
	singularity: singularity_image
	shell:
		"""
		#workaround!  Snakemake was causing the vcf's to appear newer than the indices, which tabix didn't like
		touch {sample}/2.polish/pilon/sub_runs/*/*.vcf.gz.tbi

		#get properly sorted intervals
		ls {sample}/2.polish/pilon/ranges/ | tr ':-' '\t' | sort -k1,1 -k2,2g | awk '{{print $1,":",$2,"-",$3}}' | tr -d ' ' > sorted_ranges.tmp

		#get header
		(zcat {sample}/2.polish/pilon/sub_runs/*/*.vcf.gz | head -1000 | grep ^#

		#get corrections within each range (omitting DUP records, which bcftools can't understand)
		cat sorted_ranges.tmp | xargs -n 1 -I foo sh -c "
		tabix {sample}/2.polish/pilon/sub_runs/foo/{sample}_foo.vcf.gz foo | grep -v '0/0' | grep -v DUP" |

		#sort by position
		sort -k1,1 -k2,2g) |

		#compress and store
		bgzip > {output} || true

		#index
		tabix -p vcf {output}
		"""

rule pilon_consensus:
	input:
		choose_pilon_input(),
		rules.pilon_aggregate_vcf.output
	output:
		"{sample}/2.polish/pilon/{sample}_pilon.fa"
	singularity: singularity_image
	shell:
		"""
		bcftools consensus -f {input} -o {output}
		"""
			#&& rm -rf {sample}/2.polish/pilon/sub_runs {sample}/2.polish/pilon/ranges

rule polish_final:
	input: choose_polish
	output: "{sample}/2.polish/{sample}_polished.fa"
	shell:
		"""
		cut -f1 -d ':' {input} > {output}
		"""

#Circularization
#########################################################

rule circularize_mapreads:
	input:
		rules.polish_final.output,
		'{sample}/1.assemble/assemble_{g}/{sample}_{g}.correctedReads.fasta.gz'.format(
			g = config['genome_size'].split(",")[0],
			sample = sample)
	output:
		"{sample}/3.circularization/0.aligned_corrected.bam",
		"{sample}/3.circularization/0.aligned_corrected.bam.bai",
	threads: 8
	resources:
		time=6,
		mem=32
	singularity: singularity_image
	shell:
		"""
		minimap2 {input} -ax map-ont -t {threads} | samtools sort --threads {threads} > {output[0]}
		samtools index {output[0]}
		"""

checkpoint extract_bigtigs:
	input:
		rules.polish_final.output,
		rules.polish_final.output[0] + ".fai",
	output: directory("{sample}/3.circularization/1.candidate_genomes/")
	singularity: singularity_image
	params:
		min_size = 1700000
	shell:
		"""
		#mkdir {output}
		cat {input[1]} | awk '{{if ($2 > {params.min_size}) print $1}}' | xargs -n 1 -I foo sh -c "
			samtools faidx {input[0]} foo > {sample}/3.circularization/1.candidate_genomes/foo.fa
		"
		"""

rule circularize_bam2reads:
	input:
		rules.circularize_mapreads.output,
		"{sample}/3.circularization/1.candidate_genomes/{tig}.fa"
	output:
		"{sample}/3.circularization/2.circularization/spanning_tig_circularization/{tig}/{tig}_terminal_reads.fq.gz"
	singularity: singularity_image
	shell:
		"""
		(samtools idxstats {sample}/3.circularization/0.aligned_corrected.bam | grep {wildcards.tig} | awk '{{if ($2 > 50000) print $1, ":", $2-50000, "-", $2; else print $1, ":", 1, "-", $2 }}' | tr -d ' ';
		 samtools idxstats {sample}/3.circularization/0.aligned_corrected.bam | grep {wildcards.tig} | awk '{{if ($2 > 50000) print $1, ":", 1, "-", 50000; else print $1, ":", 1, "-", $2 }}' | tr -d ' ') |
		xargs -I foo sh -c 'samtools view -h {sample}/3.circularization/0.aligned_corrected.bam foo | samtools fastq - || true' | bgzip > {output}
		"""

rule circularize_assemble:
	input:
		rules.circularize_bam2reads.output
	output: "{sample}/3.circularization/2.circularization/spanning_tig_circularization/{tig}/{tig}.contigs.fasta"
	params:
		directory="{sample}/3.circularization/2.circularization/spanning_tig_circularization/{tig}",
	singularity: singularity_image
	resources:
		time=12,
		mem=50
	threads: 8
	shell:
		"""
		canu -useGrid=False -assemble -p {wildcards.tig} -d {params.directory}  \
		-nanopore-corrected {input} genomeSize=100000
		"""

		# """
		# canu -useGrid=False -assemble -p {params.prefix} -d {params.directory}  \
		# -nanopore-corrected {input} genomeSize=100000 gridOptions='{params.gridopts}' \
		# batMemory=128 batThreads=1
		# """

rule circularize_spantig_pre:
	input:
		"{sample}/3.circularization/1.candidate_genomes/{tig}.fa", #rules.extract_bigtigs.output
		rules.circularize_assemble.output,
		rules.circularize_assemble.output[0] + '.fai'
	output:
		"{sample}/3.circularization/2.circularization/spanning_tig_circularization/{tig}/potential_circularization_alignments.tsv"
	singularity: singularity_image
	params:
		directory = "{sample}/3.circularization/2.circularization/spanning_tig_circularization/{tig}",
		prefix="spanning_tigs_to_ref"
	threads: 4
	resources:
		time=4,
		mem=16
	shell:
		"""
		nucmer -b 5000 {input[0]} {input[1]} -p {params.directory}/{params.prefix} #-t {threads} #reverted nucmer back down to 3, no more multithreading :(

		delta-filter -q {params.directory}/{params.prefix}.delta > {params.directory}/{params.prefix}.filt.delta

		show-coords -Tq {params.directory}/{params.prefix}.filt.delta | cut -f8,9 | sed '1,3d' | sort | \
		uniq -c | tr -s ' ' '\\t' | cut -f2-99 | grep -v ^1 | cut -f2,3 > {params.directory}/potential_circularizations.tsv || true

		show-coords -Tql {params.directory}/{params.prefix}.filt.delta | grep -f {params.directory}/potential_circularizations.tsv | cat > {output} || true
		"""

rule circularize_spantig:
	input: rules.circularize_spantig_pre.output
	output: "{sample}/3.circularization/2.circularization/spanning_tig_circularization/{tig}/contig_spanned.txt"
	params:
		margin=10000
	script:
		"scripts/spancircle.py"

rule circularize_span_trim:
	input:
		rules.circularize_spantig.output,
		rules.extract_bigtigs.output[0] + '{tig}.fa',
		rules.extract_bigtigs.output[0] + '{tig}.fa.fai',
		rules.circularize_assemble.output,
		rules.circularize_assemble.output[0] + '.fai'
	output:
		"{sample}/3.circularization/3.circular_sequences/sh/{tig}_span_trim.sh"
	params:
		outfa = "{sample}/3.circularization/3.circular_sequences/{tig}_spanned.fa"
	run:
		span_out = open(input[0], 'r').readlines()
		cmd = ''
		if span_out == ['done\n'] or span_out[0].strip() == 'no circularizations': #no circularization occurred
			print('Nothing to do')
		else:
			trim = span_out[0].strip()
			trim_cmd = 'samtools faidx ' + input[1] + ' ' + trim + " > " + params.outfa + "\n"
			cmd += trim_cmd

			if len(span_out) == 3:
				extend = span_out[1].strip()
				extend_cmd = 'samtools faidx ' + input[3] + ' ' + extend + " | grep -v '>'" + " >> " + params.outfa + "\n"
	#			print(extend_cmd)
				cmd += extend_cmd

		open(output[0], 'w').write(cmd + '\n')

def aggregate_span_trim(wildcards):
	checkpoint_output = checkpoints.extract_bigtigs.get(**wildcards).output[0]
	result = expand(rules.circularize_span_trim.output, #"{sample}/3.circularization/3.circular_sequences/sh/{tig}_span_trim.sh",
		sample=wildcards.sample,
		tig=glob_wildcards(os.path.join(checkpoint_output, '{tig}.fa')).tig)
	return(result)

rule circularize_overcirc:
	input:
		"{sample}/3.circularization/1.candidate_genomes/{tig}.fa"
	output: "{sample}/3.circularization/2.circularization/overcircularized/overcirc_{tig}.txt"
	params:
		delta = '{sample}/3.circularization/2.circularization/overcircularized/{tig}'
	threads: 8
	singularity: singularity_image
	script:
		"scripts/encircle.py"

rule circularize_overcirc_trim:
	input:
		rules.circularize_overcirc.output,
		rules.extract_bigtigs.output[0] + "{tig}.fa",
		rules.extract_bigtigs.output[0] + '{tig}.fa.fai',
		rules.circularize_assemble.output,
		rules.circularize_assemble.output[0] + '.fai'
	output:
		"{sample}/3.circularization/3.circular_sequences/sh/{tig}_span_overcirc.sh"
	params:
		outfa = "{sample}/3.circularization/3.circular_sequences/{tig}_overcirc.fa"
	run:
		span_out = open(input[0], 'r').readlines()
		cmd = ''
		if span_out == ['done\n']: #no circularization occurred
			print('Nothing to do')
		else:
			trim = span_out[0].strip()
			trim_cmd = 'samtools faidx ' + input[1] + ' ' + trim + " > " + params.outfa + "\n"
			cmd += trim_cmd

		open(output[0], 'w').write(cmd + '\n')

def aggregate_overcirc_trim(wildcards):
	checkpoint_output = checkpoints.extract_bigtigs.get(**wildcards).output[0]
	result = expand(rules.circularize_overcirc_trim.output, #"{sample}/3.circularization/3.circular_sequences/sh/{tig}_span_overcirc.sh",
		sample=wildcards.sample,
		tig=glob_wildcards(os.path.join(checkpoint_output, '{tig}.fa')).tig)
	return(result)

rule circularize_final:
	input:
		rules.polish_final.output,
		rules.polish_final.output[0] + '.fai',
		aggregate_overcirc_trim,
		aggregate_span_trim
	output:
		'{sample}/3.circularization/4.{sample}_circularized.fasta'
	threads: 1
	singularity: singularity_image
	shell:
		"""
		cat {sample}/3.circularization/3.circular_sequences/sh/* | bash
		ls {sample}/3.circularization/3.circular_sequences/ | grep .fa$ | cut -f1 -d '_' > circs.tmp || true
		(cat {input[1]} | grep -vf circs.tmp |
		cut -f1 | xargs samtools faidx {input[0]}; cat {sample}/3.circularization/3.circular_sequences/*.fa) |
		sed 's/\([ACTG]\)\\n/\1/g' | fold -w 120 | cut -f1 -d ':' > {output}
		#rm circs.tmp
		"""

#Misassembly detection
#########################################################

def skip_circularization_or_not():
	if config['skip_circularization'] == 'True' or config['skip_circularization'] == True:
		return(rules.polish_final.output)
	else:
		return(rules.circularize_final.output)

rule final:
	input: skip_circularization_or_not()[0].replace('.fasta', '.corrected.fasta') #perform one last round of misassembly breakage
	output: "{sample}/5.final/{sample}_final.fa"
	shell: "cp {input} {output}"

#Utility functions
#########################################################

rule align_paf:
	input:
		'{ref}.f{asta}',
		'{sample}/0.basecall/{sample}.fq'.format(sample = sample)
	output:
		"{ref}.f{asta}.paf"
	threads: 8
	singularity: singularity_image
	shell:
		"minimap2 -t {threads} -x map-ont {input} > {output}"

rule align_bam:
	input:
		'{ref}.f{asta}', #the asta bit makes this work for .fa and .fasta files
		'{sample}/0.basecall/{sample}.fq'.format(sample = sample)
	output:
		"{ref}.f{asta}.bam"
	threads: 16
	resources:
		time=6,
		mem=16
	singularity: singularity_image
	shell:
		"minimap2 -t {threads} -ax map-ont {input} | samtools sort --threads {threads} > {output}"

rule bam_idx:
	input:
		'{some}.bam'
	output:
		'{some}.bam.bai'
	singularity: singularity_image
	shell:
		"samtools index {input}"

rule faidx:
	input: '{something}.f{asta}'
	output: '{something}.f{asta}.fai'
	singularity: singularity_image
	shell: "samtools faidx {input}"
