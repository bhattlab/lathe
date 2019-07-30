'''
Long read assembly and post-processing workflow for assembling genomes from metagenomes.

Author: Eli Moss
'''

import os
import glob
from pathlib import Path
import git
import snakemake

localrules: no_merge, basecall_staging, pilon_ranges, pilon_aggregate_vcf, assemble_final, faidx, extract_bigtigs, circularize_final, polish_final

#set up singularity image used for the bulk of rules
singularity_image = "shub://elimoss/lathe:longread"

#extract samplename from config
sample = config['sample_name']

#find fast5 files. These are expected to be below the directory provided in the config. All fast5's below this directory will be processed.
fast5_files = glob.glob(os.path.join(config['fast5_directory'], '**', '*.fast5'), recursive=True) #list provided directory recursively
#create a dictionary which relates fast5 basenames to the full path of each file
fast5_basename_to_path = {}
for f in fast5_files:
	fast5_basename_to_path[os.path.splitext(os.path.basename(f))[0]] = f

#perform a check on the Lathe git repo and exit if not up to date
onstart:
	import git
	print("Checking for updates or modifications")
	repo_dir = os.path.dirname(workflow.snakefile)
	print(repo_dir)
	repo = git.Repo(repo_dir)
	assert not repo.bare
	repo_git = repo.git
	stat = repo_git.diff('origin/master')
	print([stat])
	if stat != "":
		print('WARNING: Differences to latest version detected. Please reset changes and/or pull repo.')
	#elif stat.split("\n")[1]
	sys.exit("Test complete")


rule all:
	input:
		"{sample}/5.final/{sample}_final.fa".format(sample = sample), #request final output
		'{sample}/0.basecall/nanoplots/Weighted_LogTransformed_HistogramReadlength.png'.format(sample = sample) #request QC data

#Basecalling and assembly
#########################################################
rule basecall_staging:
	input: lambda wildcards: fast5_basename_to_path[wildcards.fast5_basename]
	output: '{sample}/0.basecall/data_links/{fast5_basename}/{fast5_basename}.fast5'
	shell:
		"ln -s {input} {output}"

rule basecall:
	#Call bases from raw fast5 data with guppy basecaller. Called once per input folder of fast5 files.
	input: "{sample}/0.basecall/data_links/{fast5_basename}/{fast5_basename}.fast5"
	output: '{sample}/0.basecall/raw_calls/{fast5_basename}/sequencing_summary.txt'
	threads: 4
	resources:
		time=12,
		mem=16
	params:
		in_dir = "{sample}/0.basecall/data_links/{fast5_basename}/",
		out_dir = '{sample}/0.basecall/raw_calls/{fast5_basename}/'
	singularity: singularity_image
	shell:
		"guppy_basecaller --cpu_threads_per_caller {threads} -i {params.in_dir} -s {params.out_dir} " +
		"--flowcell {fc} --kit {k}" .format(fc=config['flowcell'], k = config['kit'])

rule basecall_final:
	#Collate called bases into a single fastq file
	input: expand('{{sample}}/0.basecall/raw_calls/{foo}/sequencing_summary.txt', foo = [os.path.basename(f).split('.')[0] for f in fast5_files])
	output: '{sample}/0.basecall/{sample}.fq'
	shell:
		"find {sample}/0.basecall/raw_calls/*/*.fastq | xargs cat > {output}"

rule nanoplot:
	#Run nanoplot on the collated .fq file containing basecalled reads. This produces helpful stats
	#such as read length and basecall qualities
	input: rules.basecall_final.output
	output: '{sample}/0.basecall/nanoplots/Weighted_LogTransformed_HistogramReadlength.png'
	resources:
		time=4
	threads: 12
	singularity: "docker://quay.io/biocontainers/nanoplot:1.20.0--py35_0"
	shell: "NanoPlot --fastq {input} -t {threads} --color green  -o " + "{s}/0.basecall/nanoplots".format(s = sample)

rule assemble_canu:
	#Run canu. This can be run either in distributed fashion on the cluster or in a single job. If run in distributed fashion (usedgrid is set to 'True')
	#then canu must be installed in the user's environment. In this case, the singualarity image is not used. Canu arguments are passed through from
	#the config within the 'canu_args' key. 'gridOptions' are also passed through from the config and instruct Canu on any additional parameters
	#needed to run on the cluster.
	input: rules.basecall_final.output
	output:
		'{sample}/1.assemble/assemble_{genome_size}/{sample}_{genome_size}.contigs.fasta',
		#'{sample}/1.assemble/assemble_{genome_size}/{sample}_{genome_size}.correctedReads.fasta.gz'
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

rule assemble_flye:
	input: rules.basecall_final.output
	output:
		"{sample}/1.assemble/assemble_{genome_size}/assembly.fasta"
	threads: 16
	resources:
		mem=100,
		time=100
	singularity: "docker://quay.io/biocontainers/flye:2.4.2--py27he860b03_0"
	shell:
		"flye -t {threads} --nano-raw {input} --meta -o {wildcards.sample}/1.assemble/assemble_{wildcards.genome_size}/ -g {wildcards.genome_size}"

rule misassemblies_detect:
	#Detect regions in the assembly which are not spanned by more than a single long read, indicating likely misassemblies.
	#This is done by first generating a list of all non-overlapping windows of specified length tiled across the assembly. For each,
	#aligned reads are retrieved and the start and endpoints of their alignments are determined. These are filtered for those alignments
	#which span completely across the window. This approach is premised on misasssemblies causing clipped read alignments which can then be
	#identified as alignments which do not span the misassembled locus.
	input:
		"{sample}/{sequence}.fa{sta}",
		"{sample}/{sequence}.fa{sta}.fai",
		"{sample}/{sequence}.fa{sta}.bam",
		"{sample}/{sequence}.fa{sta}.bam.bai"
	output: "{sample}/{sequence}.fa{sta}.misassemblies.tsv" #"{sample}/1.assemble/assemble_{genome_size}/misassemblies/misassemblies.tsv"
	params:
		window_width = 2000,
		min_tig_size = 50000
	resources:
		mem=24,
		time=24
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
	#Given a list of misassembled loci identified in the rule misassemblies_detect, break the assembly at each location.
	#This must accommodate instances where multiple breaks occur in a single contig, producing multiple contigs with successive
	#numbers appended to the name of the original unbroken contig.
	input:
		"{sample}/{sequence}.fa{sta}.misassemblies.tsv",
		"{sample}/{sequence}.fa{sta}.fai",
		"{sample}/{sequence}.fa{sta}"
	output:
		"{sample}/{sequence}.corrected.fa{sta}"
	singularity: singularity_image
	shell:
		"""
		if [ -s {input[0]} ] #if the input is nonempty...
		then
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
				END {{ print(prev_tig,prev_coord,prev_len)
			}}' | sed "s/\(.*\) \(.*\)\ \(.*\)/\\1:\\2-\\3/g" |  xargs samtools faidx {input[2]} \
			| cut -f1 -d ':' | awk '(/^>/ && s[$0]++){{$0=$0\"_\"s[$0]}}1;' > {output[0]}

			#cut -f1 {input[0]} > {sample}/{wildcards.sequence}.tigs.toremove
			#grep -vf {sample}/{wildcards.sequence}.tigs.toremove {input[1]} | cut -f1 | xargs samtools faidx {input[2]} >> {output[0]}
			#rm {sample}/{wildcards.sequence}.tigs.toremove
		else
			cp {input[2]} {output}
		fi
		"""

def choose_assembler():
	if config['assembler'] == 'canu':
		return(expand('{{sample}}/1.assemble/assemble_{genome_size}/{sample}_{genome_size}.contigs.fasta', genome_size = config['genome_size'].split(",")))
	elif config['assembler'] == 'flye':
		return(expand("{{sample}}/1.assemble/assemble_{g}/assembly.corrected.fasta", g = config['genome_size'].split(",")))

rule merge:
	#Conservatively merge the two subassemblies
	input: choose_assembler()
	output: "{sample}/1.assemble/{sample}_merged.fasta"
	resources:
		time=6,
		mem=24
	singularity: "shub://elimoss/lathe:quickmerge"
	shell:
		"merge_wrapper.py {input} -ml 10000 -c 5 -hco 10; mv merged_out.fasta {output}"

rule no_merge:
	input: choose_assembler() # "{{sample}}/1.assemble/assemble_{g}/{{sample}}_{g}.contigs.corrected.fasta".format(g = config['genome_size'])
	output: "{sample}/1.assemble/{sample}_nomerge.fasta"
	shell:
		"ln -s ../../{input} {output}"

def choose_merge():
	if len(config['genome_size'].split(",")) == 1:
		return(rules.no_merge.output)
	else:
		return(rules.merge.output)

rule contig_size_filter:
	#Optional convenience function. Filter out contigs below a certain size. Not used by default,
	#but requested by users with specific use cases.
	input:
		choose_merge()[0],
		choose_merge()[0] + '.fai'
	output: '{sample}/1.assemble/{sample}_merged_mincontig_{contig_cutoff}.fasta',
	shell: "sort -k2,2gr {input[1]} | awk '{{if ($2 > {wildcards.contig_cutoff}) print $1}}' | xargs samtools faidx {input[0]} > {output}"

def choose_contig_cutoff(wildcards):
	#If a contig cutoff is specified in the config, then perform a size cutoff on the assembled contigs by
	#requesting the output of the relevant rule.
	if 'min_contig_size' in config and int(config['min_contig_size'] > 0):
		return(rules.contig_size_filter.output)
	else:
		return(choose_merge()[0])

rule assemble_final:
	#Request the final output of the assembly phase of the pipeline. Some minor cleanup and a uniquefication of the contig names
	#are performed.
	input: choose_contig_cutoff
	output: "{sample}/1.assemble/{sample}_raw_assembly.fa"
	shell:
		"cat {input} | cut -f1 -d '_' | fold -w 120 | awk '(/^>/ && s[$0]++){{$0=$0\"_\"s[$0]}}1;' > {output}"


#Polishing
#########################################################

def choose_polish(wildcards):
	#Determine which type of polishing should be performed. Options are short read (requires a dataset), long read,
	#both, and none. These are controlled by the skip_polishing and short_read variables in the config. Long read
	#is the default. Short read is performed instead if short reads are provided. No polishing is specified with the
	#skip_polishing config variable.
	if 'skip_polishing' in config and (config['skip_polishing'] == True or config['skip_polishing'] == 'True'):
		return(rules.assemble_final.output)
	elif config['short_reads'] != '':
		return(rules.pilon_consensus.output)
	else:
		return(rules.medaka_aggregate.output)

def choose_pilon_input():
	#allows consensus refinement with both long read and short read polishing if the polish_both variable is given the value
	#True in the config file.
	if 'polish_both' in config and config['polish_both'] == True:
		return(rules.medaka_aggregate.output)
	else:
		return(rules.assemble_final.output)

def get_racon_input(wildcards):
	#this method will choose a previous iteration of racon or the unpolished assembly, depending on the value of the iteration wildcard
	if int(wildcards.iteration) == 1:
		return(rules.assemble_final.output[0] + '.paf', rules.assemble_final.output[0])
	else:
		result = "{sample}/2.polish/racon/{sample}_racon_{iter}.fa".format(sample = wildcards.sample, iter = str(int(wildcards.iteration) - 1))
		return(result + '.paf', result)

rule racon:
	#Perform long read polishing with racon. This rule is used in multiple iterations.
	input:
		rules.basecall_final.output,
		get_racon_input
	output:
		'{sample}/2.polish/racon/{sample}_racon_{iteration}.fa'
	threads: 12
	singularity: "docker://quay.io/biocontainers/racon:1.3.2--he941832_0"
	resources:
		mem=48,
		time=24
	shell:
		"""
		racon -m 8 -x -6 -g -8 -w 500 -t {threads} {input} > {output}
		"""

checkpoint medaka_ranges:
	input:
		'{sample}/2.polish/racon/{sample}_racon_4.fa', # request four iterations of racon, as specified by medaka docs
		'{sample}/2.polish/racon/{sample}_racon_4.fa.fai' # request four iterations of racon, as specified by medaka docs
	output:
		directory('{sample}/2.polish/medaka/ranges'),
		directory('{sample}/2.polish/medaka/sub_fa')
	singularity: singularity_image
	shell:
		"""
		mkdir -p {output}
		cut -f1 {input[1]} | xargs -n 1 -I foo sh -c 'touch {output[0]}/foo; samtools faidx {input[0]} foo > {output[1]}/foo.fa'
		"""

rule medaka_consensus:
	input:
		"{sample}/2.polish/medaka/ranges/{range}",
		"{sample}/2.polish/racon/{sample}_racon_4.fa.bam",
		"{sample}/2.polish/racon/{sample}_racon_4.fa.bam.bai"
	output:
		"{sample}/2.polish/medaka/subruns/{range}_probs.hdf"
	singularity: singularity_image
	threads: 1
	resources:
		mem=lambda wildcards, attempt: attempt * 8,
		time=lambda wildcards, attempt: attempt * 12
	shell:
		"""
		medaka consensus {input[1]} {output} --model r941_flip213 --threads {threads} --regions {wildcards.range}
		"""

rule medaka_stitch:
	#Perform long read polishing with medaka. This is used following multiple iterations of racon polishing. The
	#number of iterations is specified in the input filename below.
	input:
		rules.medaka_consensus.output,
		"{sample}/2.polish/medaka/sub_fa/{range}.fa"
	output: '{sample}/2.polish/medaka/subruns/{range}_medaka.fa'
	threads: 1
	resources:
		mem=16,
		time=12
	singularity: singularity_image
	shell:
		"""
		medaka stitch {input[0]} {input[1]} {output}
		"""

#		medaka_consensus -i {input[0]} -d {input[1]} -o {sample}/2.polish/medaka -t {threads} -m r941_flip213
#		cut -f1 -d ':' {sample}/2.polish/medaka/consensus.fasta > {output}
#		"""

def aggregate_medaka_subsetruns(wildcards):
	#The individually polished sub-regions of the overall assembly must be collected into a consensus assembly. This rule
	#gathers all the outputs that have been generated for input into the aggregation rule.
	checkpoint_output = checkpoints.medaka_ranges.get(**wildcards).output[0]
	result = expand(rules.medaka_stitch.output,
		sample=wildcards.sample,
		range=glob_wildcards(os.path.join(checkpoint_output, '{range,tig.+}')).range)
	return(result)


rule medaka_aggregate:
	input: aggregate_medaka_subsetruns
	output: "{sample}/2.polish/medaka/{sample}_medaka.fa"
	params:
		indir = "{sample}/2.polish/medaka/sub_fa/*"
	shell:
		"""
		cat {params.indir} | cut -f1 -d ':' > {output}
		"""

rule align_short_reads:
	#Align short reads to the unpolished or longread-polished assembly, depending on the output of
	#the choose_pilon_input rule above.
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
	#Generate a large collection of intervals within the assembly for individual parallelized polishing with Pilon.
	#This is done with bedtools makewindows, but the result is stored as a collection of empty files so Snakemake can
	#handle subsequent parallelization of the Pilon jobs.
	input:
		choose_pilon_input(),
		choose_pilon_input()[0] + '.fai'
	output: directory('{sample}/2.polish/pilon/ranges')
	singularity: singularity_image
	shell:
		"""
		mkdir {output}
		bedtools makewindows -w 100000 -g {input[1]} | awk '{{print $1,\":\", $2+ 1, \"-\", $3}}'  | tr -d ' ' |
		xargs -n 1 -I foo touch {output}/foo
		"""

rule pilon_subsetrun:
	#Run Pilon on a small subset of the assembly. Short reads mapping to the window are selected and downsampled to 50x target
	#coverage, greatly reducing the computational cost of polishing the region.
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
	#group: "pilon_subsetruns" #this causes problems with slurm integration
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
			--genome {params.fa} --targets {wildcards.range} \
			--unpaired {params.bam} --output {sample}_{wildcards.range} --outdir {params.subrun_folder} \
			--vcf --nostrays --mindepth 1
		bgzip {params.subrun_folder}/{sample}_{wildcards.range}.vcf
		tabix -fp vcf {params.subrun_folder}/{sample}_{wildcards.range}.vcf.gz
		"""

def aggregate_pilon_subsetruns(wildcards):
	#The individually polished sub-regions of the overall assembly must be collected into a consensus assembly. This rule
	#gathers all the outputs that have been generated for input into the aggregation rule.
	checkpoint_output = checkpoints.pilon_ranges.get(**wildcards).output[0]
	result = expand(rules.pilon_subsetrun.output,
		sample=wildcards.sample,
		range=glob_wildcards(os.path.join(checkpoint_output, '{range,[^.]*}'))[0]
		)
	return(result)

rule pilon_aggregate_vcf:
	#In order to aggregate the changes made during short read polishing, this rule collects the contents of all
	#VCF files generated in the subset runs in preparation to generate a corrected consensus fasta. In order that
	#consensus sequence generation doesn't trip up, this VCF must be carefully deduplicated and sorted before
	#being compressed and indexed.
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
	#Generate short read polished consensus sequence from all changes contained within the aggregated VCF file.
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

rule polish_final:
	#Generate the final output of the polishing phase, whether that's with short reads, long reads, both or neither.
	#Colons and any following characters are omitted from contig names.
	input: choose_polish
	output: "{sample}/2.polish/{sample}_polished.fasta"
	shell:
		"""
		cut -f1 -d ':' {input} > {output}
		"""

#Circularization
#########################################################
checkpoint extract_bigtigs:
	#Only genome-scale contigs are tested for circularity. This rule defines what size that is and extracts those contigs
	#to individual fasta files.
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
	#Extracts reads mapping to the genome-scale contigs to be tested for circularity. These reads are reformatted to fastq and compressed.
	input:
		rules.polish_final.output[0] + '.bam',
		"{sample}/3.circularization/1.candidate_genomes/{tig}.fa",
		rules.polish_final.output[0] + '.bam.bai',
	output:
		"{sample}/3.circularization/2.circularization/spanning_tig_circularization/{tig}/{tig}_terminal_reads.fq.gz"
	singularity: singularity_image
	shell:
		"""
		(samtools idxstats {input[0]} | grep {wildcards.tig} | awk '{{if ($2 > 50000) print $1, ":", $2-50000, "-", $2; else print $1, ":", 1, "-", $2 }}' | tr -d ' ';
		 samtools idxstats {input[0]} | grep {wildcards.tig} | awk '{{if ($2 > 50000) print $1, ":", 1, "-", 50000; else print $1, ":", 1, "-", $2 }}' | tr -d ' ') |
		xargs -I foo sh -c 'samtools view -h {input[0]} foo | samtools fastq - || true' | paste - - - - | sort | uniq | tr '\t' '\n' | bgzip > {output}
		"""

rule circularize_assemble:
	#Reads extracted in circularize_bam2reads are assembled with Canu in order to recover a contig spanning the two
	#ends of a circular genome contig
	input:
		rules.circularize_bam2reads.output
	output: "{sample}/3.circularization/2.circularization/spanning_tig_circularization/{tig}/assembly.fasta"
	params:
		directory="{sample}/3.circularization/2.circularization/spanning_tig_circularization/{tig}",
	singularity: "docker://quay.io/biocontainers/flye:2.4.2--py27he860b03_0"
	resources:
		time=4,
		mem=32
	threads: 4
	shell:
		"""
		flye -t {threads} --nano-raw {input} -o {params.directory} -g 1m
		"""
		#needed for canu:
		#canu -useGrid=False -assemble -p {wildcards.tig} -d {params.directory}  \
		#-nanopore-corrected {input} genomeSize=100000
		#{tig}.contigs
		#singularity_image #

rule circularize_spantig_pre:
	#Prepare to determine if the contig assembled in circularize_assemble actually spans the two ends of the putative genome contig.
	#Preparation entails performing an alignment of the potentially spanning contig to the putative genome contig and post-processing
	#the result.
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
	#Run the script which determines if the putative genome contig is spanned by the smaller contig assembled from terminal reads,
	#indicating circularity.
	input: rules.circularize_spantig_pre.output
	output: "{sample}/3.circularization/2.circularization/spanning_tig_circularization/{tig}/contig_spanned.txt"
	params:
		margin=10000
	script:
		"scripts/spancircle.py"

rule circularize_span_trim:
	#Trim circularized genome contigs at the determined wrap-around point
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
	#Collect the genome sequences produced by spanning contig circularization
	checkpoint_output = checkpoints.extract_bigtigs.get(**wildcards).output[0]
	result = expand(rules.circularize_span_trim.output, #"{sample}/3.circularization/3.circular_sequences/sh/{tig}_span_trim.sh",
		sample=wildcards.sample,
		tig=glob_wildcards(os.path.join(checkpoint_output, '{tig}.fa')).tig)
	return(result)

rule circularize_overcirc:
	#Test putative genome contigs for circularity by self-alignment and determination of 'overcircularization',
	#assembly beyond the circular wraparound point. This produces redundant sequences at the termini of the putative
	#genome contig which can be determined by corner-cutting off-diagonal alignments in a self-alignment dotplot.
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
	#Trim the genome contig circularized by overcircularization detection according to the determined wraparound point.
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
			print('No over-circularization')
		else:
			trim = span_out[0].strip()
			trim_cmd = 'samtools faidx ' + input[1] + ' ' + trim + " > " + params.outfa + "\n"
			cmd += trim_cmd

		open(output[0], 'w').write(cmd + '\n')

def aggregate_overcirc_trim(wildcards):
	#Collect the circular genome sequences obtained by overcircularization detection.
	checkpoint_output = checkpoints.extract_bigtigs.get(**wildcards).output[0]
	result = expand(rules.circularize_overcirc_trim.output, #"{sample}/3.circularization/3.circular_sequences/sh/{tig}_span_overcirc.sh",
		sample=wildcards.sample,
		tig=glob_wildcards(os.path.join(checkpoint_output, '{tig}.fa')).tig)
	return(result)

rule circularize_final:
	#Collect the circular genome sequences and add them back into the total assembly fasta files
	#This is done by generating a .sh file for each putative genome which either yields the circularized
	#sequence or no sequence, depending on whether circularization was successful.
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
		find {sample}/3.circularization/3.circular_sequences/sh/ -type f | xargs cat | bash
		ls {sample}/3.circularization/3.circular_sequences/ | grep .fa$ | cut -f1-2 -d '_' > circs.tmp || true
		(cat {input[1]} | grep -vf circs.tmp |
		cut -f1 | xargs samtools faidx {input[0]}; ls {sample}/3.circularization/3.circular_sequences/* | grep .fa$ | xargs cat) |
		tr -d '\\n' | sed 's/\\(>[contigscaffold]*[_0-9]*\\)/\\n\\1\\n/g' | fold -w 120 | cut -f1 -d ':' | grep -v '^$' > {output} || true
		rm circs.tmp
		"""

#Misassembly detection
#########################################################

def skip_circularization_or_not():
	#Allows all circularization to be skipped when unneeded
	if config['skip_circularization'] == 'True' or config['skip_circularization'] == True:
		return(rules.polish_final.output)
	else:
		return(rules.circularize_final.output)

rule final:
	#Perform one last round of misassembly detection and removal, followed by generation of the final output assembly file
	input: skip_circularization_or_not()[0].replace('.fasta', '.corrected.fasta') #perform one last round of misassembly breakage
	output: "{sample}/5.final/{sample}_final.fa"
	shell: "cp {input} {output}"

#Utility functions
#########################################################

rule align_paf:
	#Align long reads to a fasta, storing the result in .paf format.
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
	#Align long reads to a fasta, storing the result in .bam format.
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
	#index a .bam file
	input:
		'{some}.bam'
	output:
		'{some}.bam.bai'
	singularity: singularity_image
	shell:
		"samtools index {input}"

rule faidx:
	#index a .fasta file
	input: '{something}.f{asta}'
	output: '{something}.f{asta}.fai'
	singularity: singularity_image
	shell: "samtools faidx {input}"
