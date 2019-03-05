
'''
rule circlator:
	input:
		rules.unique_tignames.output,
		'{{sample}}/1.assemble_{g}/{{sample}}_{g}.correctedReads.fasta.gz'.format(g = config['genome_size'].split(",")[0]),
	output: '{sample}/2.circlator/06.fixstart.fasta'
	threads: 16
	resources:
		time=100,
		mem=200
	singularity: singularity_image
	shell:
		"""
		conda activate circlator
		rmdir {sample}/2.circlator/
		circlator all --verbose --threads {threads} \
		--merge_min_id 85 --merge_breaklen 4000 --merge_reassemble_end 50000 --assembler canu \
		--data_type nanopore-corrected --bwa_opts '-x ont2d'  \
		{input} {sample}/2.circlator
		"""
#--split_all_reads split all reads may be used with canu assembly #"--merge_min_length_merge 500 "
'''


'''
rule pilon:
	input:
		rules.circlator_fixstart.output,
		rules.align_short_reads.output,
		rules.align_short_reads.output[0] + '.bai'
	output: "{sample}/3a.pilon/{sample}_pilon.fasta"
	resources:
		mem=8000,
		time=48
	params:
		java_mem = 8000
	shell:
		"java -Xmx{params.java_mem}G -jar $(which pilon | sed 's/\/pilon//g')/../share/pilon*/pilon*.jar --genome {input[0]} " +
		"--frags {input[1]} --output {sample}_pilon --outdir {sample}/3a.pilon/ --tracks" #

'''

'''
rule extract_candidate_genome_tigs:
	input:
		rules.circlize.output,
		rules.circlize.output[0] + '.fai'
	output:
		dynamic("{sample}/2b.circlator_big_tigs/{tigname}.tigname")
	shell:
		#"sort -k2,2gr {input[1]} | awk '{{if ($2 > 2000000) print $1}}' > {output}"
		"sort -k2,2gr {input[1]} | awk '{{if ($2 > 2000000) print $1}}' | " +
		"xargs -n 1 -I foo sh -c 'echo foo > {sample}/2b.circlator_big_tigs/foo.tigname'"

rule circlator_progcheck:
	input:
		rules.merge.output
	output: '{sample}/2.circlator/0.progcheck_ok'
	shell: "circlator progcheck && touch {output}"

rule circlator_mapreads:
	input:
		rules.merge.output,
		'{{sample}}/1.assemble_{g}/{{sample}}_{g}.correctedReads.fasta.gz'.format(g = config['genome_size'].split(",")[0]),
		rules.circlator_progcheck.output
	output:	"{sample}/2.circlator/1.mapreads.bam"
	threads: 24
	resources:
		time=24,
		mem=48
	shell: "circlator mapreads {input[0]} {input[1]} {output} --threads {threads} --bwa_opts '-x ont2d' --verbose"

rule circlator_bam2reads:
	input: rules.circlator_mapreads.output
	output: "{sample}/2.circlator/2.bam2reads.fasta"
	resources:
		mem=8,
		time=6
	shell: "circlator bam2reads --discard_unmapped {input} {wildcards.sample}/2.circlator/2.bam2reads"

rule circlator_assemble:
	input: rules.circlator_bam2reads.output
	output: "{sample}/2.circlator/3.assemble/canu.contigs.fasta"
	threads: 24
	resources:
		time=24,
		mem=180
	shell:
		"circlator assemble --threads {threads} --assembler canu --data_type nanopore-corrected {input} {wildcards.sample}/2.circlator/3.assemble/"

rule circlator_merge:
	input:
		rules.merge.output,
		'{{sample}}/1.assemble_{g}/{{sample}}_{g}.correctedReads.fasta.gz'.format(g = config['genome_size'].split(",")[0]),
		rules.circlator_assemble.output
	output:
		"{sample}/2.circlator/4.merge.fasta"
	threads: 4
	resources:
		time=80,
		mem=180
	shell: "circlator merge {input[0]} {input[2]} {wildcards.sample}/2.circlator/4.merge " +
			"--min_id 85 --breaklen 4000 --reassemble_end 50000 --threads {threads} --reads {input[1]} " +
			"--assembler canu --data_type nanopore-corrected"

rule circlator_clean:
	input: rules.circlator_merge.output
	output: "{sample}/2.circlator/5.clean.fasta"
	resources:
		time=24,
		mem=20
	shell: "circlator clean {input} {wildcards.sample}/2.circlator/5.clean"

rule circlator_fixstart:
	input: rules.circlator_clean.output
	output: "{sample}/2.circlator/6.fixstart.fasta"
	resources:
		mem=24,
		time=24
	shell: "circlator fixstart {input} {wildcards.sample}/2.circlator/6.fixstart"

rule circlize_2mb:
	input:
		"{sample}/2.circlator/6.fixstart.fasta"
		#rules.extract_candidate_genome_tigs.output
		"{sample}/2b.circlator_big_tigs/{tigname}.tigname"
	output: '{sample}/2b.circlator_big_tigs/{tigname}/6.fixstart.fasta'
	threads: 16
	resources:
		time=24,
		mem=60
	shell:
		'rmdir {sample}/2b.circlator_big_tigs/{wildcards.tigname}; ' +
		'circlator all --verbose --threads {threads} ' +
		' --merge_min_id 85 --merge_min_length 1600 --merge_breaklen 4000 --assembler canu ' + #--split_all_reads split all reads may be used with canu assembly
		"--data_type nanopore-corrected --bwa_opts '-x ont2d' --merge_reassemble_end 50000 " + #"--merge_min_length_merge 500 " +
		'{input[0]} --b2r_discard_unmapped --b2r_only_contigs {input[1]} ' +
		'{sample}/1.assemble_{g}/{sample}.correctedReads.fasta.gz {sample}/2b.circlator_big_tigs/{wildcards.tigname}/'

rule circlize_final:
	input:
		dynamic(rules.circlize_2mb.output)
	output:
		'{sample}/2c.circlator_final/{sample}_circularized.fa',
		'{sample}/2c.circlator_final/circularized_contigs.tsv'
	shell:
		"cut -f1 {sample}/2.circlator/sanitized.fa.fai > {sample}/2.circlator/sanitized_tigs.list; " + #aggregate sequences non-redundantly
		"ls {sample}/2b.circlator_big_tigs | grep -v txt > {sample}/2.circlator/recirclized_tigs.list; " +
		"(grep -vf {sample}/2.circlator/recirclized_tigs.list {sample}/2.circlator/sanitized_tigs.list | " +
		"xargs samtools faidx {sample}/2.circlator/sanitized.fa; " +
		"cat {input}) > {output[0]}; " +

		#aggregate circularization logs
		" (grep -P '\\t1' {sample}/2.circlator/04.merge.circularise.log; " +
		"grep -P '\\t1' {sample}/2b.circlator_big_tigs/*/04.merge.circularise.log | cut -f2,99 -d ':') " +
		"| sed -e 's/:[^tig\\t]*/_/g' | sed 's/_\\t/\\t/g' | sort -u > {output[1]}"
'''


rule tombo_annotate:
	input:
		#lambda wildcards: "{s}/0.basecall/raw_calls/{f1}/{f2}/sequencing_summary.txt".format(s = sample, f1 = wildcards.dir.split('___')[0], f2 = wildcards.dir.split('___')[1])
		#this horrible input specifies the precise basecalling run needed for this batch of reads.  It is simpler to require that all basecalling be complete, but less correct.
		rules.basecall_final.output
	output: '{sample}/4.tombo/0.annotate/{dir}_ANNOTATED'
	threads: 1
	resources:
		time=lambda wildcards, attempt: 6 * (2**(attempt - 1)),
		mem=8
	params:
		dir_to_annotate = lambda wildcards: fast5_abspath_run_subfolder[wildcards.dir.replace('___', '/')],
		f1 = lambda wildcards: wildcards.dir.split('___')[0],
		f2 = lambda wildcards: wildcards.dir.split('___')[1]
	singularity: singularity_image
	shell:
		"source activate tombo && " +
		"tombo preprocess annotate_raw_with_fastqs --overwrite --fast5-basedir {params.dir_to_annotate} " +
		"--fastq-filenames $(ls {sample}/0.basecall/raw_calls/{params.f1}/{params.f2}/workspace/pass/*.fastq) " +
	 	"--sequencing-summary-filenames {sample}/0.basecall/raw_calls/{params.f1}/{params.f2}/sequencing_summary.txt " +
	 	"--processes {threads} && touch {output}"

rule tombo_resquiggle:
	input:
		rules.final_assembly.output,
		'{sample}/4.tombo/0.annotate/{dir}_ANNOTATED'
	output:
		'{sample}/4.tombo/1.resquiggle/{dir}_RESQUIGGLED'
	threads: 1
	resources:
		time=lambda wildcards, attempt: 12 * (2**(attempt - 1)),
		mem=8
	params:
		dir_to_squiggle = lambda wildcards: fast5_abspath_run_subfolder[wildcards.dir.replace('___', '/')],
		index_filename = lambda wildcards: '../.' + wildcards.dir.split('___')[1] + '.RawGenomeCorrected_000.tombo.index'
	singularity: singularity_image
	shell:
		"source activate tombo; " +
		"tombo resquiggle --overwrite {params.dir_to_squiggle} " +
		"{input[0]} --processes {threads} --num-most-common-errors 5 && " +
		"ln {params.dir_to_squiggle}/{params.index_filename} {output}"

rule tombo_detect:
	input: expand('{{sample}}/4.tombo/1.resquiggle/{dir}_RESQUIGGLED', dir = [d.replace('/', '___') for d in fast5_abspath_run_subfolder])
	output:
		#'{sample}/4.tombo/readstats.{methyl}.tombo.per_read_stats',
		protected('{sample}/4.tombo/stats.{methyl}.tombo.stats')
	threads: 48
	resources:
		mem=lambda wildcards, attempt: 100 * attempt,
		time=lambda wildcards, attempt: 48 * attempt
	params:
		fast5_dirs = [d for d in fast5_run_subfolder_abspath]
	singularity: singularity_image
	shell:
		"source activate tombo && " +
		"tombo detect_modifications alternative_model --fast5-basedirs {params.fast5_dirs} " +
	 	"--statistics-file-basename {sample}/4.tombo/stats " +
	 	#"--per-read-statistics-basename {sample}/4.tombo/readstats " +
	 	"--alternate-bases {wildcards.methyl} --processes {threads} --dna "

rule tombo_wig:
	input:
		rules.tombo_detect.output,
	output:
		"{sample}/4.tombo/{sample}_{methyl}_sites.fraction_modified_reads.plus.wig"
	singularity: singularity_image
	shell:
		"source activate tombo; " +
		"tombo text_output browser_files --browser-file-basename {sample}/4.tombo/{sample}_{wildcards.methyl}_sites " +
		"--file-type fraction --statistics-filename {input} "

rule tombo_bed_fasta:
	input:
		rules.tombo_wig.output,
		rules.final_assembly.output,
	output:
		"{sample}/4.tombo/{sample}_{methyl}_sites.bed",
		"{sample}/4.tombo/{sample}_{methyl}_sites.fa"
	params:
		meth_frac = 0.5,
		site_radius = 2
	singularity: singularity_image
	shell:
		"cat {input[0]} | wig2bed | " +
		"awk '{{if ($5 > {params.meth_frac}) print $1, $2-{params.site_radius}, $3+{params.site_radius}, $4, $5}}' | " +
		"tr ' ' '\\t' | awk '{{if ($2 < 0) print $1, 0, $3, $4, $5; else print $0}}' | tr ' ' '\\t' > {output[0]}; " +
		"bedtools getfasta -s -fi {input[1]} -bed {output[0]} | tr -d '()' > {output[1]}"

rule plasmid_kmer_methyl:
	input:
		rules.tombo_bed_fasta.output,
		rules.final_assembly.output[0] + '.fai'
	output: dynamic("{sample}/tmp_kmer_methyl_{methyl}/jellyfish/tsv/{tig}.tsv")
	params:
		kmer_length = 5
	singularity: singularity_image
	shell:
		"mkdir -p {sample}/tmp_kmer_methyl_{wildcards.methyl}/fa {sample}/tmp_kmer_methyl_{wildcards.methyl}/jellyfish/db {sample}/tmp_kmer_methyl_{wildcards.methyl}/jellyfish/tsv;\n" +
		"cut -f1 {input[2]} | xargs -n 1 -I foo sh -c \"grep -A1 foo {input[1]} > {sample}/tmp_kmer_methyl_{wildcards.methyl}/fa/foo_sites.fa\";\n" +
		"cut -f1 {input[2]} | xargs -n 1 -I foo sh -c \"  \
jellyfish count -C {sample}/tmp_kmer_methyl_{wildcards.methyl}/fa/foo_sites.fa -m {params.kmer_length} -s 10000 -o {sample}/tmp_kmer_methyl_{wildcards.methyl}/jellyfish/db/foo.jf; \
jellyfish dump {sample}/tmp_kmer_methyl_{wildcards.methyl}/jellyfish/db/foo.jf --tab --column | sort -k1,1 > {sample}/tmp_kmer_methyl_{wildcards.methyl}/jellyfish/tsv/foo.tsv\";\n "

rule plasmid_kmer_methyl_join:
	input: dynamic("{sample}/tmp_kmer_methyl_{{methyl}}/jellyfish/tsv/{{tig}}.tsv".format(sample = sample))
	output: "{sample}/5.assign_plasmids/kmer/methyl_kmers_{methyl}.tsv"
	run:
		"join_several {input} > {output}"
	#script:
	#	"scripts/join_several2.py"

rule plasmid_kmer_tig:
	input:
		rules.final_assembly.output[0],
		rules.final_assembly.output[0] + '.fai'
	output: dynamic("{sample}/tmp_kmer_wholetig/jellyfish/whole_tsv/{tig}.tsv")
	params:
		kmer_length = 5
	singularity: singularity_image
	shell:
		"mkdir -p {sample}/tmp_kmer_wholetig/jellyfish/whole_tsv {sample}/tmp_kmer_wholetig/jellyfish/db {sample}/tmp_kmer_wholetig/fa; \n" +
		"cat {input[1]} | cut -f1 | xargs -n1 -I foo -P 8 sh -c 'samtools faidx {input[0]} foo > {sample}/tmp_kmer_wholetig/fa/foo.fa; \n\
jellyfish count -C {sample}/tmp_kmer_wholetig/fa/foo.fa -m {params.kmer_length} -s 10000 -o {sample}/tmp_kmer_wholetig/jellyfish/db/foo.jf; \n\
jellyfish dump {sample}/tmp_kmer_wholetig/jellyfish/db/foo.jf --tab --column | sort -k1,1 > {sample}/tmp_kmer_wholetig/jellyfish/whole_tsv/foo.tsv'; \n"

rule plasmid_kmer_tig_join:
	input: dynamic("{sample}/tmp_kmer_wholetig/jellyfish/whole_tsv/{{tig}}.tsv".format(sample = sample))
	output:"{sample}/5.assign_plasmids/kmer/tig_kmers.tsv",
	run:
		"join_several {input} > {output}"
	#script:
	# 	"scripts/join_several2.py"

rule plasmid_dendro_prebin:
	input:
		"{sample}/5.assign_plasmids/kmer/methyl_kmers_5mC.tsv",
		"{sample}/5.assign_plasmids/kmer/methyl_kmers_6mA.tsv",
		rules.plasmid_kmer_tig_join.output
	output:
		"{sample}/5.assign_plasmids/pre_binning_steps_done"
	shell:
		"touch {output}"

rule plasmid_list:
	input:
		rules.final_assembly.output,
		rules.final_assembly.output[0] + '.fai'
	output:
		"{sample}/5.assign_plasmids/plasmids.tsv"
	shell:
		"cut -f2 {input[1]} | xargs -n 1 -I foo grep foo {input[2]} | awk '{{if ($2 < 1000000) print $0}}' > {output}"

rule plasmid_dendro:
	input:
		"{sample}/5.assign_plasmids/kmer/methyl_kmers_5mC.tsv",
		"{sample}/5.assign_plasmids/kmer/methyl_kmers_6mA.tsv",
		"{sample}/5.assign_plasmids/kmer/tig_kmers.tsv",
		"{sample}/5.assign_plasmids/bin_tig_mapping.tsv",
		"{sample}/5.assign_plasmids/plasmids.tsv",
		"{sample}/5.assign_plasmids/species_assignments.tsv",
	output:
		"{sample}/5.assign_plasmids/dendro.pdf"
	script:
	 	'scripts/plasmid_dendro.R'


#could be integrated into binning workflow:
#ls *.fai | xargs -n 1 -I foo sh -c "cat foo | sed 's/^/foo\t/g'" | sed 's/.fa.fai//g' | cut -f1,2 > ../bin_tig_mapping.tsv


rule nanopolish_index_fofn:
	input:
		expand('{{sample}}/0.basecall/raw_calls/{foo}/sequencing_summary.txt', foo = fast5_abspath_run_subfolder.keys())
	output:
		"{sample}/4.nanopolish/seq_summ_fofn.list"
	run:
		f = open(output[0], 'w')
		f.write("\n".join(input))
		f.close()

rule nanopolish_index:
	input:
		rules.nanopolish_index_fofn.output,
		config['fast5_parent_dir'],
		rules.basecall_final.output
	output:
		rules.basecall_final.output[0] + '.index'
	resources:
		time=48
	shell:
		"nanopolish index --verbose -f {input[0]} -d {input[1]} {input[2]}" #.format(blar = " -d ".join(fast5_dirs))

rule nanopolish_ranges:
	input: rules.final_assembly.output
	output: dynamic('{sample}/4.nanopolish/ranges/{range}')
	shell:
		"nanopolish_makerange.py {input} | xargs -n 1 -I foo touch {sample}/4.nanopolish/ranges/foo"

rule nanopolish:
	input:
		'{sample}/4.nanopolish/ranges/{range}',
		rules.basecall_final.output,
		rules.final_assembly.output[0] + '.bam',
		rules.final_assembly.output,
		rules.final_assembly.output[0] + '.bam.bai',
		rules.nanopolish_index.output,
	output:
		"{sample}/4.nanopolish/range_vcfs/{sample}_{range}.vcf"
	threads: 4
	resources:
		time= lambda wildcards, attempt: 48 * attempt, #this may be reduced if it incurs scheduling delays
		mem=16
	shell:
		"nanopolish variants --fix-homopolymers " + #--max-haplotypes=10000 " +
		"--consensus -o {output} -w {wildcards.range} -r {input[1]} " +
		"-b {input[2]} -g {input[3]} -t {threads} --min-candidate-frequency 0.1"

rule nanopolish_final:
	input:
		rules.final_assembly.output,
		dynamic(rules.nanopolish.output)
	output:
		"{sample}/4.nanopolish/{sample}_nanopolish.fa"
	shell:
		"nanopolish vcf2fasta -g {input[0]} {sample}/4.nanopolish/range_vcfs/* > {output}"

rule nanopolish_call_methylation:
	input:
		rules.basecall_final.output,
		'{sample}/3.circlize/06.fixstart.fasta.bam',
		rules.final_assembly.output,
		rules.nanopolish_index.output
	output:
		'{sample}/4.nanopolish/{sample}-{methylation_type}.tsv'
	threads: 24
	resources:
		time=24,
		mem=24
	shell:
		"nanopolish call-methylation --reads {input[0]} --bam {input[1]} --genome {input[2]}" +
		"--threads {threads} --methylation {wildcards.methylation_type} > {output}"
