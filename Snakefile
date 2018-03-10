#!/usr/bin/env python

localrules: bwa_index_setup, postprocess_raw, postprocess_final

rule targets:
    input:
        config['reads1'],
        config['reads2'],
        config['assembly']

rule bwa_index_setup:
    input:
        config['assembly']
    output:
        "{samp}/idx/{samp}.fa".format(samp = config['sample'])
    resources:
        mem=1,
        time=1
    threads: 1
    shell:
        "cp {asm} {samp}/idx/{samp}.fa".format(asm = config['assembly'], samp = config['sample'])


rule bwa_index:
    input:
        "{samp}/idx/{samp}.fa".format(samp=config['sample'])
    output:
        "{samp}/idx/{samp}.amb".format(samp=config['sample']),
        "{samp}/idx/{samp}.ann".format(samp=config['sample']),
        "{samp}/idx/{samp}.bwt".format(samp=config['sample']),
        "{samp}/idx/{samp}.pac".format(samp=config['sample']),
        "{samp}/idx/{samp}.sa".format(samp=config['sample'])
    log:
        "{samp}/logs/bwa_index.log".format(samp=config['sample'])
    resources:
        mem=8,
        time=24
    threads: 1
    params:
        prefix="{samp}/idx/{samp}".format(samp=config['sample']),
        algorithm="bwtsw"
    wrapper:
        "0.22.0/bio/bwa/index"


rule bwa_align:
    input:
        "{samp}/idx/{samp}.fa".format(samp=config['sample']),
        config['reads1'],
        #config['reads2'],
        "{samp}/idx/{samp}.amb".format(samp=config['sample']),
        "{samp}/idx/{samp}.ann".format(samp=config['sample']),
        "{samp}/idx/{samp}.bwt".format(samp=config['sample']),
        "{samp}/idx/{samp}.pac".format(samp=config['sample']),
        "{samp}/idx/{samp}.sa".format(samp=config['sample'])
    log:
        "{samp}/logs/bwa_mem.log".format(samp=config['sample'])
    output:
        "{samp}/{samp}.bam".format(samp=config['sample'])
    resources:
        mem=32,
        time=48
    threads: 24
    shell:
        "bwa mem -t {threads}" + " {samp}/idx/{samp} {r1} {r2} ".format(
            samp=config['sample'],
            r1 = config['reads1'],
            r2 = config['reads2']
            ) +
        " | samtools sort --threads {threads}" + " > {samp}/{samp}.bam".format(samp=config['sample'])

rule metabat:
    input:
        "{samp}/idx/{samp}.fa".format(samp=config['sample']),
        "{samp}/{samp}.bam".format(samp=config['sample'])
    output:
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    log:
        "{samp}/logs/metabat.log".format(samp=config['sample'])
    resources:
        mem=64,
        time=24
    threads: 24
    shell:
        "~/moss/tools/classification_and_binning/metabat/runMetaBat.sh --seed 1 -t {threads} --unbinned " \
        + "{input};" +
        " mv {samp}.fa.metabat-bins--unbinned/* {samp}/bins/; \
        rmdir {samp}.fa.metabat-bins--unbinned/; \
        mv {samp}.fa.depth.txt {samp}/".format(samp=config['sample'])

rule checkm:
    input:
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    output:
        "{samp}/checkm/checkm.tsv".format(samp=config['sample'])
    log:
        "{samp}/logs/checkm.log".format(samp=config['sample'])
    resources:
        mem=64,
        time=24
    threads: 24
    shell:
        "module add prodigal; checkm lineage_wf -t {threads}" +" -x fa --tab_table -f {samp}/checkm/checkm.tsv {samp}/bins/ {samp}/checkm".format(samp=config['sample'])

rule aragorn:
    input:
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    output:
        "{samp}/rna/trna/bin.1.fa.txt".format(samp=config['sample'])
    log:
        "{samp}/logs/aragorn.log".format(samp=config['sample'])
    resources:
        mem=8,
        time=1
    threads: 24
    shell:
        "module add aragorn; \
        ls {samp}/bins/ | grep -v '.fa.fai$'| xargs -n 1 -I foo ".format(samp=config['sample']) +
        "-P {threads} " +
        "sh -c 'aragorn -t {samp}/bins/foo -o {samp}/rna/trna/foo.txt';".format(samp=config['sample'])
#        grep Total {samp}/rna/trna/*.txt | sed 's/\.fa_trna\.txt\:Total tRNA genes = /\t/g' | sed 's/trna\///g' |sed 's/\.asm_/.asm\t/g' | sed 's/_canu_/_canu\t/g' | sed 's/spades_/spades\t/g' | sed 's/^ab./ab\t/g' | sed 's/^em./em\t/g' > trna.tsv

rule barrnap:
    input:
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    output:
        "{samp}/rna/rrna/bin.1.fa.txt".format(samp=config['sample']) #rrna.tsv
    log:
        "{samp}/logs/barrnap.log".format(samp=config['sample'])
    resources:
        mem=8,
        time=1
    threads: 24
    shell:
        "module add barrnap; \
        ls {samp}/bins/ | grep -v '.fa.fai$'| xargs -n 1 -I foo ".format(samp=config['sample']) +
        "-P {threads} " + "sh -c 'barrnap {samp}/bins/foo > {samp}/rna/rrna/foo.txt'".format(samp=config['sample'])

rule quast:
    input:
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    output:
        "{samp}/quast/".format(samp=config['sample'])
    log:
        "{samp}/logs/quast.log".format(samp=config['sample'])
    resources:
        mem=8,
        time=2
    threads: 24
    shell:
        "ls {samp}/bins/ | grep -v '.fa.fai$' | xargs -n 1 -I foo -P ".format(samp=config['sample']) + \
		"{threads} " + \
		"sh -c 'quast.py -o {samp}/quast/foo {samp}/bins/foo \
        --contig-thresholds 0,10000,50000,100000,250000,500000,1000000,2000000,3000000 --fast '".format(samp=config['sample'])

rule prokka:
    input:
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    output:
        "{samp}/prokka/bin.1.fa/{samp}_bin.1.fa.gff".format(samp = config['sample'])
    log:
        "{samp}/logs/prokka.log".format(samp=config['sample'])
    resources:
        mem=32,
        time=6
    threads: 48
    shell:
        "module add prokka; \
		ls {samp}/bins/ | grep -v '.fa.fai$' | xargs -n 1 -I foo -P ".format(samp=config['sample']) + \
		"{threads} " + \
        "prokka {samp}/bins/foo --outdir {samp}/prokka/foo --prefix {samp}_foo --centre X --compliant --force ".format(samp=config['sample'])

rule kraken:
        input:
            "{samp}/idx/{samp}.fa".format(samp=config['sample'])
        output:
            "{samp}/classify/{samp}.krak".format(samp=config['sample'])
        log:
            "{samp}/logs/kraken.log".format(samp=config['sample'])
        resources:
            mem=320,
            time=6
        threads: 24
        shell:
            "kraken --db {krak} ".format(krak = config['krakendb']) +
            " --fasta-input {input} --output {output} --preload --threads {threads}"

rule kraken_translate:
    input:
        "{samp}/classify/{samp}.krak".format(samp=config['sample'])
    output:
        "{samp}/classify/{samp}.tsv".format(samp=config['sample'])
    log:
        "{samp}/logs/kraken.log".format(samp=config['sample'])
    resources:
        mem=8,
        time=1
    shell:
        "kraken-translate {input} --mpa-format " +
        "--db {krak} ".format(krak=config['krakendb']) +
        "> {output}"

rule fasta_index:
    input:
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    output:
        "{samp}/bins/bin.1.fa.fai".format(samp=config['sample'])
    log:
        "{samp}/logs/faidx.log".format(samp=config['sample'])
    resources:
        mem=8,
        time=1
    threads: 24
    shell:
        "ls {samp}/bins/*.fa | xargs -n 1 -P ".format(samp=config['sample']) +
        " {threads} samtools faidx "

rule label_bins:
    input:
        "{samp}/classify/{samp}.tsv".format(samp=config['sample']),
        "{samp}/bins/bin.1.fa.fai".format(samp=config['sample'])
    output:
        "{samp}/classify/bin_species_calls.tsv".format(samp=config['sample'])
    log:
        "{samp}/logs/assign_species.log".format(samp=config['sample'])
    script:
        "assign_species.py"


rule postprocess_raw:
	input:
		rules.prokka.output,
		rules.quast.output,
		rules.checkm.output,
		rules.aragorn.output,
		rules.barrnap.output,
		rules.label_bins.output
	output:
		temp("{samp}/final/prokka.tmp".format(samp=config['sample'])),
		temp("{samp}/final/quast.tmp".format(samp=config['sample'])),
		temp("{samp}/final/checkm.tmp".format(samp=config['sample'])),
		temp("{samp}/final/trna.tmp".format(samp=config['sample'])),
		temp("{samp}/final/rrna.tmp".format(samp=config['sample'])),
		temp("{samp}/final/classify.tmp".format(samp=config['sample']))
	log:
		"{samp}/logs/postprocess.log".format(samp=config['sample'])
	resources:
		mem=2,
		time=1
	shell:
		"(echo 'Sample Bin Genes' | tr ' ' '\t'; \
			find {samp}/prokka/ -name '*.gff' | xargs grep -c CDS | cut -f4 -d '/' | sed 's/.fa.gff:/\t/g' | sed 's/{samp}_/{samp}\t/g' | sort -k2,2g) > {samp}/final/prokka.tmp; ".format(samp=config['sample']) + \
		"find {samp}/quast/ -name 'transposed_report.tsv' | xargs cat | sort -u | sed 's/^bin/{samp}\tbin/g' | sed 's/^Assembly/Sample\tBin/g' > {samp}/final/quast.tmp; ".format(samp=config['sample']) + \
		"cat {rules.checkm.output} " + " | sed 's/^bin/{samp}\tbin/g' | sed 's/^Bin Id/Sample\tBin/g' > {samp}/final/checkm.tmp; ".format(samp=config['sample']) + \
		"(echo 'Sample Bin tRNA' | tr ' ' '\t'; \
			grep Total {samp}/rna/trna/* | sed 's/\/rna\/trna\//\t/g' | sed 's/.fa.txt:Total tRNA genes = /\t/g') > {samp}/final/trna.tmp; ".format(samp=config['sample']) + \
		"(echo 'Sample Bin rna.16S rna.23S rna.5S' | tr ' ' '\t'; \
			paste \
			<(grep -c 16S 6753_d/rna/rrna/* | sed 's/\/rna\/rrna\//\t/g' | sed 's/.fa.txt:/\t/g') \
			<(grep -c 23S 6753_d/rna/rrna/* | sed 's/\/rna\/rrna\//\t/g' | sed 's/.fa.txt:/\t/g' | cut -f3) \
			<(grep -c 5S 6753_d/rna/rrna/* | sed 's/\/rna\/rrna\//\t/g' | sed 's/.fa.txt:/\t/g' | cut -f3) \
		) > {samp}/final/rrna.tmp; ".format(samp=config['sample']) + \
		"(echo 'Sample Bin Majority.Class Majority.Fraction Size.Mb Final.Class' | tr ' ' '\t'; \
			cat {rules.label_bins.output}" + " | sed 's/\.fa//g' | sed 's/^/{samp}\t/g') > {samp}/final/classify.tmp;".format(samp=config['sample'])

rule postprocess_final:
	input:
		rules.postprocess_raw.output
	output:
		"{samp}/final/{samp}.tsv".format(samp=config['sample'])
	log:
		"{samp}/logs/postprocess.log".format(samp=config['sample'])
	script:
		"join_final_tables.R"
