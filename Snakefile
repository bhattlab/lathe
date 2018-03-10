#!/usr/bin/env python

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
        time=1,
        ntasks=1
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
        time=24,
        ntasks=4
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
        time=48,
        ntasks=24
#    config:
#        "--partition 'nih_s10'"
    shell:
        "bwa mem -t 24 {samp}/idx/{samp} {r1} {r2} ".format(
            samp=config['sample'],
            r1 = config['reads1'],
            r2 = config['reads2']
            ) +
        " | samtools sort --threads 23 > {samp}/{samp}.bam".format(samp=config['sample'])

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
        time=24,
        ntasks=24
    shell:
        "~/moss/tools/classification_and_binning/metabat/runMetaBat.sh --seed 1 -t 24 --unbinned " \
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
        time=24,
        ntasks=24
    shell:
        "module add prodigal; checkm lineage_wf -t 24 -x fa --tab_table -f {samp}/checkm/checkm.tsv {samp}/bins/ {samp}/checkm".format(samp=config['sample'])

rule aragorn:
    input:
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    output:
        "{samp}/rna/trna/bin.1.fa.txt".format(samp=config['sample'])
    log:
        "{samp}/logs/aragorn.log".format(samp=config['sample'])
    resources:
        mem=8,
        time=1,
        ntasks=24
    shell:
        "module add aragorn; \
        ls {samp}/bins/ | xargs -n 1 -I foo -P 24 sh -c 'aragorn -t {samp}/bins/foo -o {samp}/rna/trna/foo.txt';".format(samp=config['sample'])
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
        time=1,
        ntasks=24
    shell:
        "module add barrnap; \
        ls {samp}/bins/ | xargs -n 1 -I foo -P 24 sh -c 'barrnap {samp}/bins/foo > {samp}/rna/rrna/foo.txt'".format(samp=config['sample'])

rule quast:
    input:
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    output:
        "{samp}/quast/".format(samp=config['sample'])
    log:
        "{samp}/logs/quast.log".format(samp=config['sample'])
    resources:
        mem=8,
        time=2,
        ntasks=24
    shell:
        "ls {samp}/bins/ | xargs -n 1 -I foo -P 24 sh -c 'quast.py -o {samp}/quast/foo {samp}/bins/foo \
        --contig-thresholds 0,10000,50000,100000,250000,500000,1000000,2000000,3000000 --fast '".format(samp=config['sample'])

rule prokka:
    input:
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    output:
        "{samp}/prokka/{samp}.gff".format(samp=config['sample'])
    log:
        "{samp}/logs/prokka.log".format(samp=config['sample'])
    resources:
        mem=32,
        time=6,
        ntasks=24
    shell:
        "module add prokka; \
        prokka --outdir {samp}/prokka/ --prefix {samp} --centre X --compliant --force --cpus 24 ".format(samp=config['sample']) +
        "{input}"

rule kraken:
        input:
            "{samp}/idx/{samp}.fa".format(samp=config['sample'])
        output:
            "{samp}/classify/{samp}.krak".format(samp=config['sample'])
        log:
            "{samp}/logs/kraken.log".format(samp=config['sample'])
        resources:
            mem=320,
            time=6,
            ntasks=24
        config:
            "--partition = nih_s10"
        shell:
            "kraken --db {config['krakendb']} --fasta-input {input} --output {output} --preload --threads 24"

rule kraken_translate:
    input:
        "{samp}/classify/{samp}.krak".format(samp=config['sample'])
    output:
        "{samp}/classify/{samp}.tsv".format(samp=config['sample'])
    log:
        "{samp}/logs/kraken.log".format(samp=config['sample'])
    resources:
        mem=8,
        time=1,
        ntasks=1
    shell:
        "kraken-translate {input} --mpa-format --db {config['krakendb']} > {output}"

rule label_bins:
    input:
        "{samp}/classify/{samp}.tsv".format(samp=config['sample']),
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    output:
        "{samp}/classify/bin_species_calls.tsv"
    script:
        "scripts/assign_species.py"


rule postprocess:
    input:
        "{samp}/prokka/{samp}.gff".format(samp=config['sample']),
        "{samp}/quast/".format(samp=config['sample']),
        "{samp}/checkm/checkm.tsv".format(samp=config['sample']),
        "{samp}/rna/rrna/bin.1.fa.txt".format(samp=config['sample']), #rrna.tsv
        "{samp}/rna/trna/bin.1.fa.txt".format(samp=config['sample'])
    output:
        "{samp}/final/{samp}.tsv".format(samp=config['sample'])
    resources:
        mem=8,
        time=1,
        ntasks=1
    shell:
        "echo 'do final data collation.'"


    #collate assembly stats per organism for figure
#	(echo Gut+Library | paste - <(find quast_results/ | grep transposed_report.tsv | head -1 | xargs cat | grep Assembly) | tr '+' '\t'; find quast_results/ | grep transposed_report.tsv | xargs cat | grep -v Assembly | sort -k1,1 | sed 's/\.asm_/.asm\t/g' | sed 's/_canu_/_canu\t/g' | sed 's/spades_/spades\t/g' | sed 's/^ab./ab\t/g' | sed 's/^em./em\t/g') > assembly_stats.tsv

    #count genes per bin with prokka
    #module add prokka
#    ls species_combined_links | xargs -n 1 -I foo ssubN_Xhrs prokka_foo 40 100 prokka --outdir ../outs/foo ../../bin_links/foo --centre X --compliant --force
#    find | grep gff$ | xargs grep -c CDS | tr ':' '\t' | sed 's/^\.\///g' | tr '/' '_' | sed 's/_PROKKA_.*gff//g' |  sed 's/outs_//g' | sed 's/\.fa//g' | sed 's/\.asm_/.asm\t/g' | sed 's/_canu_/_canu\t/g' | sed 's/spades_/spades\t/g' | sed 's/^ab./ab\t/g' | sed 's/^em./em\t/g' > gene_counts.tsv

 	#count genes locally
#    nohup sh -c "ls ../../species_combined_links | grep -v .fa.fai | xargs -n 1 -I foo -P 32 prokka --centre XXX --compliant --force --outdir ../outs/foo ../../species_combined_links/foo " &
