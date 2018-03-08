#!/usr/bin/env python

rule targets:
    input:
        config['reads1'],
        config['reads2'],
        config['reads_se'],
        config['assembly']

rule bwa_index_setup:
    input:
        config['assembly']
    output:
        "{samp}/idx/{samp}.fa".format(samp = config['sample'])
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
    params:
        prefix="{samp}/idx/{samp}".format(samp=config['sample']),
        algorithm="bwtsw"
    wrapper:
        "0.22.0/bio/bwa/index"


rule bwa_align:
    input:
        "{samp}/idx/{samp}.fa".format(samp=config['sample']),
        config['reads1'],
        config['reads2'],
        "{samp}/idx/{samp}.amb".format(samp=config['sample']),
        "{samp}/idx/{samp}.ann".format(samp=config['sample']),
        "{samp}/idx/{samp}.bwt".format(samp=config['sample']),
        "{samp}/idx/{samp}.pac".format(samp=config['sample']),
        "{samp}/idx/{samp}.sa".format(samp=config['sample'])
    log:
        "{samp}/logs/bwa_mem.log".format(samp=config['sample'])
    output:
        "{samp}/{samp}.bam".format(samp=config['sample']),
        "{samp}/{samp}.bam.bai".format(samp=config['sample'])
    shell:
        "bwa mem -t 64 {samp}/idx/{samp} {r1} {r2} ".format(
            samp=config['sample'],
            r1 = config['reads1'],
            r2 = config['reads2']) +
        " | samtools sort --threads 63 > {samp}/{samp}.bam".format(samp=config['sample'])

rule run_metabat:
    input:
        "{samp}/idx/{samp}.fa".format(samp=config['sample']),
        "{samp}/{samp}.bam".format(samp=config['sample'])
    output:
        "{samp}/bins/bin.1.fa".format(samp=config['sample'])
    #config:
    #    "-m 4"
    #    "-n metabat"
    #    "-c 8"
    log:
        "{samp}/logs/metabat.log".format(samp=config['sample'])
    shell:
        "~/moss/tools/classification_and_binning/metabat/runMetaBat.sh -t 64 {input}" + " -o {samp}/bins/".format(samp=config['sample'])
#    wrapper:
#        "0.2.0/bio/samtools/sort"

'''
rule plot:
    input:
        "raw/{dataset}.csv"
    output:
        "plots/{dataset}.pdf"
    config: #for ssub
    	-n jobname
    shell:
        "somecommand {input} {output}"
'''
