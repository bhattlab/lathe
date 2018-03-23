#!/usr/bin/env python

read_length = int(snakemake.params['read_length'])
idxstats = snakemake.input[0]
output = snakemake.output[0]
bin = snakemake.wildcards['bin']
samp = snakemake.wildcards['samp']

with open(idxstats, 'r') as statsf:
	stats = [l.strip().split("\t") for l in statsf.readlines()]

total_contig_bases = sum([int(i[1]) for i in stats])
total_reads = sum([int(i[2]) for i in stats])
total_read_bases = total_reads * read_length

avg_cov = total_read_bases / total_contig_bases

with open(output, 'w') as outf:
	outf.write("\t".join([samp, bin, str(avg_cov)]) + "\n")
