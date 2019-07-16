#!/usr/bin/env python

import sys
import os
import argparse
import pysam

fastaf = snakemake.input[0]
out = open(snakemake.output[0], 'w')

#knobs
smooth_gap_width = 150000
contig_edge_margin = 150000
min_smoothed_aln_len = 10000
min_aln_len = 5000

#align
print("Aligning...".format(t=str(snakemake.threads)))
os.system("nucmer -p {delta} -b 4000 -l 2000 --maxmatch {fa} {fa}".format( #-t {threads} back to mummer3, no more multithreading : (
	threads = str(snakemake.threads),
	delta = snakemake.params['delta'],
	fa = fastaf
))
os.system("show-coords -T {delta}.delta -L 2000 | sed '1,5d' > coords.tsv".format(
	delta = snakemake.params['delta']
)) #the 1,5d gets rid of the header and identity hit as well
print("Trimming genome...")

max_tiglen = 0

with open('coords.tsv') as coords:
	lines = coords.readlines()
	smoothed_lines = []
	if len(lines) > 0:
		aln_start_line = lines[0].strip().split("\t")
		prev_line = lines[0].strip().split("\t")

		#goal: identify corner-cutting parallel off-diagonals.
	#	print(prev_line)

		for l in lines[1:] + [lines[0]]: #just so we don't miss the last item
			s = l.strip().split("\t")
			#print(s)

			tigname = s[-1] #store the name of the tig.  this repeats nugatorily.

			if int(s[1]) > max_tiglen: #keep track of the largest alignment coordinate found, used as the tig length
				max_tiglen = int(s[1])

			if int(s[0]) > int(s[1]): #ignore inversions
	#				print('ignoring')
				continue

			if int(s[1]) - int(s[0]) < min_aln_len:
		#		print('ignoring')
				continue

			if abs(int(s[0]) - int(prev_line[1])) < smooth_gap_width and \
				abs(int(s[2]) - int(prev_line[3])) < smooth_gap_width:
			#	print("joining")
				pass
			else:
	#				print("terminating")
				newline = aln_start_line
				newline[1] = prev_line[1]
				newline[3] = prev_line[3]
				if int(newline[1]) - int(newline[0]) > min_smoothed_aln_len:
	#				print("storing")
					smoothed_lines.append(newline)
				aln_start_line = s
			prev_line = s


#print("output")
#for s in smoothed_lines:
#	print(s)

if len(smoothed_lines) > 0:
	#is it a corner cutting parallel off-diagonal?
	if int(smoothed_lines[0][0]) < contig_edge_margin and int(smoothed_lines[0][3]) > max_tiglen - contig_edge_margin:
		if int(smoothed_lines[-1][2]) < contig_edge_margin and int(smoothed_lines[-1][1]) > max_tiglen - contig_edge_margin:
			#out.write("It's overcircularized!")

			out.write(tigname + ":" + smoothed_lines[0][0] + '-' + smoothed_lines[-1][0] + "\n")

out.write("done\n")
