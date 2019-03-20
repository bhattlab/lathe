#!/usr/bin/env python

import fileinput
import os

prev = []
prev_is_terminal = False
verbose = True

margin = snakemake.params['margin']

with open(snakemake.output[0], 'w') as oot:
	for l in open(snakemake.input[0]):
		s = l.strip().split("\t")
		if prev == []:
			next
		if s == ['no circularizations']:
			oot.write("no circularizations\n")
			break

		target_tig = s[9]
		spanner = s[10]

		#test for terminal alignment
		r_start = min([int(i) for i in s[0:2]])
		r_end = max([int(i) for i in s[0:2]])
		q_start = min([int(i) for i in s[2:4]])
		q_end = max([int(i) for i in s[2:4]])

		r_len = int(s[7])
		q_len = int(s[8])

		if r_start < margin or r_end > (r_len - margin): #ref is terminal
			if q_start < margin or q_end > (q_len - margin): #query is terminal
				if prev_is_terminal and spanner == prev[10] and target_tig == prev[9]: #spanned event

					#trim when the center of the spanning contig aligns to both ends of the genome (consistent with overcircularization)
					if q_start < prev_q_end:
						trim_bases = prev_q_end - q_start
			#			oot.write('Need to trim ' + str(trim_bases) + " bases.\n")
						oot.write("{c}:1-{e}\n".format(c = target_tig, e = str(r_len-trim_bases)))

					#extend with spanning contig when it closes a gap
					else:
						insert_bases = [q_start, prev_q_end]
						insert_bases.sort()
						insert_string = spanner + ":" + '-'.join([str(i) for i in insert_bases])
						#print(insert_string)
						oot.write(target_tig + "\n")
						oot.write(insert_string + "\n")

					prev_is_terminal = False
					next

				prev_is_terminal=True
			else:
				prev_is_terminal = False
		else:
			prev_is_terminal = False

		#assign prevs
		prev = s
		prev_r_start = r_start
		prev_r_end   = r_end
		prev_q_start = q_start
		prev_q_end = q_end

		prev_r_len = r_len
		prev_q_len = q_len

	oot.write('done\n')
