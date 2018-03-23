#!/usr/bin/env python

import fileinput
import sys
import os


max_gap_width = int(snakemake.params['gap_width'])
input_coords_files = snakemake.input

for input_coords_file in input_coords_files:
	with open(input_coords_file) as input_coords:
		outputf = os.path.splitext(input_coords_file)[0] + '.bed'
		with open(outputf, 'w') as out:
			prev = []
			for l in input_coords.readlines():
				s = l.strip().split("\t")
				if len(s) != 9:
					continue
				ref_start = int(s[0])
				ref_end = int(s[1])
				query_start = int(s[2])
				query_end = int(s[3])
				ref = s[7]
				query = s[8]
			#	print(s)
				if len(prev) == 0:  #no previous contig to merge with
					prev.append(s)
			#		print('appended')
				elif prev[-1][7:9] == s[7:9] and \
					abs(ref_start - int(prev[-1][1])) <= max_gap_width and \
						abs(query_start - int(prev[-1][3])) <= max_gap_width:
						prev.append(s)
			#			print("small gap")
			#			print('appended')
			#		else:
			#			print(str(ref_start - int(prev[-1][1])))
			#			print(prev[-1][1])
			#			print("large gap encountered in " + query)
			#			pass
			#
			#			toprint = prev[0]
			#			toprint[1] = prev[-1][1]
			#			print("\t".join(toprint))
			#			prev = [s]

				else:# prev[-1][-1] != s[-1]: #merge and begin again
			#		print('output!')
					toprint = prev[0]
					toprint[1] = prev[-1][1]
#					print("\t".join(toprint))
					out.write("\t".join([toprint[i] for i in [7,0,1,8]]) + "\n")
					prev = [s]
			#		print("\n")

			toprint = prev[0]
			toprint[1] = prev[-1][1]

#			print("\t".join(toprint))
			out.write("\t".join([toprint[i] for i in [7,0,1,8]]) + "\n")
			prev = [s]
