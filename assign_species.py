#!/usr/bin/env python

import sys
import operator
import os

#input: 6753_d/classify/6753_d.tsv, 6753_d/bins/bin.1.fa.fai
#    output: 6753_d/classify/bin_species_calls.tsv

krakf = snakemake.input[0]
binfaifs = os.listdir(os.path.dirname(snakemake.input[1]))
binfolder = os.path.dirname(snakemake.input[1])
outf = snakemake.output[0]

winning_margin = 60 #int(snakemake.input[1])

tig_species = {}
called_species = []

with open(krakf, 'r') as krak:
	with open(outf, 'w') as out:
		for binfaif in binfaifs:
			if os.path.splitext(binfaif)[-1] != '.fai':
				continue
			species_votes = {}

			with open(os.path.join(binfolder, binfaif), 'r') as binfai:
				for l in krak.readlines():
					s = l.rstrip().split("\t")
					tig_species[s[0]] = '_'.join(s[1].split('|')[-1].split('_')[2:4])
				for l in binfai.readlines():
					s = l.rstrip().split("\t")
					if s[0] in tig_species:
						species_vote = tig_species[s[0]]
					else:
						species_vote = 'unclassified'
					votes = int(s[1])

					species_votes.setdefault(species_vote, 0)
					species_votes[species_vote] = species_votes[species_vote] + votes

			verbose = False
			if verbose:

				for idx in species_votes:
					print(idx + "\t" + str(
						round(100*species_votes[idx] /
							float(
								sum(
									species_votes.values()
									)
								)
							, 3)
						)
					)

			winning_species = max(species_votes.items(), key=operator.itemgetter(1))[0]
			winning_fraction = 100*max(species_votes.values())/float(sum(species_votes.values()))
			total_votes = sum(species_votes.values())
			bin = binfaif.split("/")[-1].replace('.fai', '')

			if winning_fraction > winning_margin:
				final_class = winning_species
			else:
				final_class = winning_species.split('_')[0]


			tmp = final_class
	#		final_class += "_" + str(called_species.count(final_class))
			called_species.append(tmp)

			out.write('\t'.join([
				bin,
				winning_species,
				str(round(winning_fraction, 3)),
				str(round(total_votes/float(1000000),2)),
				final_class
			]) + "\n")
