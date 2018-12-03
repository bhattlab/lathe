#!/usr/bin/env python

import sys

# no fancy argument parsing, but ssub in scg_tools could be used as a model for how to do this

# first argument: the output of show-coords -T filter.delta
alignments = [l.strip().split("\t") for l in open(sys.argv[1]).readlines()]

# second argument: assembly1.fa.fai (produced by samtools faidx assembly2.fa)
t1 =[l.strip().split("\t") for l in open(sys.argv[2]).readlines()]

# third argument: assembly2.fa.fai (produced by samtools faidx assembly2.fa)
t2 =[l.strip().split("\t") for l in open(sys.argv[3]).readlines()] 

# prepare a dictionary for each of the two fasta indices with
# contig names as keys and contig lengths as values
t1d = {t[0]: int(t[1]) for t in t1}
t2d = {t[0]: int(t[1]) for t in t2}

# how far away from the edge of the contig (expressed as a fraction of contig length)
# can an alignment end and still be considered to extend to the edge? 
margin = 0.1 # int(sys.argv[4]) # provide this as a command line argument instead



def within_range_start(x, l, m):
	# x: location of an alignment endpoint within a contig
	# l: length of the contig
	# m: margin, as a fraction of the contig length, for considering a coordinate within range of the edge
	if x < (l * m):
		return(True)
	return(False)

def within_range_end(x, l, m):
	# see description of within_range_start.  The only difference is (1-m) in place of m; this is because we are
	# asking if the coordinate is within range of the opposite end of the contig.
	if x > l * (1 - m):
		return(True)
	return(False)

def is_spanning(a1, a2, b1, b2, la, lb, m):
	# a contig is defined as spanning if each end is within the margin of at least one side
	if (within_range_start(a1, la, m) or within_range_start(b1, lb, m)) and \
		(within_range_end(a2, la, m) or within_range_end(b2, lb, m)):
		return(True)


for aln in alignments:
	#parse the fields in each alignment in the input into separate variables
	tig1 = aln[7]
	tig2 = aln[8]

	l1 = int(t1d[tig1])
	l2 = int(t2d[tig2])

	x1 = min([int(aln[0]), int(aln[1])])
	x2 = max([int(aln[0]), int(aln[1])])
	
	y1 = min(int(aln[2]), int(aln[3]))
	y2 = max(int(aln[2]), int(aln[3]))

	
	#print the contig pair if one spans the other
	if is_spanning(x1, x2, y1, y2, l1, l2, margin):
		print('\t'.join([tig1, tig2]))