#!/usr/bin/env python

import sys

t1 =[l.strip().split("\t") for l in open(sys.argv[2]).readlines()]
t2 =[l.strip().split("\t") for l in open(sys.argv[3]).readlines()] 

t1d = {t[0]: int(t[1]) for t in t1}
t2d = {t[0]: int(t[1]) for t in t2}


alignments = [l.strip().split("\t") for l in open(sys.argv[1]).readlines()]

margin = 0.1 #int(sys.argv[4]) #10000 is a good choice


def within_range_start(x, l, m):
	if x < (l * m):
		return(True)
	return(False)

def within_range_end(x, l, m):
	if x > l * (1 - m):
		return(True)
	return(False)

def is_spanning(a1, a2, b1, b2, la, lb, m):
	if (within_range_start(a1, la, m) or within_range_start(b1, lb, m)) and \
		(within_range_end(a2, la, m) or within_range_end(b2, lb, m)):
		return(True)


for aln in alignments:
	tig1 = aln[7]
	tig2 = aln[8]

	l1 = int(t1d[tig1])
	l2 = int(t2d[tig2])

	x1 = min([int(aln[0]), int(aln[1])])
	x2 = max([int(aln[0]), int(aln[1])])
	
	y1 = min(int(aln[2]), int(aln[3]))
	y2 = max(int(aln[2]), int(aln[3]))

	if is_spanning(x1, x2, y1, y2, l1, l2, margin):
		print('\t'.join([tig1, tig2]))