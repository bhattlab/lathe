
######## Snakemake header ########
import sys; sys.path.append("/home/elimoss/tools/anaconda3/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X \x00\x00\x006753_c/is_paired_to_assembly.bamq\x06X$\x00\x00\x006753_c/is_paired_to_assembly.bam.baiq\x07e}q\x08X\x06\x00\x00\x00_namesq\t}q\nsbX\x06\x00\x00\x00outputq\x0bcsnakemake.io\nOutputFiles\nq\x0c)\x81q\r(X)\x00\x00\x006753_c/peakcalling/cross_correlations.tsvq\x0eX\x1d\x00\x00\x006753_c/peakcalling/offset.txtq\x0fX#\x00\x00\x006753_c/peakcalling/read_shifted.bamq\x10X*\x00\x00\x006753_c/peakcalling/read_shifted_sorted.bamq\x11X.\x00\x00\x006753_c/peakcalling/read_shifted_sorted.bam.baiq\x12X$\x00\x00\x006753_c/peakcalling/offset_over1k.txtq\x13X*\x00\x00\x006753_c/peakcalling/read_shifted_over1k.bamq\x14X1\x00\x00\x006753_c/peakcalling/read_shifted_over1k_sorted.bamq\x15X5\x00\x00\x006753_c/peakcalling/read_shifted_over1k_sorted.bam.baiq\x16e}q\x17h\t}q\x18sbX\x06\x00\x00\x00paramsq\x19csnakemake.io\nParams\nq\x1a)\x81q\x1b}q\x1ch\t}q\x1dsbX\t\x00\x00\x00wildcardsq\x1ecsnakemake.io\nWildcards\nq\x1f)\x81q X\x06\x00\x00\x006753_cq!a}q"(h\t}q#X\x06\x00\x00\x00sampleq$K\x00N\x86q%sX\x06\x00\x00\x00sampleq&h!ubX\x07\x00\x00\x00threadsq\'K\x01X\t\x00\x00\x00resourcesq(csnakemake.io\nResources\nq))\x81q*(K\x01K\x01e}q+(h\t}q,(X\x06\x00\x00\x00_coresq-K\x00N\x86q.X\x06\x00\x00\x00_nodesq/K\x01N\x86q0uh-K\x01h/K\x01ubX\x03\x00\x00\x00logq1csnakemake.io\nLog\nq2)\x81q3}q4h\t}q5sbX\x06\x00\x00\x00configq6}q7(X\n\x00\x00\x00samplenameq8X\x06\x00\x00\x006753_cq9X\x06\x00\x00\x00reads1q:X\x1a\x00\x00\x00short_reads/6753_c_1.fq.gzq;X\x06\x00\x00\x00reads2q<X\x1a\x00\x00\x00short_reads/6753_c_2.fq.gzq=X\x0f\x00\x00\x00insertion_fastaq>X\x11\x00\x00\x00two_insertions.faq?X\x0e\x00\x00\x00assembly_fastaq@X\x1f\x00\x00\x00assembly/chromium_c_bins_4_6.faqAuX\x04\x00\x00\x00ruleqBX\n\x00\x00\x00peak_shiftqCub.'); from snakemake.logging import logger; logger.printshellcmds = True
######## Original script #########
import sys
import pysam
import numpy
import os

inputf = snakemake.input[0] #sys.argv[1]
samfile_in = pysam.AlignmentFile(inputf, 'rb')

cors_outf = snakemake.output[0] #sys.argv[2]
cors_out = open(cors_outf, 'w')


fwd_depths = []
rev_depths = []

#extract coverages by read orientation
print("Counting stranded coverage")
for pup in samfile_in.pileup(): #walk across positions with aligned reads
	print(pup)
	fwd = [r for r in pup.pileups if not r.alignment.is_reverse]
	fwd_depths.append(len(fwd))
	rev = [r for r in pup.pileups if r.alignment.is_reverse]
	rev_depths.append(len(rev))

print("Calculating cross-correlations")
all_cors = numpy.correlate(fwd_depths, rev_depths, "same")
cors_out.write("\n".join([str(i) for i in all_cors])) #store the cross-correlations to a file
cors_out.close()
samfile_in.close()

for min_offset in [0, 1000]:
	if min_offset == 0:
		offset_outf = snakemake.output[1] #sys.argv[3]
		offset_out = open(offset_outf, 'w')

		bam_outf = snakemake.output[2] #sys.argv[4]
		bam_out = pysam.AlignmentFile(bam_outf, 'wb', template=samfile_in)

	if min_offset == 1000:
		offset_outf = snakemake.output[5] #sys.argv[3]
		offset_out = open(offset_outf, 'w')

		bam_outf = snakemake.output[6] #sys.argv[4]
		bam_out = pysam.AlignmentFile(bam_outf, 'wb', template=samfile_in)


	#check cross-correlations
	cors = all_cors[min_offset:2000] #generate cross-correlations for one-sided overlaps of the two vectors of depths

	#best offset is given by the degree of overlap between the
	#two matrices at each cross-correlation calculation.
	best_offset = len(cors)/2 - numpy.argmax(cors)
	best_offset = best_offset/2 #divide by two to get the amount by which
								#forward and reverse reads should each be shifted


	print("Outputting shifted read BAM files")
	for read in samfile_in.fetch(): #walk across one bp at a time
		if read.is_reverse:
			read.reference_start = max([0,read.reference_start - best_offset])
		else:
			read.reference_start = read.reference_start + best_offset

		bam_out.write(read)

	#store the best offset to a file
	offset_out.write(str(best_offset))

	offset_out.close()
	bam_out.close()

	pysam.sort("-o", os.path.splitext(bam_outf)[0] + "_sorted.bam", bam_outf)
	pysam.index(os.path.splitext(bam_outf)[0] + "_sorted.bam")
