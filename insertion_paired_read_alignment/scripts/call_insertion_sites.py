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

for min_offset in [0, 1000]:
	if min_offset == 0:
		offset_outf = snakemake.output[1] #sys.argv[3]
		offset_out = open(offset_outf, 'w')

		bam_outf = snakemake.output[3] #sys.argv[4]
		bam_out = pysam.AlignmentFile(bam_outf, 'wb', template=samfile_in)

	if min_offset == 1000:
		offset_outf = snakemake.output[2] #sys.argv[3]
		offset_out = open(offset_outf, 'w')

		bam_outf = snakemake.output[4] #sys.argv[4]
		bam_out = pysam.AlignmentFile(bam_outf, 'wb', template=samfile_in)


	#check cross-correlations
	cors = all_cors[min_offset:2000] #generate cross-correlations for one-sided overlaps of the two vectors of depths

	#best offset is given by the degree of overlap between the
	#two matrices at each cross-correlation calculation.
	best_offset = numpy.argmax(cors) + min_offset #len(cors)/2 - numpy.argmax(cors) 
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
	offset_out.write(str(best_offset) + "\n")

	offset_out.close()
	bam_out.close()

	pysam.sort("-o", os.path.splitext(bam_outf)[0] + "_sorted.bam", bam_outf)
	pysam.index(os.path.splitext(bam_outf)[0] + "_sorted.bam")

samfile_in.close()
