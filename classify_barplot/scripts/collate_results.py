import os

inputs = snakemake.input
output = []

sample_condtp = {}
with open(snakemake.config['datasets']) as asmf:
	basedir = os.path.dirname(snakemake.config['datasets'])
	for l in asmf.readlines():
		s = l.strip().split("\t")
		if len(s) == 1 or s[0] == 'Sample':
			continue

		sample = s[0]
		timepoint = s[1]
		condition = s[2]

		sample_condtp[sample] = "\t".join([condition, timepoint])


for inf in inputs:
	with open(inf, 'r') as i:
		for l in i.readlines():
			s = l.split("\t")
			sample = os.path.splitext(os.path.basename(inf))[0]
			newline = [sample, sample_condtp[sample]] + s[1:5] + ["_".join(s[6:])]
			output.append(newline)

with open(snakemake.output[0], 'w') as out:
	out.write("\n".join(["\t".join(o) for o in output]))
