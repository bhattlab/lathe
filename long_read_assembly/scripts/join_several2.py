#!/usr/bin/env python
import sys
import os

files = snakemake.input

basic_cmd = "join -e 0 -o auto -a1 -a2 {f0} {f1} -t '\t' "

cmd = 'echo ' + ' '.join(files) + "| tr ' ' '\\t' ;"
cmd += basic_cmd.format(f0=files[0], f1 = files[1]) 
remaining_joins = [basic_cmd.format(f0='-',f1=f) for f in files[2:]]
if len(remaining_joins) > 0:
	cmd += "|\n" + " |\n".join(remaining_joins)

cmd += ' > ' + snakemake.output[0]
os.system(cmd)