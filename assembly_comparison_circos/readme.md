# Circos Assembly Comparison Visualization Workflow

# Inputs

### Alter config.yaml to provide the following:
 * **Reference**: Sequence against which all others will be aligned. Forms the "backbone" of the circos plot. Multi-fasta.

 * **Dark contigs**: contigs which will be shaded darker in the final plot. Leave empty for all dark. Two tab-delimited columns: assembly file (matches the first column of Assemblies below), contig name.

 * **Highlight sequences**: Sequence(s) which will be aligned against the reference and called out with a triangle in the innermost track. Multi-fasta. Highlights can be colored by group according to the first ';'-delimited sequence in the fasta sequence identifiers. Only two groups supported.

 * **Assemblies**: List of assemblies annotated with condition and time point from which to find contigs mapping to the chosen ref. Three tab delimited columns: assembly, condition, timepoint. Conditions and timepoints can be any string, so long as it matches that from other assemblies of the same group. No underscores allowed! The order of conditions (from the bottom of each track to the top) will reflect what is found in this file. The time points are sorted alphabetically.

 * **Highlight intensities**: <# time points + 1> tab-delimited columns.  First is highlight sequence identifiers.  Second through last are decimals ranging 0 to 1 inclusive indicating the strength of highlight desired. Filled and empty circles are used for values above and below 0.5.  Values below 0 receive an 'x'.
