# workflows
Workflows for bioinformatic data processing and analysis, with heavy emphasis on metagenomic DNA sequencing.

## bin_wf

Snakemake workflow for aligning, binning, classifying and evaluating a
metagenomic assembly.

## circos_wf
Snakemake workflow for visualizing assemblies of a particular genome across conditions and time points.  Calls out pre-identified sequences, highlights selected contigs.

#### Inputs

	**Reference**: Sequence against which all others will be aligned. Forms the "backbone" of the circos plot. Multi-fasta.

	**Dark contigs**: contigs which will be shaded darker in the final plot. Leave empty for all dark. Two tab-delimited columns: assembly name, contig name.

	**Highlight sequences**: Sequence(s) which will be aligned against the reference and called out with a triangle in the innermost track. Multi-fasta.

	**Assemblies**: List of assemblies annotated with condition and time point from which to find contigs mapping to the chosen ref. Three tab delimited columns: assembly, condition, timepoint. Conditions and timepoints can be any string, so long as it matches that from other assemblies of the same group.

