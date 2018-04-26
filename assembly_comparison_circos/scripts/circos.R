require(circlize)

output = snakemake@output[[1]]
pdf(output, width = 8, height = 8)#, res = 300, units = 'in', type = 'cairo')

dark.contigs = read.table(snakemake@input[[2]])
colnames(dark.contigs) = c('assembly', 'condition', 'timepoint', 'contig')

#rRNA
highlight = read.table(snakemake@input[[3]])
highlight.intensities = read.table(snakemake@input[[4]])

#contig alignments
contig.alignments = unlist(snakemake@input)[5:length(snakemake@input)]
timepoint = sapply(contig.alignments, function(x) strsplit(strsplit(basename(x), '\\.')[[1]][1], '_')[[1]][2])
condition = sapply(contig.alignments, function(x) strsplit(strsplit(basename(x), '\\.')[[1]][1], '_')[[1]][1])

contig.alignments = data.frame(contig.alignments, timepoint, condition)

dark.contigs = merge(dark.contigs, contig.alignments)


#genome
genome = read.table(snakemake@input[[1]])

#aesthetic knobs
highlight_size = 1.2
tick.interval = 0.5e6
track.height = 0.1
contig.alpha = 0.3
bed1.color = adjustcolor("#409159", alpha = contig.alpha)
bed2.color = adjustcolor("#E69F00", alpha = contig.alpha)
bed3.color = adjustcolor("#5AC6E2", alpha = contig.alpha)
bed.colors.faint = c(bed1.color, bed2.color, bed3.color)

bed1.color = adjustcolor("#409159", alpha = 0.8)
bed2.color = adjustcolor("#E69F00", alpha = 0.9)
bed3.color = adjustcolor("#5AC6E2", alpha = 0.9)
bed.colors.bold = c(bed1.color, bed2.color, bed3.color)


locus.alpha = 1

color.highlight = rgb(1, 0.3, 0.1, alpha = locus.alpha)
color.highlight2 = rgb(0, 0.5, 1, alpha = locus.alpha)
highlight.colors = c(color.highlight, color.highlight2)

contig_height = 0.5
min.contig.size = 0
contig.shave.width = 1000

#init
genome = genome[genome[,2] > min.contig.size,]
par(mar =c(2, 2, 2, 2))
circos.par("track.height" = track.height,canvas.xlim =c(-1, 1), canvas.ylim =c(-1, 1), "start.degree" = 90)
cytoband = data.frame(genome[,1], rep(0, nrow(genome)), genome[,2], genome[,1], rep('foo', nrow(genome)))
cytoband[,1] = as.character(cytoband[,1])
circos.initializeWithIdeogram(cytoband=cytoband, plotType = NULL)
genome.length = sum(cytoband[,3])
breaks = seq(0, genome.length, by = tick.interval)

#CONTIGS
for (tp in unique(timepoint)){

	#contigs
	beds = contig.alignments$contig.alignments[contig.alignments$timepoint == tp]

	contig.list = lapply(beds, function(bed){
		contigs1 = read.table(as.character(bed), sep = "")
		contigs1 = contigs1[contigs1[,3] - contigs1[,2] > min.contig.size, ]
		contigs1[,2] = contigs1[,2] + contig.shave.width
		contigs1[,3] = contigs1[,3] - contig.shave.width
		contigs1
	})

	circos.genomicTrackPlotRegion(contig.list, bg.border = NA, stack = TRUE, track.height = track.height, ylim = c(0,0.5), panel.fun =function(region, value, ...) {
		i = getI(...)
		col = bed.colors.faint[i]
		circos.genomicRect(region, value, col = col, border = "white", ytop = i + contig_height, ybottom = i - contig_height, ...)

		#endarken darker contigs
		to_darken = which(as.character(value[,1]) %in% as.character(dark.contigs$contig[dark.contigs$timepoint == tp]))
		region.dark = region[to_darken,]
		value.dark = value[to_darken,]

		col = bed.colors.bold[i]
		if (length(value.dark) > 0){
			circos.genomicRect(region.dark, value.dark, col = col, border = "white", ytop = i + contig_height, ybottom = i - contig_height, ...)
		}

		#draw axis on the first track
		if (tp == unique(timepoint)[1]){
			circos.axis(h = "top", major.at = breaks, labels = c('', '', '1', '', '2', '', '3', '', '4', '', '5', '', '6'), minor.ticks = 0,
								major.tick.percentage = 0.2, labels.away.percentage = 0.05, labels.cex = 1.5)
		}
	})
}

#deduplicate highlight sequence beds
#highlight_interval = 30000
#breaks = seq(1, max(highlight[,2]), highlight_interval)
#highlight[,2] = breaks[cut(highlight[,2], breaks, labels = FALSE)+1]
#highlight[,3] = highlight[,2] + 1000
#highlight = unique(highlight)

highlight = highlight[order(highlight[,2]),]
highlight = cbind(highlight[,c(1,2)], highlight[,4])


if (length(strsplit(as.character(highlight[1,3]), ';')[[1]]) > 1){
	#color the highlights by the first ';' delimited portion of their original sequence names.
	highlight$type = sapply(highlight[,3], function(x) strsplit(as.character(x), ';')[[1]][1])
} else {
	#don't.
	highlight$type = rep('foo', nrow(highlight))
}


print(head(highlight))
print(head(highlight.intensities))
print(timepoint)
colnames(highlight) = c('Ref.contig', 'Coord', 'Highlight.seq', 'Group.member')
colnames(highlight.intensities) = c('Highlight.seq', unique(sort(timepoint)))

#intensities over 1 are reduced to 1
highlight.intensities[,-c(1)][highlight.intensities[,-c(1)]>1] = 1

#HIGHLIGHTED SEQUENCES
for (tp in sort(unique(timepoint))){
	highlight.merge = merge(highlight, highlight.intensities)
	highlight.merge = highlight.merge[,c(-1)]
	highlight.merge = highlight.merge[c(1,2,which(colnames(highlight.merge) == tp),3)]

	highlight.merge$Group.member = factor(highlight.merge$Group.member)
	levels(highlight.merge$Group.member) = c(TRUE, FALSE)

	bedlist = list(
		highlight.merge[highlight.merge$Group.member == TRUE,],
		highlight.merge[highlight.merge$Group.member == FALSE,]
	)

	circos.genomicTrackPlotRegion(bedlist, stack = FALSE, bg.border = NA, track.height = track.height, ylim = c(0,1), panel.fun =function(region, value, ...) {
		i = getI(...)
		color = highlight.colors[i]

		data = region[region[,2] >= 0,]
		circos.points(data[,1], rep(length(unique(timepoint)) + 2.5, nrow(data)), pch = 1, col = color, cex = highlight_size)
		circos.points(data[,1], rep(length(unique(timepoint)) + 2.5, nrow(data)), pch = 16, col = color, cex = (data[,2]) * highlight_size)

		no_data = region[region[,2] < 0,]
		circos.points(no_data[,1], rep(length(unique(timepoint)) + 2.5, nrow(no_data)), pch = 4, col = color, cex = highlight_size)
	})
}


circos.clear()
dev.off()
