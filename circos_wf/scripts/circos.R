require(circlize)

output = snakemake@output[[1]]
pdf(output, width = 8, height = 8)#, res = 300, units = 'in', type = 'cairo')

#rRNA
#bed16s = tryCatch(read.table(paste(base.dir, '16s_loci/16s.bed', sep = "")), error = function(e) '')
#bed23s = tryCatch(read.table(paste(base.dir, '16s_loci/23s.bed', sep = "")), error = function(e) '')

#genome
genome = read.table(snakemake@input[[1]])

#aesthetic knobs
tick.interval = 0.5e6
rRNA_width = 1000
track.height = 0.1
contig.alpha = 0.5
bed1.color = rgb(0, 0.4, 0.5, alpha = contig.alpha)#bluish
bed2.color = rgb(0.1, 0.1, 0.1, alpha = contig.alpha)#grey

locus.alpha = 0.8
color.16s = rgb(1, 0.5, 0, alpha = locus.alpha)
color.23s = rgb(0, 0.5, 1, alpha = locus.alpha)
contig_height = 0.5
min.contig.size = 25000
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
for (pair in seq(4, length(snakemake@input), by = 2)){

	#contigs
	bed1.f = snakemake@input[[pair]]
	bed2.f = snakemake@input[[pair+1]]

	contigs1 = read.table(bed1.f, sep = "")
	contigs1 = contigs1[contigs1[,3] - contigs1[,2] > min.contig.size, ]
	contigs1[,2] = contigs1[,2] + contig.shave.width
	contigs1[,3] = contigs1[,3] - contig.shave.width
	contigs2 = read.table(bed2.f, sep = "")
	contigs2 = contigs2[contigs2[,3] - contigs2[,2] > min.contig.size, ]
	contigs2[,2] = contigs2[,2] + contig.shave.width
	contigs2[,3] = contigs2[,3] - contig.shave.width


	circos.genomicTrackPlotRegion(list(contigs1, contigs2), bg.border = NA, stack = TRUE, track.height = track.height, ylim = c(0,0.5), panel.fun =function(region, value, ...) {
		i = getI(...)
		col = c(bed1.color, bed2.color)[i]
		circos.genomicRect(region, value, col = col, border = "white", ytop = i + contig_height, ybottom = i - contig_height, ...)
		if (pair == 4){
			circos.axis(h = "top", major.at = breaks, labels = c('', paste(breaks[2:length(breaks)]/1e6, '', sep ='')), minor.ticks = 0,
								major.tick.percentage = 0.2, labels.away.percentage = 0.05, labels.cex = 1)
		}
	})
}

#INSERTIONS AND 16S
if (FALSE) { #}(bed16s != '' && bed23s != ''){

	#deduplicate beds
	if (nrow(bed16s) > 1){
		bed16s = bed16s[order(bed16s[,2]),]
		bed16s = bed16s[order(bed16s[,1]),]
		bad_idx = c(0)
		for (i in 1:(nrow(bed16s)-1)){
			if (bed16s[i,1] == bed16s[i+1,1] & bed16s[i,3] > bed16s[i+1,2] - 10000){
				bad_idx = c(bad_idx, i)
			}
		}
		bad_idx = bad_idx[-1]
		if (length(bad_idx) > 0){
			bed16s = bed16s[-c(bad_idx),]
		}
	}

	if (nrow(bed16s) > 1){
		bed23s = bed23s[order(bed23s[,2]),]
		bed23s = bed23s[order(bed23s[,1]),]
		bad_idx = c(0)
		for (i in 1:(nrow(bed23s)-1)){
			if (bed23s[i,1] == bed23s[i+1,1] & bed23s[i,3] > bed23s[i+1,2] - 10000){
				bad_idx = c(bad_idx, i)
			}
		}
		bad_idx = bad_idx[-1]
		if (length(bad_idx) > 0){
			bed23s = bed23s[-c(bad_idx),]
		}
	}

	bed16s = bed16s[order(bed16s[,2]),]
	bed16s = cbind(bed16s[,c(1,2)], bed16s[,2] + rRNA_width, color.16s)

	bed23s = bed23s[order(bed23s[,2]),]
	bed23s = cbind(bed23s[,c(1,2)], bed23s[,2] + rRNA_width, color.23s)

	bedlist = list(bed16s, bed23s)


	circos.genomicTrackPlotRegion(bedlist, stack = FALSE, bg.border = NA, track.height = track.height, ylim = c(0,1), panel.fun =function(region, value, ...) {
		i = getI(...)
		color = c(color.16s, color.23s)[i]
		circos.genomicRect(region, value, col = color, alpha = locus.alpha, border = color, ytop = 2.3, ybottom = 1.35, ...)
		circos.points(region[,1], rep((i/6) + 0.75, nrow(region)), pch = 16, col = color, cex = 1.2)
	})
}


circos.clear()
dev.off()
