require(circlize)
to_file = TRUE
#assembly or reference
base.base.dir = '~/scg4_moss/projects/10x/13.athena/circos/easy_10/'
for (f in list.files(base.base.dir, pattern = "kleb|*genomic.fna|*.merged.fa$")){
	print(f)
	base.dir = paste(base.base.dir, f, "/", sep = "")#'~/scg4_moss/projects/10x/13.athena/circos/easy_10/'
	#base.dir = '~/scg4_moss/projects/10x/13.athena/assembly_closure/klebsiella/'
	if(to_file){
		png(paste('~/Google Drive/Stanford/Bhatt lab/Bhatt Lab Drive/Project Collaboration/10X-GVHD Paper/Figures/Figure 2/', f, '.png', sep = ""), units = 'in', width = 8, height = 8, res = 300)
	}
	#rRNA
	bed16s = tryCatch(read.table(paste(base.dir, '16s_loci/16s.bed', sep = "")), error = function(e) '')
	bed23s = tryCatch(read.table(paste(base.dir, '16s_loci/23s.bed', sep = "")), error = function(e) '')

	#genome
	genome = read.table(paste(base.dir, 'merge.fasta.fai', sep = ""))#paste(base.dir, 'GCF_000163075.1_ASM16307v1_genomic.fna.fai', sep = ''))

	#contigs
	bed1.f = paste(base.dir, 'klebathena.bed', sep = '')
	bed2.f = paste(base.dir, 'klebref.bed', sep = '')
	bed3.f = paste(base.dir, 'nextera.bed', sep = '')

	#aesthetic knobs
	tick.interval = 0.5e6
	rRNA_width = 1000
	track.height = 0.3
	contig.alpha = 0.5
	bed1.color = rgb(0, 0.4, 0.5, alpha = contig.alpha)#bluish
	bed2.color = rgb(0.1, 0.1, 0.1, alpha = contig.alpha)#grey
	bed3.color = rgb(0.8, 0.1, 0, alpha = contig.alpha) #reddish

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
	letter_idx = 1
	alpha = c('A', 'B', 'C', 'D')

	contigs1 = read.table(bed1.f, sep = "")
	contigs1 = contigs1[contigs1[,3] - contigs1[,2] > min.contig.size, ]
	contigs1[,2] = contigs1[,2] + contig.shave.width
	contigs1[,3] = contigs1[,3] - contig.shave.width
	contigs2 = read.table(bed2.f, sep = "")
	contigs2 = contigs2[contigs2[,3] - contigs2[,2] > min.contig.size, ]
	contigs2[,2] = contigs2[,2] + contig.shave.width
	contigs2[,3] = contigs2[,3] - contig.shave.width
	contigs3 = read.table(bed3.f, sep = "")
	contigs3 = contigs3[contigs3[,3] - contigs3[,2] > min.contig.size, ]
	contigs3[,2] = contigs3[,2] + contig.shave.width
	contigs3[,3] = contigs3[,3] - contig.shave.width


	circos.genomicTrackPlotRegion(list(contigs3, contigs1, contigs2), bg.border = NA, stack = TRUE, track.height = track.height, ylim = c(0,0.5), panel.fun =function(region, value, ...) {
		i = getI(...)
		col = c(bed3.color, bed1.color, bed2.color)[i]
		circos.genomicRect(region, value, col = col, border = "white", ytop = i + contig_height, ybottom = i - contig_height, ...)
		circos.axis(h = "top", major.at = breaks, labels = c('', paste(breaks[2:length(breaks)]/1e6, '', sep ='')), minor.ticks = 0,
								major.tick.percentage = 0.2, labels.away.percentage = 0.05, labels.cex = 1)
	})


	#circos.genomicTrackPlotRegion(contigs3, bg.border = NA, stack = FALSE, track.height = 0.05, ylim = c(0,0.5), panel.fun =function(region, value, ...) {
	#	circos.genomicRect(region, value, col = bed3.color, border = "white", ytop = 0.75, ybottom = 0, ...)
	#})


		#INSERTIONS AND 16S
		if (bed16s != '' && bed23s != ''){

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

	#export as PDF 8in. x 8in.
	if(to_file){
		dev.off()
	}
}
