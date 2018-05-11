require(ggplot2)
require(grid)
require(doBy)


f = snakemake@input[[1]]
coords = read.table(f, sep = '\t', comment.char = '', quote = '', skip = 4)
colnames(coords) = c('Start.Query', 'End.Query', 'Start.Ref', 'End.Ref', 'Length.Query', 'Length.Ref', 'Identity', 'Query', 'Ref')

coords = coords[coords$Length.Query > 5000,]

#reverse all QUERIES that have mostly reverse-oriented hits in order to clarify layout
orientations = unique(data.frame(paste(coords$Query), sapply(coords$Query, function(x) sum((coords$Start.Query < coords$End.Query & coords$Start.Ref < coords$End.Ref)[coords$Query == x])/sum(coords$Query == x))))
names(orientations) = c('Query', 'Fraction.Forward')
toflip = orientations$Query[orientations$Fraction.Forward > 0.5]
#flip
for (flipped in toflip){
	max.coord = max(c(coords$Start.Query[paste(coords$Query)==flipped], coords$End.Query[paste(coords$Query)==flipped]))
	coords[coords$Query == flipped,]$Start.Query = max.coord - coords[coords$Query == flipped,]$Start.Query
	coords[coords$Query == flipped,]$End.Query = max.coord - coords[coords$Query == flipped,]$End.Query
}

#reverse all REFERENCES that have mostly reverse-oriented hits in order to clarify layout
orientations = unique(data.frame(paste(coords$Ref), sapply(coords$Ref, function(x) sum(((coords$Start.Query < coords$End.Query & coords$Start.Ref > coords$End.Ref) | (coords$Start.Query > coords$End.Query & coords$Start.Ref < coords$End.Ref))[coords$Ref == x])/sum(coords$Ref == x))))
names(orientations) = c('Ref', 'Fraction.Forward')
toflip = orientations$Ref[orientations$Fraction.Forward > 0.5]
#flip
for (flipped in toflip){
	max.coord = max(c(coords$Start.Ref[paste(coords$Ref)==flipped], coords$End.Ref[paste(coords$Ref)==flipped]))
	coords[coords$Ref == flipped,]$Start.Ref = max.coord - coords[coords$Ref == flipped,]$Start.Ref
	coords[coords$Ref == flipped,]$End.Ref = max.coord - coords[coords$Ref == flipped,]$End.Ref
}


#sort the query contigs by appearance in the reference alignment
coords$Query = factor(coords$Query)
centroids = sapply(coords$Query, function(x) median(coords$Start.Ref[coords$Query == x]))
coords = coords[order(coords$Ref, centroids, decreasing = F),]
coords$Query = factor(coords$Query, levels = unique(coords$Query))

#sort the reference contigs by appearance in the query alignment
coords$Ref = factor(coords$Ref)
centroids = sapply(coords$Ref, function(x) median(coords$Start.Query[coords$Ref == x]))
coords = coords[order(coords$Query, centroids, decreasing = F),]
coords$Ref = factor(coords$Ref, levels = unique(coords$Ref))

#format into ggplot-compatible form
toplot = data.frame(c(
	coords$Start.Query, coords$End.Query),
	c(coords$Start.Ref, coords$End.Ref),
	c(paste(coords$Query), paste(coords$Query)),
	c(paste(coords$Ref), paste(coords$Ref)),
	c(rownames(coords), rownames(coords)
))
names(toplot) = c('Query', 'Ref', 'Query.Name', 'Ref.Name', 'Line.Seg')

#apply ordering determined earlier
toplot$Query.Name = factor(toplot$Query.Name, levels = levels(coords$Query))
toplot$Ref.Name = factor(toplot$Ref.Name, levels = levels(coords$Ref))

#edit SPAdes-style contig names
levels(toplot$Query.Name) = gsub('_length.*', '', levels(toplot$Query.Name))

blarg = toplot

#placemarker

toplot = blarg
#modify alignments to make faceting unnecessary
#increment contig alignments along query axis
maxima = summaryBy(Query ~ Query.Name, toplot, FUN = max)
deltas = cumsum(maxima[,2]) - maxima[,2]
names(deltas) = maxima[,1]
toplot$Query = toplot$Query + deltas[toplot$Query.Name]

#increment contig alignments along reference axis
maxima = summaryBy(Ref ~ Ref.Name, toplot, FUN = max)
deltas = cumsum(maxima[,2]) - maxima[,2]
names(deltas) = maxima[,1]
toplot$Ref = toplot$Ref + deltas[toplot$Ref.Name]
require(scales)
dotplot.color.brewer = brewer_pal(palette = 'Set3', type = 'seq')
alt.color.brewer = brewer_pal(palette = 'Set1', type = 'seq')
alt.colors = alt.color.brewer(9)
dotplot.colors = dotplot.color.brewer(12)
dotplot.colors[2] = alt.colors[2] #get rid of light yellow
dotplot.colors = c(dotplot.colors, alt.colors, dotplot.colors, alt.colors)

ggplot(toplot, aes(x = Ref, y = Query, group = Line.Seg, colour = Query.Name)) +
	geom_path(size = 0.7) +
	scale_x_continuous(expand=c(0,0), breaks = seq(0, 20e6, 1e6), labels = paste(seq(0, 20, 1), '', sep = '')) +
	scale_y_continuous(expand=c(0,0), breaks = seq(0, 20e6, 1e6), labels = paste(seq(0, 20, 1), '', sep = '')) +
	scale_color_manual(values = dotplot.colors) +
	theme(legend.position="none") +
	xlab('Reference Assembly (Mb)') +
	ylab('Athena Assembly (Mb)') +
	theme(panel.spacing = grid::unit(0.75, "mm"),
				panel.background = element_blank(),
				panel.grid.major = element_line(colour = "grey",size=0.5, linetype = 'dotted'),
				axis.title.x = element_text(size = rel(0.8)),
				axis.title.y = element_text(size = rel(0.8))
				)# +
	#facet_grid(toplot$Query.Name~toplot$Ref.Name, scales = 'free', space = 'free')# + theme(strip.background = element_blank(),strip.text.x = element_blank(),strip.text.y = element_blank())

ggsave(snakemake@output[[1]], width = 4, height = 4)
