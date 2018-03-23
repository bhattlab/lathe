require(ggplot2)
require(doBy)
require(RColorBrewer)
require(scales)
require(plyr)

na.rm = function(x) {
	return(x[!is.na(x)])
}
proper=function(s) sub("(.)", ("\\U\\1"), tolower(s), pe=TRUE)

#s = read.table('~/scg4_moss/projects/10x/9.classify/classify/class_long.tsv', comment.char = '', quote = '', sep = '\t')
s = read.table(snakemake@input[[1]], comment.char = '', quote = '', sep = '\t')
colnames(s) = c('Sample', 'Condition', 'Timepoint', 'Percent', 'Reads', 'Direct Reads', 'Tax', 'Taxid')

s$Taxid = sapply(s$Taxid, function(x) paste(strsplit(as.character(x), '_')[[1]][1:min(c(2, length(strsplit(as.character(x), '_')[[1]])))], collapse = ' '))
#taxonomic level level
s = s[s$Tax == snakemake@params[['taxlevel']] ,]

#capitalize timepoints
levels(s$Timepoint) = toupper(levels(s$Timepoint))
#capitalize conditions
levels(s$Condition) = proper(levels(s$Condition))
levels(s$Condition)[3] = gsub('Shortread', 'Short Read', levels(s$Condition)[3])

#remove human
s = s[s$Taxid != 'Homo_sapiens',]
s = s[s$Taxid != 'Homo',]

collapsed = summaryBy(Reads ~ Sample + Condition + Timepoint + Taxid, s, FUN = sum)
collapsed$Relative.Abundance = ddply(collapsed, c('Sample', 'Condition', 'Timepoint', 'Taxid'), function(x) x$Reads.sum/sum(collapsed$Reads.sum[collapsed$Sample == x[1,1]]))[,5]

toplot = collapsed

names(toplot)[ncol(toplot)] = 'Relative Abundance'
abundance.threshold = sort(toplot$`Relative Abundance`, decreasing = T)[11]

toplot = cbind(toplot, toplot$`Relative Abundance` > abundance.threshold)
#toplot = toplot[toplot$`Relative Abundance` > 0.001,]
#toplot$`Relative Abundance` = toplot$`Relative Abundance`/sum(toplot$`Relative Abundance`)

toplot$Taxid = factor(toplot$Taxid)
toplot$Taxid = factor(toplot$Taxid, levels = unique(toplot$Taxid[order(toplot$`Relative Abundance`, decreasing = T)]))

colors = c(brewer.pal(12, 'Set3'), rep(c('#C0C0C0', '#DCDCDC'), length(toplot$Taxid)))
colors[9] = '#CD6155'
ggplot(toplot) + geom_bar(aes(x = Condition, y = `Relative Abundance`, fill = Taxid), stat = 'identity') + scale_fill_manual(values = colors, breaks = levels(toplot$Taxid)[1:12], name = '') + ylim(0, 1)  + facet_grid(.~Timepoint) + theme(axis.text.x = element_text(angle = 30, hjust = 1))


#+ theme(axis.text.x = element_blank())
if (snakemake@params[['taxlevel']] == 'S'){
		ggsave('plot/species_composition.pdf', width = 8.5, height = 4)
}
if (snakemake@params[['taxlevel']] == 'G'){
	ggsave('plot/genus_composition.pdf', width = 8, height = 4)
}