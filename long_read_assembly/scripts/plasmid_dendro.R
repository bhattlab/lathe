#!/usr/bin/env R

#require(vegan)
require(plyr)
require(ggplot2)
require(factoextra)
require(reshape2)
require(stringr)

format.table = function(t){
	t.collapsed = data.frame(t(t))
	t.collapsed$tig = row.names(t.collapsed)
	t.collapsed = data.frame(merge(m, t.collapsed, all = T))
	t.collapsed = t.collapsed[complete.cases(t.collapsed),]
	t.collapsed = ddply(data.frame(t.collapsed), .(bin), function(x) apply(x[3:ncol(x)], 2, sum))
	row.names(t.collapsed) = t.collapsed$bin
	t.collapsed = t.collapsed[,-1]
	t.collapsed
}

t1 = read.table(          snakemake@input[[1]], sep = '\t', header = T, row.names = 1)
t2 = read.table(          snakemake@input[[2]], sep = '\t', header = T, row.names = 1)
t3 = read.table(          snakemake@input[[3]], sep = '\t', header = T, row.names = 1)
m  = read.table(          snakemake@input[[4]], sep = '\t')
plasmids = read.table(    snakemake@input[[5]])
species_assn = read.table(snakemake@input[[6]], header = T)


#knobs
min_bin_size = 0.5
min_plasmid_size = 40e3

colnames(m) = c('bin', 'tig')
colnames(t1) = str_extract(colnames(t1), 'tig[0-9tig_]{8,21}')
colnames(t2) = str_extract(colnames(t2), 'tig[0-9tig_]{8,21}')
colnames(t3) = str_extract(colnames(t3), 'tig[0-9tig_]{8,21}')
colnames(plasmids) = c('tig', 'size', 'myeh', 'meh', 'nah')
row.names(species_assn) = species_assn$Bin

#filter things
plasmids = plasmids[plasmids$size > min_plasmid_size,]
species_assn = species_assn[species_assn$Size.Mb > min_bin_size,]
m = m[m$bin %in% species_assn$Bin,]

plasmids.t1 = t1[,which(colnames(t1) %in% plasmids$tig)]
plasmids.t2 = t2[,which(colnames(t2) %in% plasmids$tig)]
plasmids.t3 = t3[,which(colnames(t3) %in% plasmids$tig)]

t1 = t1[,which(!colnames(t1) %in% plasmids$tig)]
t2 = t2[,which(!colnames(t2) %in% plasmids$tig)]
t3 = t3[,which(!colnames(t3) %in% plasmids$tig)]

t1.collapsed = format.table(t1)
t2.collapsed = format.table(t2)
t3.collapsed = format.table(t3)

#add circular contigs back in

t1.collapsed = rbind(t1.collapsed, t(plasmids.t1))
t2.collapsed = rbind(t2.collapsed, t(plasmids.t2))
t3.collapsed = rbind(t3.collapsed, t(plasmids.t3))

t1.m = melt(cbind(rownames(t1.collapsed), t1.collapsed))
colnames(t1.m) = c('Bin', 'Kmer', 'Count.5mC')
t2.m = melt(cbind(rownames(t2.collapsed), t2.collapsed))
colnames(t2.m) = c('Bin', 'Kmer', 'Count.6mA')
t3.m = melt(cbind(rownames(t3.collapsed), t3.collapsed))
colnames(t3.m) = c('Bin', 'Kmer', 'Count.Tig')

#pm

t.m = merge(t1.m, t2.m, all = T)
t.m = merge(t.m, t3.m, all = T)
t.m[is.na(t.m)] = 0

#normalize by total kmers found in the same sequence collection
t.m$Count.5mC = t.m$Count.5mC/(t.m$Count.Tig + 1)
t.m$Count.6mA = t.m$Count.6mA/(t.m$Count.Tig + 1)

t.m.5 = t.m[,1:3]
colnames(t.m.5)[3] = 'Count'
t.m.5$Kmer = paste(t.m.5$Kmer, '.m5C', sep = '')
t.m.6 = t.m[,c(1,2,4)]
colnames(t.m.6)[3] = 'Count'
t.m.6$Kmer = paste(t.m.6$Kmer, '.m6A', sep = '')

t.m2 = rbind(t.m.5, t.m.6)

t.norm = t(reshape(t.m2, direction = 'wide', timevar = 'Bin', idvar = 'Kmer')[,-1])
t.norm = t.norm[,apply(t.norm, 2, sum) > 0]
rownames(t.norm) = gsub('Count.', '', rownames(t.norm))
d = dist(t.norm, method = 'euclidean', diag=T, upper =T) #euclidean jaccard
h = hclust(d)
new_labs = apply(species_assn[h$labels,2:3], 1, function(x) paste(rev(x), collapse = '_'))
new_labs = new_labs[!new_labs == 'NA_NA']
h$labels[1:length(new_labs)] = paste(h$labels[1:length(new_labs)], new_labs)

fviz_dend(h, k = nrow(t.norm)/2, color_labels_by_k = F, horiz=T, palette='Set3', labels_track_height = 2)
ggsave(snakemake@output[[1]], width = 12, height = 10)
