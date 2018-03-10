prokka = read.table(snakemake@input[[1]], sep = '\t', header = T, comment.char = '', quote = '', row.names = NULL)
quast = read.table(snakemake@input[[2]], sep = '\t', header = T, comment.char = '', quote = '', row.names = NULL)
checkm = read.table(snakemake@input[[3]], sep = '\t', header = T, comment.char = '', quote = '', row.names = NULL)
trna = read.table(snakemake@input[[4]], sep = '\t', header = T, comment.char = '', quote = '', row.names = NULL)
rrna = read.table(snakemake@input[[5]], sep = '\t', header = T, comment.char = '', quote = '', row.names = NULL)
classify = read.table(snakemake@input[[6]], sep = '\t', header = T, comment.char = '', quote = '', row.names = NULL)

out = prokka

out = merge(out, quast)

out = merge(out, checkm)

out = merge(out, trna)

out = merge(out, rrna)

out = merge(out, classify)

write.table(out, snakemake@output[[1]], sep = "\t", quote = F, row.names = F)
