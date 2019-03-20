#!/usr/bin/env bash

#Count kmers in assembly contigs, each individually

cut -f1 ../3c.circlator_final/ | xargs -n 1 -I foo sh -c "samtools faidx ../2.medaka/consensus.fasta foo > whole_tigs/foo.fa"
"cat ../3.circulator_canu/06.fixstart.fasta.fai | cut -f1 | " +
		"xargs -n1 -I foo -P 8 sh -c 'samtools faidx ../2.medaka/consensus.fasta foo > tigs/foo.fa'" +
		"jellyfish count -C whole_tigs/foo -m 6 -s 10000 -o jellyfish/db/foo.jf; " +
		"jellyfish dump jellyfish/db/foo.jf --tab --column | sort -k1,1 > jellyfish/whole_tsv/foo.tsv" +
		"join_several jellyfish/whole_tsv/* | bash > jellyfish/joined_whole.tsv" +
		"(head -1 jellyfish/joined_whole.tsv; sed '1d' jellyfish/joined_whole.tsv | sort -k1,1 )  > jellyfish/joined_whole_sorted.tsv"
