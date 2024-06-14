input_file = '/Users/adrian/research/akureyri/data/transcriptomics/GSE114007_raw_counts_normal.tsv'
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
