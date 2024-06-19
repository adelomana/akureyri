
# the objective of this script is to convert raw counts into TPM and ensembl ids.

# https://divingintogeneticsandgenomics.com/post/how-to-convert-raw-counts-to-tpm-for-tcga-data-and-make-a-heatmap-across-cancer-types/

# https://search.r-project.org/CRAN/refmans/DGEobj.utils/html/convertCounts.html


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("edgeR")

#
# 1. load libraries
#
library(DESeq2)
library(DGEobj.utils)

#
# 2. user-defined variables
#
counts_file = '/Users/adrian/research/akureyri/data/transcriptomics/GSE114007_raw_counts_normal.tsv'
annotation_file = '/Users/adrian/research/akureyri/results/deseq2_pipeline/annotation.tsv'

#
# 3. read information
# 

# 3.1. counts data
counts = read.table(counts_file, sep='\t', header=TRUE, row.names=1)
dim(counts)
View(counts)

# 3.2. read gene length
anno = read.table(annotation_file, 
                  sep='\t', 
                  header=TRUE, 
                  row.names=1,
                  quote="") # R surprises
dim(anno)
View(anno)

# 
# 4. analysis
#
gene_lengths = c()
found_symbols = c()
lost_symbols = c()
symbols = rownames(counts)

for (symbol in symbols){
  gene_length = anno[anno$external_gene_name == symbol, ]$geneLength[1]

  if (is.na(gene_length)) {
    lost_symbols = c(lost_symbols, symbol)
  } else {
    gene_lengths = c(gene_lengths, gene_length)
    found_symbols = c(found_symbols, symbol)
    }
} 

length(symbols)
length(found_symbols)
length(gene_lengths)
length(lost_symbols)

annotated_counts = counts[found_symbols, ]
annotated_counts_matrix = as.matrix(annotated_counts)
View(annotated_counts_matrix)

tpm = convertCounts(annotated_counts_matrix, unit='TPM', geneLength=gene_lengths)

write.table(tpm, 'converted_tpm.tsv', sep='\t', quote=FALSE)



