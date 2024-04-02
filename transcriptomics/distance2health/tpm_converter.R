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
counts_file = '/Users/adrian/research/akureyri/data/GSE114007_normal_normalized.counts.txt'
annotation_file = '/Users/adrian/research/akureyri/results/deseq2_pipeline/annotation.tsv'

#
# 3. read information
# 

# 3.1. counts data
df = read.table(counts_file, sep='\t', header=TRUE, row.names=1)
dim(df)
drop <- c("Average.Normal","Max")
df = df[, !(names(df) %in% drop)]
dim(df)

epsilon = 2**-4.647587
y = (2**df) - epsilon
counts = round(y, digits=0)

# 3.2. read gene length
ann = read.table(annotation_file, 
                  sep='\t', 
                  header=TRUE, 
                  row.names=1,
                  quote="") # R surprises
print(dim(ann))

# 
# 4. analysis
#
gene_lengths = c()
found_symbols = c()
symbols = rownames(counts)

for (symbol in symbols[1:10]){
  gene_length = ann[ann$external_gene_name == symbol, ]$geneLength[1]
  print(c(symbol, gene_length))
  
  gene_lengths = c(gene_lengths, gene_length)
  found_symbols = c(found_symbols, symbol)
} 
gene_lengths

convertCounts(counts, unit='TPM')


