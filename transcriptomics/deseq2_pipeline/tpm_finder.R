rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(biomaRt)

#
# 1. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/akureyri/results/kallisto.100"
results_dir = '/Users/adrian/research/akureyri/results/deseq2_pipeline'

#
# 1. generate gene to transcript mapping
#

mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="hsapiens_gene_ensembl",
                        host = 'https://www.ensembl.org',
                        verbose = TRUE)
working_attributes = c('ensembl_transcript_id', 
                      'ensembl_gene_id', 
                      'external_gene_name',
                      'gene_biotype',
                      'description')
t2g = biomaRt::getBM(attributes=working_attributes, 
                     mart=mart,
                     verbose=TRUE)
dim(t2g)
View(t2g)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
paths = file.path(dirnames, 'abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[8]))
times = c('two', 'two', 'two', 'two', 'fourteen', 'fourteen',
          'fourteen', 'fourteen', 'two', 'two', 'fourteen', 'fourteen')
sample = c('test01', 'test09', 'test02', 'test10', 'test03', 'test11',
           'test04', 'test12', 'test05', 'test06', 'test07', 'test08')
cultures = c('2D', '2D', '3D', '3D', '2D', '2D',
             '3D', '3D', '2D', '3D', '2D', '3D')

metadata = data.frame(sample)
metadata$label = labels
metadata$time = times
metadata$culture = cultures
metadata$path = paths
metadata <- metadata[order(metadata$sample), ]
View(metadata)

#
# 3. read files
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

#
# 4. find abundance
#
tpm = txi$abundance
colnames(tpm) = metadata$sample
dim(tpm)
View(tpm)

#
# 5. store
#
store = paste(results_dir, '/DESeq2_TPM_values.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)

store = paste(results_dir, '/annotation.tsv', sep='')
write.table(t2g, file=store, quote=FALSE, sep='\t', col.names=NA)
