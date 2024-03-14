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
library(BiocParallel)
library(crayon) 
library(ggplot2)

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
write.table(metadata,
            file = paste(results_dir, '/metadata.tsv', sep=''), 
            quote=FALSE, sep='\t')

#
# 4. contrasts
#
threshold = 10

#
# 4.1. contrast effect of time for 2D
#
rule = metadata$culture == '2D'
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~time) 
dds$time = relevel(dds$time, ref="two")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
anti_results = res[which(res$padj > 0.05 & abs(res$log2FoldChange) < 1), ]
cat(blue(paste('contrast effect of time for 2D:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(filtred_results, 
            file=paste(results_dir, '/effect_time_2D.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/effect_time_2D.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('time')) + ggtitle('effect time for 2D')
ggsave(file.path(results_dir, 'effect_time_2D.png'))

#
# 4.2. contrast effect of time for 3D
#
rule = metadata$culture == '3D'
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~time) 
dds$time = relevel(dds$time, ref="two")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
anti_results = res[which(res$padj > 0.05 & abs(res$log2FoldChange) < 1), ]
cat(blue(paste('contrast effect of time for 3D:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(filtred_results, 
            file=paste(results_dir, '/effect_time_3D.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, '/effect_time_3D.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('time')) + ggtitle('effect time for 3D')
ggsave(file.path(results_dir, 'effect_time_3D.png'))

#
# 4.3. contrast effect of culture on Day 2
#
rule = metadata$time == 'two'
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~culture) 
dds$culture = relevel(dds$culture, ref="2D")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
anti_results = res[which(res$padj > 0.05 & abs(res$log2FoldChange) < 1), ]
cat(blue(paste('contrast effect of culture on Day 2:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(filtred_results, 
            file=paste(results_dir, '/effect_culture_day2.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, '/effect_culture_day2.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('culture')) + ggtitle('effect culture on Day 2')
ggsave(file.path(results_dir, 'effect_culture_day2.png'))

#
# 4.4. contrast effect of culture on Day 14
#
rule = metadata$time == 'fourteen'
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~culture) 
dds$culture = relevel(dds$culture, ref="2D")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
anti_results = res[which(res$padj > 0.05 & abs(res$log2FoldChange) < 1), ]
cat(blue(paste('contrast effect of culture on Day 14:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(filtred_results, 
            file=paste(results_dir, '/effect_culture_day14.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, '/effect_culture_day14.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('culture')) + ggtitle('effect culture Day 14')
ggsave(file.path(results_dir, 'effect_culture_day14.png'))

#
# 4.5. interaction
#
working_metadata = metadata
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~time+culture+time:culture) 

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~time+culture)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
cat(blue(paste('interaction:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(filtred_results, 
            file=paste(results_dir, '/interaction.tsv', sep=''), quote=FALSE, sep='\t')

# 
# 4.6. plot everything
#
working_metadata = metadata
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~time+culture+time:culture) 

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

plotPCA(rlog(dds), ntop=500, intgroup=c('time', 'culture')) + ggtitle('all')
ggsave(file.path(results_dir, 'all.png'))
