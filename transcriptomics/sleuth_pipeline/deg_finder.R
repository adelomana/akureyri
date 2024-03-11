rm(list = ls())

library(biomaRt)
library(sleuth)

#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/akureyri/results/kallisto.100"
results_dir = '/Users/adrian/research/akureyri/results/sleuth_pipeline'

#
# 1. generate gene to transcript mapping
#

# annotation defined from sleuth walk through, https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="hsapiens_gene_ensembl",
                        host = 'https://www.ensembl.org',
                        verbose = TRUE)
t2g = biomaRt::getBM(attributes=c("ensembl_transcript_id", 
                                  "ensembl_gene_id", 
                                  "external_gene_name"), mart=mart)
t2g = dplyr::rename(t2g, target_id=ensembl_transcript_id, ens_gene=ensembl_gene_id, ext_gene=external_gene_name)

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
# View(metadata)

#
# 3. contrasts
#

# using LRT instead of Wald because Wald gives three times more.
# Using pval_aggregate = TRUE bc otherwise no DEGs. Get FC from est counts file

#
# 3.1. contrast effect of time for 2D
#
rule = metadata$culture == '2D'
s2c = metadata[rule, ]
dim(s2c)
View(s2c)
# prepare contrast
so = sleuth_prep(s2c,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 read_bootstrap_tpm = TRUE)
# contrast 
so = sleuth_fit(so, ~time, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
sleuth_table = sleuth_results(so, 
                               'reduced:full', 
                               'lrt',
                               show_all = FALSE,
                               pval_aggregate = TRUE)
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05)
dim(sleuth_significant)
plot_pca(so, color_by = 'time', title = 'effect time for 2D')
write.table(sleuth_significant, 
            file = paste(results_dir, '/effect_time_for_2D.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)

#
# 3.2. contrast effect of time for 3D
#
rule = metadata$culture == '3D'
s2c = metadata[rule, ]
dim(s2c)
View(s2c)
# prepare contrast
so = sleuth_prep(s2c,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 read_bootstrap_tpm = TRUE)
# contrast 
so = sleuth_fit(so, ~time, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
sleuth_table = sleuth_results(so, 
                              'reduced:full', 
                              'lrt',
                              show_all = FALSE,
                              pval_aggregate = TRUE)
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05)
dim(sleuth_significant)
plot_pca(so, color_by = 'time',  title = 'effect time for 3D')
write.table(sleuth_significant, 
            file = paste(results_dir, '/effect_time_for_3D.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)

#
# 3.3. contrast effect of culture at 2 days
#
rule = metadata$time == 'two'
s2c = metadata[rule, ]
dim(s2c)
View(s2c)
# prepare contrast
so = sleuth_prep(s2c,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 read_bootstrap_tpm = TRUE)
# contrast 
so = sleuth_fit(so, ~culture, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
sleuth_table = sleuth_results(so, 
                              'reduced:full', 
                              'lrt',
                              show_all = FALSE,
                              pval_aggregate = TRUE)
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05)
dim(sleuth_significant)
plot_pca(so, color_by = 'culture',  title = 'effect culture at 2 days')
write.table(sleuth_significant, 
            file = paste(results_dir, '/effect_culture_at_2days.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)

#
# 3.4. contrast effect of culture at 14 days
#
rule = metadata$time == 'fourteen'
s2c = metadata[rule, ]
dim(s2c)
View(s2c)
# prepare contrast
so = sleuth_prep(s2c,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 read_bootstrap_tpm = TRUE)
# contrast 
so = sleuth_fit(so, ~culture, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
sleuth_table = sleuth_results(so, 
                              'reduced:full', 
                              'lrt',
                              show_all = FALSE,
                              pval_aggregate = TRUE)
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05)
dim(sleuth_significant)
plot_pca(so, color_by = 'culture', title = 'effect culture at 14 days')
write.table(sleuth_significant, 
            file = paste(results_dir, '/effect_culture_at_14days.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)

#
# 3.5 interaction
#
s2c = metadata
dim(s2c)
View(s2c)
# prepare contrast
so = sleuth_prep(s2c,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 read_bootstrap_tpm = TRUE)
# contrast 
so = sleuth_fit(so, ~time+culture+time:culture, 'full')
so = sleuth_fit(so, ~time+culture, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
sleuth_table = sleuth_results(so, 
                              'reduced:full', 
                              'lrt',
                              show_all = FALSE,
                              pval_aggregate = TRUE)
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05)
dim(sleuth_significant)
write.table(sleuth_significant, 
            file = paste(results_dir, '/interaction.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)





