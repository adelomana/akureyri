#
# -1. packages installation
#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18")
# 
# # install sleuth
# BiocManager::install("rhdf5")
# devtools::install_github("pachterlab/sleuth")

# installing some added packages needed for biomart to function
# install("BiocFileCache", force=TRUE)

rm(list = ls())

library(biomaRt)
library(sleuth)
library(crayon) # so the messages are blue
library(tictoc)

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

# useful for debugging
# tempo = listAttributes(mart)
# View(tempo)

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

#
# 3. create a sleuth object
#
s2c = metadata
print(s2c)

# prepare object for sleuth
cat(blue('preparing sleuth object...'), fill=TRUE)
tic()
so = sleuth_prep(s2c, 
                 target_mapping=t2g, 
                 aggregation_column='ens_gene', 
                 read_bootstrap_tpm=TRUE, 
                 extra_bootstrap_summary=TRUE,
                 gene_mode=TRUE) # omit this and you'll get transcripts. not a threat.
toc()
# theads option does not impact elapsed time

#
# 4. store estimated counts per gene. Sleuth recommends "scaled_reads_per_base"
#
cat(blue('storing'), fill=TRUE)

tpm_table = sleuth_to_matrix(so, 'obs_norm', 'scaled_reads_per_base')
write.table(tpm_table, 
          file.path(results_dir, 'sleuth_scaled_reads_bygene.tsv'),
          quote = FALSE,
          sep = '\t')

tpm_table = sleuth_to_matrix(so, 'obs_norm', 'tpm')
write.table(tpm_table, 
            file.path(results_dir, 'sleuth_tpm_bygene.tsv'),
            quote = FALSE,
            sep = '\t')
