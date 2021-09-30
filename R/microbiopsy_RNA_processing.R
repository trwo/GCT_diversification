
library(ggplot2)
library(ggpubr)

raw_data = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/original_nature_submission/analyses/RNA/filtering/germ_cell_lcm_rna_preQC.txt', header = T, sep = '\t', stringsAsFactors = F) #just normal testis and GCT raw count data, 500 microbiopsies

#remove transcriptomes that were:
#- identified as contaminated during microdissection
#- non-neoplastic tissues that are not healthy seminiferous tubules
#- uncertainty as to the histological categorisation of the microbiopsy

raw_data = raw_data[, !colnames(raw_data) %in% c("PD43298a_Z1_2_SL02_04_SEC02_04_F2", 'PR43299c_lo0009', 'PR43299c_lo0010', 'PR43299c_lo0011', 'PR43299c_lo0012', "PD43296a_SL02_SEC02_G1", "PD43296a_Z1_2_SL02_04_SEC02_04_E1A", "PD43298a_SL02_SEC02_B3", "PD43298a_SL02_SEC02_C3", "PD43298a_Z1_2_SL02_04_SEC02_04_G2", "PR43299c_lo0002", "PR43299c_lo0006", "PR42036a_SL0D_SEC02_A5", "PR42036a_SL0D_SEC02_G3", "PR42036a_SL0D_SEC02_H3", "PR43298a_SL04_SEC04_D5",             "PR43298a_SL04_SEC04_E5" ,"PR43298a_SL04_SEC04_G9", "PR45543a_SL02_SEC03_B4", "PR45543a_SL02_SEC03_D2", "PR45543a_SL02_SEC03_F2", "PR46269c_SL02_SEC03_E12", "PR46269c_SL02_SEC03_F12", "PR46270c_SL02_SEC03_B5", "PR46270c_SL02_SEC03_C5", "PR46270c_SL02_SEC03_D5", "PR46270c_SL02_SEC03_E3", "PR46270c_SL02_SEC03_F3",             "PR46270c_SL02_SEC03_G3", "PR46270c_SL02_SEC03_H3", "PR46967a_SL02_SEC03_D8", "PR46967a_SL02_SEC03_E4", "PR46967a_SL02_SEC03_G6", "PR46968d_SL02_SEC02_03_H5", "PR46969c_SL02_SEC03_C11", "PR46969c_SL02_SEC03_D11", "PR46969c_SL02_SEC03_E11", "PR46969c_SL02_SEC03_F11", "PR46969c_SL02_SEC03_G11", "PR46969c_SL02_SEC03_H11")] #460 samples left

#how many reads, across all 55502 mapped features do we have?
summary(colSums(raw_data[, 7:ncol(raw_data)])) #highly variable library sizes
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   682   86612  246980  502107  651605 3863511 

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/figures/all_samples_read_count_QC.pdf', height = 6, width = 6, useDingbats = F)
ggplot() +
  geom_histogram(mapping = aes(colSums(raw_data[, 7:ncol(raw_data)])), bins = 100, fill = '#999999') + theme_pubr() + coord_cartesian(expand = F) + xlab('Number of reads across all mapped features') + ylab('Number of microbiopsies')
#spike in samples with under 1000 features expressed
dev.off()

#if we choose a threshold of 5 reads as the minimum to consider a gene truly expressed, how many are expressed per microbiopsy?
summary(colSums(raw_data[, 7:ncol(raw_data)] >= 5))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0    4150    8916    9077   13654   22865 

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/figures/all_samples_feature_expression_QC.pdf', height = 6, width = 6, useDingbats = F)
ggplot() +
  geom_histogram(mapping = aes(colSums(raw_data[, 7:ncol(raw_data)] >= 5)), bins = 100, fill = '#999999') + theme_pubr() + coord_cartesian(expand = F) + geom_vline(xintercept = 1000, lty = 2, col = 'red') + xlab('Features expressed per microbiopsy to a depth of at least 5 reads') + ylim(0, 15) + ylab('Number of microbiopsies')
#spike in samples with under 1000 features expressed
dev.off()

sum(colSums(raw_data[, 7:ncol(raw_data)] >= 5) >= 1000) #416

filtered_data = raw_data[, c(rep(TRUE, 6), colSums(raw_data[, 7:ncol(raw_data)] >= 5) >= 1000)]

write.table(filtered_data, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/supplementary_data/germ_cell_lcm_raw_counts_filtered_samples_20210812.txt', col.names = T, row.names = F, quote = F, sep = '\t')

#generate tpm values
x = filtered_data[, 7:ncol(filtered_data)]
row.names(x) = filtered_data$Geneid
x = x / filtered_data$Length #normalise for transcript length
tpm <- t( t(x) * 1e6 / colSums(x)) #normalise to read depth
write.table(tpm, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/supplementary_data/germ_cell_lcm_tpm_counts_filtered_samples_20210812.txt', col.names = T, row.names = T, quote = F, sep = '\t')

#add total counts to mrna data
mrna = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/supplementary_data/mRNA_summary_metadata.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '') #descriptive column contains apostrophes
raw_data = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/rna_counts/germ_cell_lcm_rna_preQC.txt', header = T, sep = '\t', stringsAsFactors = F) #raw count data
row.names(mrna) = mrna$Sample
mrna$raw_counts = colSums(raw_data[, 7:ncol(raw_data)])[row.names(mrna)]
summary(mrna$raw_counts)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 13822  125931  298080  553863  716367 3863511 

write.table(mrna, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/supplementary_data/mRNA_summary_metadata.txt', col.names = T, row.names = F, sep = '\t', quote = F)
