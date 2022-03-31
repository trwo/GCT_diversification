
#genes to remove to leave ones of more biological interest to compare
source('/lustre/scratch119/casm/team294rr/to3/testes/tumour/RNA_analysis/geneSets.R')
prot_cod = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/RNA_analysis/protein_coding_genes.txt', header = F, sep='\t', stringsAsFactors = F)[,1]

#read in data
tpm = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nat_comms_2022_R0/00_SUPPLEMENTARY_DATA/germ_cell_lcm_tpm_counts_filtered_samples_20210812.txt', header = T, sep = '\t', stringsAsFactors = F)
log2tpm <- log2(as.matrix(tpm) + 1)
log2tpm = log2tpm[(row.names(log2tpm) %in% prot_cod) & !(row.names(log2tpm) %in% unique(c(hgGenes, hkGenes, igGenes, cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))),] #remove genes whose expression will be related to normal contamination or cellular functions not related to phenotype
metadata = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nat_comms_2022_R0/00_SUPPLEMENTARY_DATA/mRNA_summary_metadata.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
row.names(metadata) = metadata$Sample
inv_mixed_rna = metadata[metadata$Diagnosis == 'Post-pubertal mixed germ cell tumour' & metadata$Normal.Tumour == 'Tumour' & metadata$Updated.Description != 'GCNIS',]
inv_mixed_rna = inv_mixed_rna[inv_mixed_rna$Case.ID %in% c('PD43296', 'PD43298', 'PD45543', 'PD45544', 'PD46269'),] #make direct contrast in each tumour to their transcriptome

#between histologies
interhisto_t_list = list()
for(i in unique(inv_mixed_rna$Case.ID)){ #per patient
  assessed_histo = c() #record of histologies already analysed
  t_prox_score_vec = c() #vector per patient comparison scores
  for(j in unique(inv_mixed_rna[inv_mixed_rna$Case.ID == i,]$Updated.Description)){
    assessed_histo = c(assessed_histo, j) #add histology to list to not compare against
    if(length(assessed_histo) != length(unique(inv_mixed_rna[inv_mixed_rna$Case.ID == i,]$Updated.Description))){ #stop when you are at the last histology as you have already compared it to all of the others
      for(k in unique(inv_mixed_rna[inv_mixed_rna$Case.ID == i & inv_mixed_rna$Updated.Description == j,]$Sample)){ #samples per histology per tumour
        choice_sample = log2tpm[, k] #vector from sample in histology of interest
        for(l in unique(inv_mixed_rna[inv_mixed_rna$Case.ID == i & !(inv_mixed_rna$Updated.Description %in% assessed_histo),]$Sample)){
          comparator_sample = log2tpm[, l] #vector from a sample in another histology
          t_prox = cor(choice_sample, comparator_sample, method = 'pearson')
          t_prox_score_vec = c(t_prox_score_vec, t_prox)
        }
      }
    }
    
  }
  interhisto_t_list = c(interhisto_t_list, list(t_prox_score_vec))
}

names(interhisto_t_list) = unique(inv_mixed_rna$Case.ID)

sapply(interhisto_t_list, median)

#within histologies
intrahisto_t_list = list()
for(i in unique(inv_mixed_rna$Case.ID)){ #per patient
  t_prox_score_vec = c() #vector per patient comparison scores
  for(j in unique(inv_mixed_rna[inv_mixed_rna$Case.ID == i,]$Updated.Description)){ #per histology
    if(nrow(inv_mixed_rna[inv_mixed_rna$Case.ID == i & inv_mixed_rna$Updated.Description == j,]) > 1){ #if histology has more than one microbiopsy
      assessed_samples = c() #record of samples already analysed
      for(k in unique(inv_mixed_rna[inv_mixed_rna$Case.ID == i & inv_mixed_rna$Updated.Description == j,]$Sample)){ #per sample
        assessed_samples = c(assessed_samples, k) #add to vector analysed
        choice_sample = log2tpm[, k]
        
        if(k != unique(inv_mixed_rna[inv_mixed_rna$Case.ID == i & inv_mixed_rna$Updated.Description == j,]$Sample)[length(unique(inv_mixed_rna[inv_mixed_rna$Case.ID == i & inv_mixed_rna$Updated.Description == j,]$Sample))]){ #so the last sample isn't compared against itself
          for(l in unique(inv_mixed_rna[inv_mixed_rna$Case.ID == i & inv_mixed_rna$Updated.Description == j & !(inv_mixed_rna$Sample %in% assessed_samples),]$Sample)){
            comparator_sample = log2tpm[, l]
            t_prox = cor(choice_sample, comparator_sample, method = 'pearson')
            t_prox_score_vec = c(t_prox_score_vec, t_prox)
          }
        }
      }
    }
  }
  intrahisto_t_list = c(intrahisto_t_list, list(t_prox_score_vec))
}

names(intrahisto_t_list) = unique(inv_mixed_rna$Case.ID)
sapply(intrahisto_t_list, median)

#control - against normal testis
norm_control = metadata[metadata$Normal.Tumour == 'Normal',]

norm_tumour_t_list = list()
for(i in unique(inv_mixed_rna$Case.ID)){ #per patient
  t_prox_score_vec = c() #vector per patient comparison scores
  for(j in inv_mixed_rna[inv_mixed_rna$Case.ID == i,]$Sample){
    choice_sample = log2tpm[, j]
    for(k in norm_control$Sample){
      comparator_sample = log2tpm[, k]
      t_prox = cor(choice_sample, comparator_sample, method = 'pearson')
      #t_prox = cosine(choice_sample, comparator_sample)[1,1]
      #t_prox = cor(choice_sample, comparator_sample, method = 'spearman')
      t_prox_score_vec = c(t_prox_score_vec, t_prox)
    }
  }
  norm_tumour_t_list = c(norm_tumour_t_list, list(t_prox_score_vec))
}

names(norm_tumour_t_list) = unique(inv_mixed_rna$Case.ID)

sapply(norm_tumour_t_list, median)

#add together
interhisto_t_df = data.frame(unlist(interhisto_t_list))
interhisto_t_df$patient = substr(row.names(interhisto_t_df), 0, 7)
interhisto_t_df$comparison = 'Inter-histology'
row.names(interhisto_t_df) = NULL
colnames(interhisto_t_df) = c('t_prox', 'patient', 'comparison')

intrahisto_t_df = data.frame(unlist(intrahisto_t_list))
intrahisto_t_df$patient = substr(row.names(intrahisto_t_df), 0, 7)
intrahisto_t_df$comparison = 'Intra-histology'
row.names(intrahisto_t_df) = NULL
colnames(intrahisto_t_df) = c('t_prox', 'patient', 'comparison')

norm_tumour_t_df = data.frame(unlist(norm_tumour_t_list))
norm_tumour_t_df$patient = substr(row.names(norm_tumour_t_df), 0, 7)
norm_tumour_t_df$comparison = 'Normal v Tumour'
row.names(norm_tumour_t_df) = NULL
colnames(norm_tumour_t_df) = c('t_prox', 'patient', 'comparison')

all_histo_t_df = rbind(interhisto_t_df, intrahisto_t_df, norm_tumour_t_df)

#assess significance using the permutation test
permutation_list = list()
set.seed(42)

#check how many permutations for label switching can be achieved per tumour
patient = c("PD43296", "PD43298", "PD45543", "PD45544", "PD46269") 

for(i in patient){
  df = all_histo_t_df[all_histo_t_df$patient == i,]
  max_perm = choose(nrow(df), nrow(df[df$comparison == 'Inter-histology',]))
  print(paste0(i, ' ', max_perm))
}

#well over 1000 in all cases
#all tumours
patient = c("PD45543", "PD43296", "PD43298", "PD45544", "PD46269")
nRand = 1000

for(i in patient){
  median_diff_vec = c()
  df = all_histo_t_df[all_histo_t_df$patient == i,]
  for(j in 1:nRand){
    df$label.switch = sample(df$comparison, length(df$comparison), replace = F)
    median_diff_vec = c(median_diff_vec, abs(median(df[df$label.switch == 'Inter-histology',]$t_prox) - median(df[df$label.switch == 'Intra-histology',]$t_prox)))
  }
  permutation_list = c(permutation_list, list(median_diff_vec))
}

names(permutation_list) = patient

#what are the difference between the inter and intra histology values that we observe?
obs_t_data = aggregate(t_prox ~ patient + comparison, all_histo_t_df, median)
patient = c("PD45543", "PD43296", "PD43298", "PD45544", "PD46269")

obs_t_med_diff = c()
for(i in patient){
  obs_t_med_diff = rbind(obs_t_med_diff, c(i, abs(obs_t_data[obs_t_data$patient == i & obs_t_data$comparison == "Intra-histology",]$t_prox - obs_t_data[obs_t_data$patient == i & obs_t_data$comparison == "Inter-histology",]$t_prox)))
}

obs_t_med_diff = data.frame(obs_t_med_diff)
names(obs_t_med_diff) = c('patient', 'observed_median_diff')
obs_t_med_diff$observed_median_diff = as.numeric(obs_t_med_diff$observed_median_diff)

#what is the probability this was observed by chance?
obs_t_med_diff$pval = NA
for(i in patient){
  obs_t_med_diff[obs_t_med_diff$patient == i,]$pval = sum(permutation_list[[i]] >= obs_t_med_diff[obs_t_med_diff$patient == i,]$observed_median_diff) / length(permutation_list[[i]]) #what proportion of label switched combinations generates a larger difference between the inter and intra-histology transcriptomic proximity?
}

# patient observed_median_diff  pval
# PD45543           0.08906052 0.015
# PD43296           0.05536096 0.000
# PD43298           0.12263423 0.000
# PD45544           0.17752213 0.000
# PD46269           0.12499585 0.000
