#Pairwise similarity scores

###########
##GENOMIC##
###########

library(ggplot2)
library(ggpubr)
library(ggbeeswarm)

dna_cuts = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/supplementary_data/GCT_DNA_supplementary_data.txt', header = T, sep = '\t', stringsAsFactors = F)

mixed_feature_cuts = dna_cuts[dna_cuts$Histology != 'GCNIS' & dna_cuts$Diagnosis == 'Post-pubertal mixed germ cell tumour',]

pluri_mixed_feature_cuts = mixed_feature_cuts[!mixed_feature_cuts$Case_ID %in% c('PD45545', 'PD46966', 'PD46969'),] #only one histology in these tumours represented or multiple but spread a thinly across multiple biopsies unevenly which would distort a genetic proximity score

#between histologies
interhisto_g_list = list()
for(i in unique(pluri_mixed_feature_cuts$Case_ID)){ #per patient
  assessed_histo = c() #record of histologies already analysed
  g_prox_score_vec = c() #vector of per patient comparison scores
  for(j in unique(pluri_mixed_feature_cuts[pluri_mixed_feature_cuts$Case_ID == i,]$Histology)){ #per histology
    assessed_histo = c(assessed_histo, j) #add histology to list not to compare against
    if(j != unique(pluri_mixed_feature_cuts[pluri_mixed_feature_cuts$Case_ID == i,]$Histology)[length(unique(pluri_mixed_feature_cuts[pluri_mixed_feature_cuts$Case_ID == i,]$Histology))]){ #no groups left for the final histology to be compared against
      for(k in unique(pluri_mixed_feature_cuts[pluri_mixed_feature_cuts$Case_ID == i & pluri_mixed_feature_cuts$Histology == j,]$Sample)){ #each patient sample per histology
        choice_sample = read.table(paste0('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/final_variant_calls/SNV/', k, '_post_swater_final_snvs.txt'), header = T, sep = '\t', stringsAsFactors = F)
        choice_sample = choice_sample[choice_sample$VAF >= 0.1,] #remove low level contamination from other histologies
        for(l in unique(pluri_mixed_feature_cuts[pluri_mixed_feature_cuts$Case_ID == i & !(pluri_mixed_feature_cuts$Histology %in% assessed_histo),]$Sample)){
          comparator_sample = read.table(paste0('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/final_variant_calls/SNV/', l, '_post_swater_final_snvs.txt'), header = T, sep = '\t', stringsAsFactors = F)
          comparator_sample = comparator_sample[comparator_sample$VAF >= 0.1,] #remove low level contamination from other histologies
          g_prox = length(intersect(choice_sample$ID, comparator_sample$ID)) / ((nrow(choice_sample) + nrow(comparator_sample)) / 2) #average sharedness of high VAF mutations
          g_prox_score_vec = c(g_prox_score_vec, g_prox) #add pairwise comparison to patient vector
        }
      }
    }
  }
  interhisto_g_list = c(interhisto_g_list, list(g_prox_score_vec))
}

names(interhisto_g_list) = unique(pluri_mixed_feature_cuts$Case_ID)
sapply(interhisto_g_list, median)

#within histologies
intrahisto_g_list = list()
for(i in unique(pluri_mixed_feature_cuts$Case_ID)){ #per patient
  g_prox_score_vec = c() #vector of per patient comparison scores
  for(j in unique(pluri_mixed_feature_cuts[pluri_mixed_feature_cuts$Case_ID == i,]$Histology)){ #per histology
    if(nrow(pluri_mixed_feature_cuts[pluri_mixed_feature_cuts$Case_ID == i & pluri_mixed_feature_cuts$Histology == j,]) > 1){ #if histology has more than one microbiopsy
      assessed_samples = c() #record of samples already analysed
      for(k in unique(pluri_mixed_feature_cuts[pluri_mixed_feature_cuts$Case_ID == i & pluri_mixed_feature_cuts$Histology == j,]$Sample)){ #per sample
        assessed_samples = c(assessed_samples, k) #add to vector of analysed samples
        choice_sample = read.table(paste0('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/final_variant_calls/SNV/', k, '_post_swater_final_snvs.txt'), header = T, sep = '\t', stringsAsFactors = F)
        choice_sample = choice_sample[choice_sample$VAF >= 0.1,] #remove low level contamination from other histologies
        if(k != unique(pluri_mixed_feature_cuts[pluri_mixed_feature_cuts$Case_ID == i & pluri_mixed_feature_cuts$Histology == j,]$Sample)[length(unique(pluri_mixed_feature_cuts[pluri_mixed_feature_cuts$Case_ID == i & pluri_mixed_feature_cuts$Histology == j,]$Sample))]){ #as the last sample cannot be compared against itself
          for(l in unique(pluri_mixed_feature_cuts[pluri_mixed_feature_cuts$Case_ID == i & pluri_mixed_feature_cuts$Histology == j & !(pluri_mixed_feature_cuts$Sample %in% assessed_samples),]$Sample)){
            comparator_sample = read.table(paste0('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/final_variant_calls/SNV/', l, '_post_swater_final_snvs.txt'), header = T, sep = '\t', stringsAsFactors = F)
            comparator_sample = comparator_sample[comparator_sample$VAF >= 0.1,] #remove low level contamination from other histologies
            g_prox = length(intersect(choice_sample$ID, comparator_sample$ID)) / ((nrow(choice_sample) + nrow(comparator_sample)) / 2) #average sharedness of high VAF mutations
            g_prox_score_vec = c(g_prox_score_vec, g_prox) #add pairwise comparison to patient vector
          }
        }
        
      }
    }
  }
  intrahisto_g_list = c(intrahisto_g_list, list(g_prox_score_vec))
}

names(intrahisto_g_list) = unique(pluri_mixed_feature_cuts$Case_ID)
sapply(intrahisto_g_list, median)

#combine results
interhisto_g_df = data.frame(unlist(interhisto_g_list))
interhisto_g_df$patient = substr(row.names(interhisto_g_df), 0, 7)
interhisto_g_df$comparison = 'Inter-histology'
row.names(interhisto_g_df) = NULL
colnames(interhisto_g_df) = c('gen_prox', 'patient', 'comparison')

intrahisto_g_df = data.frame(unlist(intrahisto_g_list))
intrahisto_g_df$patient = substr(row.names(intrahisto_g_df), 0, 7)
intrahisto_g_df$comparison = 'Intra-histology'
row.names(intrahisto_g_df) = NULL
colnames(intrahisto_g_df) = c('gen_prox', 'patient', 'comparison')

all_histo_g_df = rbind(interhisto_g_df, intrahisto_g_df)
all_histo_g_df$comparison = factor(all_histo_g_df$comparison, levels = c('Intra-histology', 'Inter-histology'))
all_histo_g_df$patient = as.factor(all_histo_g_df$patient)

#are these differences significant? no
PD43296_g = wilcox.test(gen_prox ~ comparison, all_histo_g_df[all_histo_g_df$patient == 'PD43296',]) #W = 24522, p-value = 0.1284
PD43298_g = wilcox.test(gen_prox ~ comparison, all_histo_g_df[all_histo_g_df$patient == 'PD43298',]) #W = 12651, p-value = 0.8893
PD45543_g = wilcox.test(gen_prox ~ comparison, all_histo_g_df[all_histo_g_df$patient == 'PD45543',]) #W = 4, p-value = 1
PD45544_g = wilcox.test(gen_prox ~ comparison, all_histo_g_df[all_histo_g_df$patient == 'PD45544',]) #W = 151, p-value = 0.4753
PD46269_g = wilcox.test(gen_prox ~ comparison, all_histo_g_df[all_histo_g_df$patient == 'PD46269',]) #W = 145, p-value = 0.6483

#adjust for multiple testing
p_g_vec = c(PD43296_g$p.value, PD43298_g$p.value, PD45543_g$p.value, PD45544_g$p.value, PD46269_g$p.value) 
p.adj_g_vec = p.adjust(p_g_vec, method="BH")
names(p.adj_g_vec) = unique(pluri_mixed_feature_cuts$Case_ID)

##################
##TRANSCRIPTOMIC##
##################

library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)

#genes to remove to leave ones of more biological interest to compare
source('/lustre/scratch119/casm/team294rr/to3/testes/tumour/RNA_analysis/geneSets.R')
prot_cod = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/RNA_analysis/protein_coding_genes.txt', header = F, sep='\t', stringsAsFactors = F)[,1]

#read in data
tpm = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/pre_processing/tmm_norm_gct_tpm_20210304.txt', header = T, sep = '\t', stringsAsFactors = F)
log2tpm <- log2(as.matrix(tpm) + 1)
log2tpm = log2tpm[(row.names(log2tpm) %in% prot_cod) & !(row.names(log2tpm) %in% unique(c(hgGenes, hkGenes, igGenes, riboGenes, riboRNAGenes, cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))),]
metadata = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/supplementary_data/GCT_RNA_metadata.txt', header = T, sep = '\t', stringsAsFactors = F)
row.names(metadata) = metadata$Sample_ID
inv_mixed_rna = metadata[metadata$Diagnosis == 'Post-pubertal mixed germ cell tumour' & metadata$Normal.Tumour == 'Tumour' & metadata$Updated.Description != 'GCNIS',]
inv_mixed_rna = inv_mixed_rna[inv_mixed_rna$Case_ID %in% c('PD43296', 'PD43298', 'PD45543', 'PD45544', 'PD46269'),] #make direct contrast in each tumour to their transcriptome

#between histologies
interhisto_t_list = list()
for(i in unique(inv_mixed_rna$Case_ID)){ #per patient
  assessed_histo = c() #record of histologies already analysed
  t_prox_score_vec = c() #vector per patient comparison scores
  for(j in unique(inv_mixed_rna[inv_mixed_rna$Case_ID == i,]$Updated.Description)){
    assessed_histo = c(assessed_histo, j) #add histology to list to not compare against
    if(length(assessed_histo) != length(unique(inv_mixed_rna[inv_mixed_rna$Case_ID == i,]$Updated.Description))){ #stop when you are at the last histology as you have already compared it to all of the others
      for(k in unique(inv_mixed_rna[inv_mixed_rna$Case_ID == i & inv_mixed_rna$Updated.Description == j,]$Sample_ID)){ #samples per histology per tumour
        choice_sample = log2tpm[, k] #vector from sample in histology of interest
        for(l in unique(inv_mixed_rna[inv_mixed_rna$Case_ID == i & !(inv_mixed_rna$Updated.Description %in% assessed_histo),]$Sample_ID)){
          comparator_sample = log2tpm[, l] #vector from a sample in another histology
          t_prox = cor(choice_sample, comparator_sample, method = 'pearson')
          t_prox_score_vec = c(t_prox_score_vec, t_prox)
        }
      }
    }
    
  }
  interhisto_t_list = c(interhisto_t_list, list(t_prox_score_vec))
}

names(interhisto_t_list) = unique(inv_mixed_rna$Case_ID)

sapply(interhisto_t_list, median)

#within histologies
intrahisto_t_list = list()
for(i in unique(inv_mixed_rna$Case_ID)){ #per patient
  t_prox_score_vec = c() #vector per patient comparison scores
  for(j in unique(inv_mixed_rna[inv_mixed_rna$Case_ID == i,]$Updated.Description)){ #per histology
    if(nrow(inv_mixed_rna[inv_mixed_rna$Case_ID == i & inv_mixed_rna$Updated.Description == j,]) > 1){ #if histology has more than one microbiopsy
      assessed_samples = c() #record of samples already analysed
      for(k in unique(inv_mixed_rna[inv_mixed_rna$Case_ID == i & inv_mixed_rna$Updated.Description == j,]$Sample)){ #per sample
        assessed_samples = c(assessed_samples, k) #add to vector analysed
        choice_sample = log2tpm[, k]
        
        if(k != unique(inv_mixed_rna[inv_mixed_rna$Case_ID == i & inv_mixed_rna$Updated.Description == j,]$Sample)[length(unique(inv_mixed_rna[inv_mixed_rna$Case_ID == i & inv_mixed_rna$Updated.Description == j,]$Sample))]){ #so the last sample isn't compared against itself
          for(l in unique(inv_mixed_rna[inv_mixed_rna$Case_ID == i & inv_mixed_rna$Updated.Description == j & !(inv_mixed_rna$Sample %in% assessed_samples),]$Sample)){
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

names(intrahisto_t_list) = unique(inv_mixed_rna$Case_ID)
sapply(intrahisto_t_list, median)

#control - against normal testis
norm_control = metadata[metadata$Normal.Tumour == 'Normal',]

norm_tumour_t_list = list()
for(i in unique(inv_mixed_rna$Case_ID)){ #per patient
  t_prox_score_vec = c() #vector per patient comparison scores
  for(j in inv_mixed_rna[inv_mixed_rna$Case_ID == i,]$Sample_ID){
    choice_sample = log2tpm[, j]
    for(k in norm_control$Sample_ID){
      comparator_sample = log2tpm[, k]
      t_prox = cor(choice_sample, comparator_sample, method = 'pearson')
      #t_prox = cosine(choice_sample, comparator_sample)[1,1]
      #t_prox = cor(choice_sample, comparator_sample, method = 'spearman')
      t_prox_score_vec = c(t_prox_score_vec, t_prox)
    }
  }
  norm_tumour_t_list = c(norm_tumour_t_list, list(t_prox_score_vec))
}

names(norm_tumour_t_list) = unique(inv_mixed_rna$Case_ID)

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

ggplot(all_histo_t_df) + 
  geom_quasirandom(mapping = aes(x = comparison, y = t_prox, col = comparison), width = 0.3) +
  geom_pointrange(mapping = aes(x = comparison, y = t_prox), stat = "summary", 
                  shape=19, 
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill="black") +
  facet_grid(. ~ patient, scales = 'free_x', space = 'free_x', switch = 'x') + theme_pubr() + theme(axis.text.x = element_blank(), panel.spacing = unit(0, 'mm'), strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 11)) + ylab('Transcriptomic similarity') + xlab ('') + ylim(0, 1)

PD43296_t = wilcox.test(t_prox ~ comparison, all_histo_t_df[all_histo_t_df$patient == 'PD43296' & all_histo_t_df$comparison %in% c('Inter-histology', 'Intra-histology'),]) #W = 117762, p-value < 2.2e-16
PD43298_t = wilcox.test(t_prox ~ comparison, all_histo_t_df[all_histo_t_df$patient == 'PD43298' & all_histo_t_df$comparison %in% c('Inter-histology', 'Intra-histology'),]) #W = 51363, p-value < 2.2e-16
PD45543_t = wilcox.test(t_prox ~ comparison, all_histo_t_df[all_histo_t_df$patient == 'PD45543' & all_histo_t_df$comparison %in% c('Inter-histology', 'Intra-histology'),]) #W = 7, p-value = 0.01399
PD45544_t = wilcox.test(t_prox ~ comparison, all_histo_t_df[all_histo_t_df$patient == 'PD45544' & all_histo_t_df$comparison %in% c('Inter-histology', 'Intra-histology'),]) #W = 3492, p-value < 2.2e-16
PD46269_t = wilcox.test(t_prox ~ comparison, all_histo_t_df[all_histo_t_df$patient == 'PD46269' & all_histo_t_df$comparison %in% c('Inter-histology', 'Intra-histology'),]) #W = 15259, p-value < 2.2e-16

p_t_vec = c(PD43296_t$p.value, PD43298_t$p.value, PD45543_t$p.value, PD45544_t$p.value, PD46269_t$p.value) 
p.adj_t_vec = p.adjust(p_t_vec, method="BH")
names(p.adj_t_vec) = unique(inv_mixed_rna$Case_ID)

all_histo_t_df$comparison = factor(all_histo_t_df$comparison, levels = c('Intra-histology', 'Inter-histology', 'Normal v Tumour'))

########
##PLOT##
########

library(cowplot)

p1 = ggplot(all_histo_g_df) + 
  geom_quasirandom(mapping = aes(x = gen_prox, y = comparison, col = comparison), width = 0.3, groupOnX = F) +
  geom_pointrange(mapping = aes(x = gen_prox, y = comparison),
                  stat = "summary", shape=19, 
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill = "black") +
  facet_grid(patient ~ ., scales = 'free_y', space = 'free_y') + theme_pubr() + theme(axis.text.y = element_blank(), panel.spacing = unit(0, 'mm'), strip.background = element_blank(), strip.placement = "none", strip.text = element_text(size = 11)) + xlab('Genetic similarity') + ylab ('') + scale_y_discrete(position = "right")  + scale_x_continuous(expand = c(0,0), limits = c(0,1)) + scale_color_manual(values = c("#F2AD00", "#5BBCD6"), name = 'Comparison')

p2 = ggplot(all_histo_t_df) + 
  geom_quasirandom(mapping = aes(x = t_prox, y = comparison, col = comparison), width = 0.4, groupOnX = F) +
  geom_pointrange(mapping = aes(x = t_prox, y = comparison),
                  stat = "summary", shape=19, 
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill = "black") +
  facet_grid(patient ~ ., scales = 'free_y', space = 'free_y', switch = 'y') + theme_pubr() + theme(axis.text.y = element_blank(), panel.spacing = unit(0, 'mm'), strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 11)) + xlab('Transcriptomic similarity') + ylab ('') + scale_y_discrete(position = "left")  + scale_x_reverse(expand = c(0,0), limits = c(1,0)) + scale_color_manual(values = c("#F2AD00", "#5BBCD6", '#999999'), name = 'Comparison')

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/gct_genetic_proximity_pearson_rna_20210305.pdf', height = 6, width = 8, useDingbats = F)
plot_grid(p1, p2, align = 'h', ncol = 2)
dev.off()