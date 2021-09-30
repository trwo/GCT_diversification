
setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/supplementary_data/')

patient = c("PD43296", "PD43298", "PD45544", "PD46269") #copy over relevant files
for(i in patient){
  system(paste0('cp /lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/multidimdpclust_noX_new_swater/01_Input/', i, '/md_out/', i, '_filtered_clusters_median_CCF.txt .'))
}

system('cp /lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/upto5_dpclust_noX_new_swater/03_Output/PD45543_DPoutput_3000iters_1000burnin_seed123/PD45543_3000iters_1000burnin_bestClusterInfo.txt .') #changed headers and name to match other files 27.9.21

patient = c("PD43296", "PD43298", "PD45543", "PD45544", "PD46269") #only mixed tumours with multiple histologies represented within each biopsy taken

ignore_clusters = list(PD43296 = c(15, 7, 8), PD43298 = c(5, 6, 8), PD45543 = c(NULL), PD45544 = c(3), PD46269 = c(2)) #low VAF/artefactual clusters that shouldn't be included in similarity metric

manifest = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/supplementary_data/WGS_summary_metadata.txt', header = T, sep = '\t', stringsAsFactors = F)

#between histologies
interhisto_g_list = list()
for(i in patient){ #per patient
  mut_clusters = read.table(paste0(i, '_filtered_clusters_median_CCF.txt'), header = T, sep = '\t', stringsAsFactors = F) #dataframe of mutations assigned per cluster by multiDPClust
  mut_clusters = mut_clusters[!mut_clusters$cluster.no %in% ignore_clusters[[i]],] # remove bad clusters
  assessed_histo = c() #record of histologies already analysed
  g_prox_score_vec = c() #vector of per patient comparison scores
  
  patient_meta = manifest[manifest$Case.ID == i & !is.na(manifest$Median.CCF),] #patient samples that are included in the phylogeny
  
  for(j in unique(patient_meta$Histology)){ #per histology
    
    assessed_histo = c(assessed_histo, j) #add histology to list not to compare against
    
    if(j != unique(patient_meta$Histology)[length(unique(patient_meta$Histology))]){ #no groups left for the final histology to be compared against
      
      for(k in patient_meta[patient_meta$Histology == j,]$Sample){ #each patient sample per histology
        choice_sample_clusters_vec = mut_clusters[mut_clusters[,paste0(k, "_subclonalFraction")] > 0.1,]$cluster.no #clusters found at a CCF of at least 10% in sample of interest
        
        for(l in patient_meta[patient_meta$Histology != j & !(patient_meta$Histology %in% assessed_histo),]$Sample){ #each sample from another histology
          comparator_sample_clusters_vec = mut_clusters[mut_clusters[,paste0(l, "_subclonalFraction")] > 0.1,]$cluster.no #clusters found at a CCF of at least 10% in sample of interest
          
          g_prox = sum(mut_clusters[mut_clusters$cluster.no %in% intersect(choice_sample_clusters_vec, comparator_sample_clusters_vec), ]$mut.no) / ((sum(mut_clusters[mut_clusters$cluster.no %in% choice_sample_clusters_vec, ]$mut.no) + sum(mut_clusters[mut_clusters$cluster.no %in% comparator_sample_clusters_vec, ]$mut.no))  / 2) #average sharedness of mutations between all clones in each sample
          g_prox_score_vec = c(g_prox_score_vec, g_prox) #add pairwise comparison to patient vector
        }
      }
    }
  }
  interhisto_g_list = c(interhisto_g_list, list(g_prox_score_vec))
}

names(interhisto_g_list) = patient
sapply(interhisto_g_list, median)

#within histologies
intrahisto_g_list = list()

for(i in patient){ #per patient
  
  mut_clusters = read.table(paste0(i, '_filtered_clusters_median_CCF.txt'), header = T, sep = '\t', stringsAsFactors = F) #dataframe of mutations assigned per cluster by multiDPClust
  mut_clusters = mut_clusters[!mut_clusters$cluster.no %in% ignore_clusters[[i]],] # remove bad clusters
  
  patient_meta = manifest[manifest$Case.ID == i & !is.na(manifest$Median.CCF),]
  
  g_prox_score_vec = c() #vector of per patient comparison scores
  
  for(j in unique(patient_meta$Histology)){ #per histology
    if(nrow(patient_meta[patient_meta$Histology == j,]) > 1){ #if histology has more than one microbiopsy
      
      assessed_samples = c() #record of samples already analysed
      
      for(k in patient_meta[patient_meta$Histology == j,]$Sample){ #per sample
        
        assessed_samples = c(assessed_samples, k) #add to vector of analysed samples
        choice_sample_clusters_vec = mut_clusters[mut_clusters[,paste0(k, "_subclonalFraction")] > 0.1,]$cluster.no
        
        if(k != patient_meta[patient_meta$Histology == j,]$Sample[length(patient_meta[patient_meta$Histology == j,]$Sample)]){ #as the last sample cannot be compared against itself
          for(l in patient_meta[patient_meta$Histology == j & !(patient_meta$Sample %in% assessed_samples),]$Sample){
            comparator_sample_clusters_vec = mut_clusters[mut_clusters[,paste0(l, "_subclonalFraction")] > 0.1,]$cluster.no
            
            g_prox = sum(mut_clusters[mut_clusters$cluster.no %in% intersect(choice_sample_clusters_vec, comparator_sample_clusters_vec), ]$mut.no) / ((sum(mut_clusters[mut_clusters$cluster.no %in% choice_sample_clusters_vec, ]$mut.no) + sum(mut_clusters[mut_clusters$cluster.no %in% comparator_sample_clusters_vec, ]$mut.no))  / 2) #average sharedness of mutations between all clones in each sample
            g_prox_score_vec = c(g_prox_score_vec, g_prox) #add pairwise comparison to patient vector
          }
        }
        
      }
    }
  }
  intrahisto_g_list = c(intrahisto_g_list, list(g_prox_score_vec))
}

names(intrahisto_g_list) = patient
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

#use permutation test to assess significance
permutation_list = list()
set.seed(42)

#check how many permutations for label switching can be achieved per tumour
patient = c("PD43296", "PD43298", "PD45543", "PD45544", "PD46269") 

for(i in patient){
  df = all_histo_g_df[all_histo_g_df$patient == i,]
  max_perm = choose(nrow(df), nrow(df[df$comparison == 'Inter-histology',]))
  print(paste0(i, ' ', max_perm))
}

#Only 15 possible combinations for PD45543 - here will we compare the other 14 combinations with the one we observe as 1000 permutations would simply resample the same sample space many times over

#PD45543 assessment
i = "PD45543"
df = all_histo_g_df[all_histo_g_df$patient == i,]
all_comb = data.frame(combn(nrow(df), nrow(df[df$comparison == 'Inter-histology',]), simplify = T)) #all possible combinations of label switching
which(df$comparison == 'Inter-histology') #X1 matches our observations
remaining_comb = all_comb[, c(2:15)]
median_diff_vec = c()

for(j in 1:ncol(remaining_comb)){
  df$label.switch = 'Intra-histology'
  df$label.switch[remaining_comb[,j]] = 'Inter-histology'
  median_diff_vec = c(median_diff_vec, abs(median(df[df$label.switch == 'Inter-histology',]$gen_prox) - median(df[df$label.switch == 'Intra-histology',]$gen_prox)))
}

permutation_list = c(permutation_list, list(median_diff_vec))

#remaining tumours
patient = c("PD43296", "PD43298", "PD45544", "PD46269")
nRand = 1000


for(i in patient){
  median_diff_vec = c()
  df = all_histo_g_df[all_histo_g_df$patient == i,]
  for(j in 1:nRand){
    df$label.switch = sample(df$comparison, length(df$comparison), replace = F)
    median_diff_vec = c(median_diff_vec, abs(median(df[df$label.switch == 'Inter-histology',]$gen_prox) - median(df[df$label.switch == 'Intra-histology',]$gen_prox)))
  }
  permutation_list = c(permutation_list, list(median_diff_vec))
}

names(permutation_list) = c('PD45543', patient)

#what are the difference between the inter and intra histology values that we observe?
obs_data = aggregate(gen_prox ~ patient + comparison, all_histo_g_df, median)
patient = c("PD45543", "PD43296", "PD43298", "PD45544", "PD46269")

obs_med_diff = c()
for(i in patient){
  obs_med_diff = rbind(obs_med_diff, c(i, abs(obs_data[obs_data$patient == i & obs_data$comparison == "Intra-histology",]$gen_prox - obs_data[obs_data$patient == i & obs_data$comparison == "Inter-histology",]$gen_prox)))
}

obs_med_diff = data.frame(obs_med_diff)
names(obs_med_diff) = c('patient', 'observed_median_diff')
obs_med_diff$observed_median_diff = as.numeric(obs_med_diff$observed_median_diff)

#what is the probability this was observed by chance?
obs_med_diff$pval = NA
for(i in patient){
  obs_med_diff[obs_med_diff$patient == i,]$pval = sum(permutation_list[[i]] >= obs_med_diff[obs_med_diff$patient == i,]$observed_median_diff) / length(permutation_list[[i]]) #what proportion of label switched combinations generates a larger difference between the inter and intra-histology genomic proximity?
}

#patient observed_median_diff      pval
#PD45543          0.036548715 0.4285714
#PD43296          0.000000000 1.0000000
#PD43298          0.008784081 0.3270000
#PD45544          0.006219076 0.2780000
#PD46269          0.062893469 0.0020000