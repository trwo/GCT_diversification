
#Site specific error rate estimation and filtering
##Adapted from Coorens et al 2019. PMID: 31806814

###SNV filtering

#-------------------------------------------------
# Libraries
#-------------------------------------------------

library("GenomicRanges")
library("deepSNV")
library("Rsamtools")

logbb = deepSNV:::logbb
dbetabinom = VGAM::dbetabinom

#-------------------------------------------------
# Functions
#-------------------------------------------------

estimateRho_gridml = function(x, mu) {
  # Estimate rho by MLE grid approach
  rhovec = 10^seq(-6,-0.5,by=0.05) # rho will be bounded within 1e-6 and 0.32
  mm = x[,2]
  cov = c(x[,1])
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=mm, size=cov, rho=rhoj, prob=mu, log=T)))
  rhovec[ll==max(ll)][1]
}

shearwater_probability=function(patient, save=NULL, path_prefix='', rho=10^-3){
  #Function to calculate probability of presence of mutation based on Shearwarer
  #'patient' is the name of the patient-specific subfolder
  #'path_prefix' is any prefix to the path necessary to find that subfolder
  #'rho' is the constant for the overdispersion parameter. If rho=NULL, calculate
  # it from the data (much slower)
  #'save' is path for output. If NULL, returns matrix
  # output of this function will be a matrix (muts by samples) of p values
  
  
  #A file of mutations in patient subdirectory (format: Chr_Ref_Pos_Alt)
  Muts_patient = read.table(paste0(path_prefix,patient,"/All_mutations_filtered.txt"))[,1] 
  #A file of with the sample names belonging to this patient
  case_samples=read.table(paste0(path_prefix,patient,"/samples.txt"))[,1]
  normal_panel = normal_samples
  norm_all_counts = all_counts[normal_panel,,]
  coords_proj = substr(Muts_patient,1,nchar(Muts_patient)-4)
  case_all_counts = all_counts[case_samples,,]
  
  avg_depth_list = c()
  for (i in 1:length(coords_proj)){
    mean_d = sum(case_all_counts[,i,]) / length(case_samples)
    avg_depth_list = c(avg_depth_list, mean_d)
  }
  
  Muts_patient = Muts_patient[avg_depth_list > 10] #Removes remaining globally low coverage variants which were the result of mapping issues
  coords_proj = coords_proj[avg_depth_list > 10]
  
  #Set up pval matrix
  pval_mat = matrix(1,ncol=length(case_samples),nrow=length(Muts_patient))
  rownames(pval_mat)=Muts_patient
  colnames(pval_mat)=case_samples
  
  
  Alt=substr(Muts_patient,nchar(Muts_patient),nchar(Muts_patient))
  Ref=substr(Muts_patient,nchar(Muts_patient)-2,nchar(Muts_patient)-2)
  
  for (s in case_samples){
    rho_est=rep(NA,length(coords_proj))
    test_counts = all_counts[s,coords_proj,]
    for (k in 1:length(coords_proj)) {
      n = sum(test_counts[coords_proj[k],])
      x = test_counts[coords_proj[k],Alt[k]]
      
      N_indiv = rowSums(norm_all_counts[,coords_proj[k],])
      X_indiv = norm_all_counts[,coords_proj[k],c("A","C","G","T")!=Ref[k]]
      pseudo = .Machine$double.eps    
      N=sum(N_indiv)
      X=sum(X_indiv)
      
      mu = max(X,pseudo)/max(N,pseudo)
      counts = cbind(N,X)
      if(is.null(rho)) rho = estimateRho_gridml(counts,mu)
      rdisp = (1 - rho)/rho
      
      prob0 = (X + x)/(N + n); prob0[prob0==0] = pseudo; prob0[prob0==1] = 1 - pseudo
      prob1s = x/(n + pseudo); prob1s[prob1s==0] = pseudo; prob1s[prob1s==1] = 1 - pseudo
      prob1c = X/(N + pseudo); prob1c[prob1c==0] = pseudo; prob1c[prob1c==1] = 1 - pseudo
      
      prob1s = pmax(prob1s,prob1c) # Min error rate is that of the population (one-sided test)
      nu0 = prob0 * rdisp; nu1s = prob1s * rdisp; nu1c = prob1c * rdisp; 
      
      # Likelihood-Ratio Tests
      LL = logbb(x, n, nu0, rdisp) + logbb(X, N, nu0, rdisp) - logbb(x, n, nu1s, rdisp) - logbb(X, N, nu1c, rdisp)
      pvals = pchisq(-2*LL, df=1, lower.tail=F)/2 # We divide by 2 as we are performing a 1-sided test
      # Saving the result
      pval_mat[k,s] = pvals
    } 
  }
  if(is.null(save)){
    return(pval_mat)
  }else{
    write.table(pval_mat,save)
  }
}

#-------------------------------------------------
# Input
#-------------------------------------------------

#setwd('/lustre/scratch117/casm/team294/to3/testes/tumour/allelecount/20200511')
#options(stringsAsFactors = F)

# Vector of normal reference samples, taken from a separate cohort of patients
normal_samples = read.table("/lustre/scratch119/casm/team294rr/to3/testes/tumour/reference_datasets/reference_lcm_normal_panel.txt")[,1] 

# Vector of all samples, including the reference panel and samples of interest
tumour_samples = read.table("/lustre/scratch119/casm/team294rr/to3/testes/tumour/reference_datasets/gct_samples.txt")[,1]
samples <- c(tumour_samples, normal_samples)

# "Bed" file of all mutations to be considered (across all patients)
# Format: Chr Ref Pos Alt
muts = read.table("GCT_SSM.snvs.bed") 
coords = paste(muts$V1, muts$V2, sep="_")

# Read in data from AlleleCounter
all_counts = array(0,dim=c(length(samples), length(coords), 4),
                   dimnames = list(samples, coords, c("A","C","G","T")))

print(length(samples)) #381 samples, 250 reference samples, 131 GCT genomes

for (k in 1:length(samples)){
  #Read in allele counts per sample
  if(file.exists(paste0("output/", samples[k],".txt"))){
    data=read.table(paste0("output/", samples[k],".txt"),comment.char = '',header=T)
    muts_data=paste(data$X.CHR,data$POS,sep="_")
    data=data[!duplicated(muts_data),]
    muts_data=muts_data[!duplicated(muts_data)]
    rownames(data)=muts_data
    all_counts[k,,]=as.matrix(data[coords,3:6])
  }
  print(k)
}

#create list of study patients
#setwd('/lustre/scratch117/casm/team294/to3/testes/tumour/cgpVAF_phase2')
sample.list <- read.table('patients_all.txt', header = F)[,1]

#generate mutations per patient list
patient.list <- unique(substr(sample.list, 0, 7))

for(i in patient.list){
  data <- read.table(paste0(i, '/', list.files(paste0(i, '/'), pattern = '_snp_vaf_nohash.tsv')[1]), comment.char = "", header = T)
  Muts = paste(data$Chrom, data$Pos, data$Ref, data$Alt, sep="_")
  write.table(Muts, paste0(i, '/All_mutations_filtered.txt'), col.names = F, row.names = F, quote = F)
}

#generate a list of samples per tumour
for(i in patient.list){
  write.table(unlist(strsplit(list.files(paste0(i, '/'), pattern = '_snps_geno.final.vcf'), '_snps_geno.final.vcf')), paste0(i, '/samples.txt'), col.names = F, row.names = F, quote = F)
}

# read in the list of study patients
patients = read.table("patients_all.txt")[,1]

#-------------------------------------------------
# Run
#-------------------------------------------------

for (patient in patients){
  shearwater_probability(patient=patient, save=paste0(patient, "/shearwater_snv_pval_mat.txt"))
  print(patient)
}

#get an idea of the mutation burden per sample after this filtering step
post_swater_mb <- c()
for (patient in patients){
  pval_mat = read.table(paste0(patient, "/shearwater_snv_pval_mat.txt"), row.names = 1, sep = ' ', header = T)
  qval_mat = apply(pval_mat,2,function(x) p.adjust(x,method="BH", n = length(as.matrix(pval_mat))))
  post_swater_mb <- c(post_swater_mb, colSums(qval_mat < 0.001, na.rm = T))
}

#BH correction for multiple testing
for (patient in patients){
  pval_mat = read.table(paste0(patient, "/shearwater_snv_pval_mat.txt"), row.names = 1, sep = ' ', header = T)
  qval_mat = apply(pval_mat,2,function(x) p.adjust(x,method="BH", n = length(as.matrix(pval_mat))))
  write.table(qval_mat, paste0(patient, "/shearwater_snv_qval_mat.txt"), row.names = T, sep = '\t', quote = F)
}

#filter vcf to leave only variants that pass the site specific error threshold

library(VariantAnnotation)

for (patient in patients){
  qval_mat = read.table(paste0(patient, "/shearwater_snv_qval_mat.txt"), row.names = 1, sep = '\t', header = T)
  samples = read.table(paste0(patient, "/samples.txt"), sep = '\t', header = F)[,1]
  for (sample in samples){
    data = readVcf(paste0(patient, '/', sample, '_snps_geno.final.vcf'))
    reads = cbind((geno(data)$FAZ)[,2]+(geno(data)$RAZ)[,2],(geno(data)$FCZ)[,2]+(geno(data)$RCZ)[,2],
                  (geno(data)$FGZ)[,2]+(geno(data)$RGZ)[,2],(geno(data)$FTZ)[,2]+(geno(data)$RTZ)[,2])
    Var = data.frame(Chr=as.character(seqnames(rowRanges(data))),
                     Pos=start(ranges(rowRanges(data))),
                     Ref=as.character(ref(data)))
    Alt_tmp = CharacterList(alt(data))
    Var$Alt = as.character(unlist(Alt_tmp))
    Var$NR=rowSums(reads)
    Var$NV=NA
    colnames(reads)=c("A","C","G","T")
    for (k in c("A","C","G","T")){
      Var$NV[Var$Alt==k] = reads[Var$Alt==k,k]
    }
    Var$VAF=Var$NV/Var$NR
    Var$ID = paste(Var$Chr, Var$Pos, Var$Ref, Var$Alt, sep = '_')
    qval_patient_sig_muts = row.names(qval_mat[qval_mat[,sample] < 0.001,])
    swater_filtered_muts = Var[Var$ID %in% qval_patient_sig_muts,]
    write.table(swater_filtered_muts, paste0(patient, '/', sample, '_post_swater_final_snvs.txt'), sep = '\t', quote = F, row.names = F)
  }
}



###Indel filtering

Similar to above but, rather than using allelecounter, we used cgpVAF/exonerate for the pileup. It was also run on a per patient basis rather than a single pileup across the cohort.


#-------------------------------------------------
# Libraries
#-------------------------------------------------

library("GenomicRanges")
library("deepSNV")
library("Rsamtools")

logbb = deepSNV:::logbb
dbetabinom = VGAM::dbetabinom

#-------------------------------------------------
# Functions
#-------------------------------------------------

estimateRho_gridml = function(x, mu) {
  # Estimate rho by MLE grid approach
  rhovec = 10^seq(-6,-0.5,by=0.05) # rho will be bounded within 1e-6 and 0.32
  mm = x[,2]
  cov = c(x[,1])
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=mm, size=cov, rho=rhoj, prob=mu, log=T)))
  rhovec[ll==max(ll)][1]
}

shearwater_indel_probability=function(patient, save=NULL, path_prefix='', rho=10^-3){
  #Function to calculate probability of presence of mutation based on Shearwater
  #'patient' is the name of the patient-specific subfolder
  #'path_prefix' is any prefix to the path necessary to find that subfolder
  #'rho' is the constant for the overdispersion parameter. If rho=NULL, calculate
  # it from the data (much slower)
  #'save' is path for output. If NULL, returns matrix
  # output of this function will be a matrix (muts by samples) of p values
  
  Muts_patient = read.table(paste0(path_prefix, patient, "/All_mutations_filtered.txt"))[,1] #A file of mutations in patient subdirectory (format: Chr_Ref_Pos_Alt)
  Muts_patient = unique(Muts_patient) #to avoid duplicating variants
  case_samples = read.table(paste0(path_prefix, patient,"/samples.txt"))[,1] #A file of with the sample names belonging to this patient
  NR = read.table(paste0(path_prefix, 'NR_all_samples_all_muts.txt'), row.names = 1, header = T, sep = '\t') #total read depth at each variant locus called in this patient across all samples
  NR[NR == 0] = 1 #to prevent NAs
  NV = read.table(paste0(path_prefix, 'NV_all_samples_all_muts.txt'), row.names = 1, header = T, sep = '\t') #variant read depth per variant per sample in patient
  matched_normals = read.table(paste0(path_prefix, 'matched_normals.txt'), header = F, sep = '\t')[,1]
  matched_normal = matched_normals[grepl(patient, matched_normals)] #picks only this patient's matched normal
  
  NR = NR[Muts_patient,] #only mutations relevant to patient
  NV = NV[Muts_patient,]
  
  Depth_filter = rowMeans(NR[, case_samples]) > 10 #Removes remaining globally low coverage variants which were the result of mapping issues
  
  NR = NR[Depth_filter,]
  NV = NV[Depth_filter,]
  
  Min_var_depth = rowSums(NV[, case_samples] >= 5) > 0 #removes low VAF sequencing artefacts. This ensures only variants called in a sample to at least 5 variant reads are kept.
  
  NR = NR[Min_var_depth,]
  NV = NV[Min_var_depth,]
  
  Muts_patient = row.names(NR) #new, slimmer list to run shearwater-like approach on
  
  #Select the normal panel not belonging to this patient with the pileup of sites where variants called in patient's tumour samples
  normal_panel = normal_samples
  norm_NR = NR[, normal_panel]
  norm_NV = NV[, normal_panel]
  
  #Set up pval matrix
  pval_mat = matrix(1,ncol=length(case_samples),nrow=length(Muts_patient))
  rownames(pval_mat)=Muts_patient
  colnames(pval_mat)=case_samples
  
  coords_proj = paste(matrix(unlist(strsplit(Muts_patient, '_')), ncol = 4, byrow = T)[,1], matrix(unlist(strsplit(Muts_patient, '_')), ncol = 4, byrow = T)[,2], sep = '_')
  
  Alt = matrix(unlist(strsplit(Muts_patient, '_')), ncol = 4, byrow = T)[,4]
  Ref = matrix(unlist(strsplit(Muts_patient, '_')), ncol = 4, byrow = T)[,3]
  
  for (s in case_samples){
    rho_est = rep(NA,length(coords_proj))
    test_counts_NR = NR[, s]
    test_counts_NV = NV[, s]
    
    for (k in 1:length(coords_proj)) {
      n = test_counts_NR[k]
      x = test_counts_NV[k]
      N_indiv = norm_NR[k,]
      X_indiv = norm_NV[k,]
      pseudo = .Machine$double.eps    
      N = sum(N_indiv)
      X = sum(X_indiv)
      mu = max(X,pseudo)/max(N,pseudo)
      counts = cbind(N,X)
      
      if(is.null(rho)) rho = estimateRho_gridml(counts,mu)
      rdisp = (1 - rho)/rho
      
      prob0 = (X + x)/(N + n); prob0[prob0==0] = pseudo; prob0[prob0==1] = 1 - pseudo
      prob1s = x/(n + pseudo); prob1s[prob1s==0] = pseudo; prob1s[prob1s==1] = 1 - pseudo
      prob1c = X/(N + pseudo); prob1c[prob1c==0] = pseudo; prob1c[prob1c==1] = 1 - pseudo
      
      prob1s = pmax(prob1s,prob1c) # Min error rate is that of the population (one-sided test)
      nu0 = prob0 * rdisp; nu1s = prob1s * rdisp; nu1c = prob1c * rdisp; 
      
      # Likelihood-Ratio Tests
      LL = logbb(x, n, nu0, rdisp) + logbb(X, N, nu0, rdisp) - logbb(x, n, nu1s, rdisp) - logbb(X, N, nu1c, rdisp)
      pvals = pchisq(-2*LL, df=1, lower.tail=F) / 2 # We divide by 2 as we are performing a 1-sided test
      # Saving the result
      pval_mat[k,s] = pvals
    } 
  }
  pval_mat = pval_mat[((rowSums(norm_NV) / rowSums(norm_NR)) < 0.2) & (NV[, matched_normal]/NR[, matched_normal] < 0.2),] #remove variants seen at a VAF of >0.2 in the normal reference panel or patient's matched normal, these end up being artefactual, even if they are seen at a significantly higher VAF in the tumour samples
  if(is.null(save)){
    return(pval_mat)
  }else{
    write.table(pval_mat, save)
  }
}

#-------------------------------------------------
# Input
#-------------------------------------------------

#setwd('/lustre/scratch117/casm/team294/to3/testes/tumour/shearwater_like_indels/20200523')
#options(stringsAsFactors = F)

#create list of patients
patient.list <- read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/cgpVAF/SNVs/patients_all.txt', header = F)[,1]

# Vector of normal reference samples, no GCT
normal_samples = read.table("/lustre/scratch117/casm/team294/to3/testes/tumour/shearwater_like_indels/20200513/rand_norm_reference.txt")[,1] 

# Vector of all tumour samples, used to create a per patient sample list
tumour_samples = read.table("/lustre/scratch117/casm/team294/to3/testes/tumour/shearwater_like_indels/20200513/gct_samples.txt")[,1]

for (i in unique(substr(tumour_samples, 0, 7))){
  patient.samples = tumour_samples[grepl(paste0(i), tumour_samples)]
  patient.samples = data.frame(patient.samples)
  write.table(patient.samples, paste0(i, '/samples.txt'), sep = '\t', col.names = F, row.names = F, quote = F)
}

#generate a list of mutations per patient in the format Chr_Ref_Pos_Alt
patient_muts_input = unlist(list.files('/lustre/scratch117/casm/team294/to3/testes/tumour/cgpVAF_phase2/', pattern = paste0('_indel_vaf_nohash.tsv')))

for (i in unique(substr(tumour_samples, 0, 7))){
  data <- read.table(paste0(i, '/', i, '.indels.bed'), header = F, sep = '\t')
  Muts = paste(data$V1,data$V2,data$V3,data$V4,sep="_")
  write.table(Muts, paste0(i, '/All_mutations_filtered.txt'), sep = '\t', col.names = F, row.names = F, quote = F)
}

# read in the list of study patients
patient.list <- read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/cgpVAF/SNVs/patients_all.txt', header = F)[,1]

#-------------------------------------------------
# Run
#-------------------------------------------------

for (patient in patient.list){
  shearwater_indel_probability(patient=patient,save=paste0(patient, "/shearwater_pval_indel_mat.txt"), rho = NULL) 
  print(patient)
}

#get an idea of the mutation burden per sample after this filtering step
post_swater_mb <- c()
for (patient in patient.list){
  pval_mat = read.table(paste0(patient, "/shearwater_pval_indel_mat.txt"), sep = ' ', row.names = 1, header = T)
  qval_mat = apply(pval_mat, 2, function(x) p.adjust(x, method="BH", n = length(as.matrix(pval_mat))))
  post_swater_mb <- c(post_swater_mb, colSums(qval_mat <= 0.001, na.rm = T))
}

#BH correction for multiple testing
for (patient in patient.list){
  pval_mat = read.table(paste0(patient, "/shearwater_pval_indel_mat.txt"), row.names = 1, sep = ' ', header = T)
  qval_mat = apply(pval_mat,2,function(x) p.adjust(x,method="BH",n = length(as.matrix(pval_mat))))
  write.table(qval_mat, paste0(patient, "/shearwater_qval_indel_mat.txt"), row.names = T, sep = '\t', quote = F)
}

#filter vcf to leave only variants that pass the site specific error threshold

library(VariantAnnotation)

for (patient in patient.list){
  qval_mat = read.table(paste0(patient, "/shearwater_qval_indel_mat.txt"), row.names = 1, sep = '\t', header = T)
  #qval_mat[is.na(qval_mat)] = 0.5 #replaces NA values should no longer be required
  samples = read.table(paste0(patient, "/samples.txt"), sep = '\t', header = F)[,1]
  for (sample in samples){
    data = readVcf(paste0('/lustre/scratch119/casm/team294rr/to3/testes/tumour/file_repo/SSMs/matched/', sample, '_indels_geno.final.vcf'))
    reads = cbind((geno(data)$MTR)[,2], (geno(data)$WTR)[,2])
    Var = data.frame(Chr = as.character(seqnames(rowRanges(data))), Pos=start(ranges(rowRanges(data))),
                     Ref=as.character(ref(data)))
    Alt_tmp = CharacterList(alt(data))
    Var$Alt = as.character(unlist(Alt_tmp))
    Var$NR = rowSums(reads)
    Var$NV = reads[,1]
    Var$VAF=Var$NV/Var$NR
    Var$ID = paste(Var$Chr, Var$Pos, Var$Ref, Var$Alt, sep = '_')
    qval_patient_sig_muts = row.names(qval_mat[qval_mat[,sample] <= 0.001,])
    swater_filtered_muts = Var[Var$ID %in% qval_patient_sig_muts,]
    write.table(swater_filtered_muts, paste0(patient, '/', sample, '_post_swater_final_indels.txt'), sep = '\t', quote = F, row.names = F)
  }
}
