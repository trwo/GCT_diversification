
#Adapted from code by Tim Coorens (https://github.com/TimCoorens/ClonalNephrogenesis).

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
  
  for (i in coords_proj){
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
patients = read.table("/lustre/scratch117/casm/team294/to3/testes/tumour/cgpVAF_phase2/patients_all.txt")[,1]

#-------------------------------------------------
# Run
#-------------------------------------------------

for (patient in patients){
  shearwater_probability(patient=patient, save=paste0("/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/variant_calls/sub_shearwater/", patient, "/shearwater_snv_pval_mat.txt"))
  print(patient)
}

#get an idea of the mutation burden per sample after this filtering step
post_swater_mb <- c()
for (patient in patients){
  pval_mat = read.table(paste0("/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/variant_calls/sub_shearwater/", patient, "/shearwater_snv_pval_mat.txt"), row.names = 1, sep = ' ', header = T)
  qval_mat = apply(pval_mat,2,function(x) p.adjust(x,method="BH", n = length(as.matrix(pval_mat))))
  post_swater_mb <- c(post_swater_mb, colSums(qval_mat < 0.001, na.rm = T))
}

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/variant_calls/sub_shearwater/')

#BH correction for multiple testing
for (patient in patients){
  pval_mat = read.table(paste0(patient, "/shearwater_snv_pval_mat.txt"), row.names = 1, sep = ' ', header = T)
  qval_mat = apply(pval_mat,2,function(x) p.adjust(x,method="BH", n = length(as.matrix(pval_mat))))
  write.table(qval_mat, paste0(patient, "/shearwater_snv_qval_mat.txt"), row.names = T, sep = '\t', quote = F)
}


#setwd('/lustre/scratch117/casm/team294/to3/testes/tumour/cgpVAF_phase2/') taken files from here

#for(patient in patients){
#  system(paste0('cp ', patient, '/*_snps_geno.final.vcf /lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/variant_calls/sub_shearwater/', patient))
#}

#for(patient in patients){
#  system(paste0('cp ', patient, '/samples.txt /lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/variant_calls/sub_shearwater/', patient))
#}

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
