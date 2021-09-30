
#Adapted from code by Tim Coorens  (https://github.com/TimCoorens/Polymerase/tree/master/Signatures)

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/signature_extraction/hdp_new_swater/novel_sig_cleanup')

options(stringsAsFactors = F)
library(hdp)
library(RColorBrewer)
library(lsa)
library(lattice)

mut.cols = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)

#Load HDP signatures
mean_sigs = read.table("hdp_prior_sigs.txt", row.names = 1, header = T)
hdp_sigs = mean_sigs
#Load reference signatures, need formatting
ref = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/reference_datasets/COSMIC_v3.2_SBS_GRCh37.txt', header = T, sep = '\t', stringsAsFactors = F, row.names = 1)

#Assess cosine similarities for all unassigned signatures (X0, N1, N2)
cosine_matrix=data.frame(matrix(nrow=ncol(hdp_sigs[, c('N1', 'N2', 'N3', 'N4')]), ncol=ncol(ref)))
rownames(cosine_matrix)=colnames(hdp_sigs[, c('N1', 'N2', 'N3', 'N4')])
colnames(cosine_matrix)=colnames(ref)

for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n,m] <- cosine(x=hdp_sigs[,rownames(cosine_matrix)[n]],
                                 y=ref[,colnames(cosine_matrix)[m]])
  }
}

write.table(cosine_matrix, "Cosine_similarities_unassigned_to_priors.txt",sep="\t",quote=F)

hdp_sigs = hdp_sigs[, c('N1', 'N2', 'N3', 'N4')]

#plot output
pdf("cosine_similarities_unassigned_to_priors.pdf", height=5, width=15)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

#select which signatures to decompose into reference signatures
sigs_to_decompose=rowSums(cosine_matrix>0.85)
#N3 is SBS35
#N4 is SBS17a

excluded_sigs = c('N3', 'N4') #already a close match
sigs_close_cosine = c()
for(n in 1:nrow(cosine_matrix)){
  if(!row.names(cosine_matrix)[n] %in% excluded_sigs){
    sigs_close_cosine = c(sigs_close_cosine, colnames(cosine_matrix)[cosine_matrix[n,]>0.7])
    print(paste0(rownames(cosine_matrix)[n],": ",paste(colnames(cosine_matrix)[cosine_matrix[n,]>0.7],collapse=",")))
  }
}

#[1] "N1: SBS3,SBS5,SBS40"
#"N2: "

#First iteration; decomposed hdp sigs into all suspected sigs 
gdsigs = sigs_close_cosine

signatures=t(ref[,gdsigs])
sample_list= c('N1', 'N2') #only keep N2 in to prevent errors from vector <-> dataframe conversions later
profiles=hdp_sigs[,sample_list]

signature_fraction = matrix(NA,nrow=nrow(signatures),ncol=length(sample_list))
rownames(signature_fraction) = rownames(signatures)
colnames(signature_fraction) = sample_list
maxiter <- 1000

for (j in 1:length(sample_list)) {
  freqs = profiles[,j]
  freqs[is.na(freqs)] = 0
  # EM algorithm with to estimate the signature contribution
  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  # Saving the signature contributions for the sample
  print(j/length(sample_list))
  signature_fraction[,j] = alpha
}

#Plot initial deconvolution and save results
pdf("Deconvolution_hdp_sigs_after_prior_R1.pdf", height=5, width=10)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot((signature_fraction[nrow(signature_fraction):1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

write.table(signature_fraction[,'N1'], "hdp_unassigned_sig_broken_down_into_cosmic_sigs.txt", sep="\t", col.names=T, row.names = T, quote=F)

sigs_deconv_R2=list()

for(n in 1:length(sample_list)){
  sigs_deconv_R2[[n]]=rownames(signature_fraction)[signature_fraction[,n]>0.2]
}

names(sigs_deconv_R2)=colnames(signature_fraction)

#No signature has one major contributor only

sigs_to_deconv=names(sigs_deconv_R2)[unlist(lapply(sigs_deconv_R2,length))>1]

ref_sigs_R2=sort(unique(unlist(sigs_deconv_R2)))
signature_fractionR2=matrix(NA,ncol=length(sigs_to_deconv),nrow=length(ref_sigs_R2))
rownames(signature_fractionR2)=ref_sigs_R2
colnames(signature_fractionR2)=sigs_to_deconv

#repeat the deconvolution with the identified constitutive signatures
n = 1
for(s in sigs_to_deconv){
  gdsigs <- sigs_deconv_R2[[s]]
  signatures <- t(ref[,gdsigs])
  
  signature_fraction = matrix(NA,nrow=nrow(signatures),ncol=length(sample_list))
  rownames(signature_fraction) = rownames(signatures)
  colnames(signature_fraction) = sample_list
  maxiter <- 1000
  
  freqs = profiles[,s]
  freqs[is.na(freqs)] = 0
  
  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  # Saving the signature contributions for the sample
  signature_fractionR2[gdsigs,n] = alpha
  n=n+1
  reconsbs <- rep(0,96)
  for (g in gdsigs) {
    reconsbs=reconsbs+(ref[,g]*alpha[g])
  }
  cosine_reconst=cosine(x=reconsbs, y=hdp_sigs[,s])
  print(paste0(s,": ",cosine_reconst))
  pdf(paste0("HDP_",s,"_reconstitution2.pdf"))
  par(mfrow=c(length(alpha)+2,1))
  par(mar=c(1,2,4,1))
  barplot(hdp_sigs[,s], col=mut.cols, main=paste0("HDP ",s),names.arg="")
  barplot(reconsbs, col=mut.cols, main=paste0("Reconstituted ",s," cosine similarity to original: ", round(cosine_reconst, digits=2)))
  for (g in gdsigs) {
    barplot(ref[,g], col=mut.cols, main=paste0("PCAWG ", g, " accounts for ", round(alpha[g], digits=2)))
  }
  dev.off()
}

#complete sig list
sigs_deconv_R2$N2 = 'N2'
sigs_deconv_R2$N3 = 'SBS35'
sigs_deconv_R2$N4 = 'SBS17a'

saveRDS(sigs_deconv_R2,"novel_hdp2refsigs.Rdata")

ref_sigs_R2 = c(ref_sigs_R2, "SBS18", "SBS1", "SBS17a", "SBS17b", 'SBS35', "SBS91") #lost SBS31, gained SBS17a

#Combine the hdp signature that did not get deconvolved and reference signatures into final table
final_sigs = cbind(hdp_sigs[,c("N2")],ref[,ref_sigs_R2])
#Rename the HDP components that didn't get deconvoluted
colnames(final_sigs)[1] = "N2"

#Order them - optional
final_sigs = final_sigs[, c('N2', 'SBS1', 'SBS5', 'SBS17a', 'SBS17b', 'SBS18', 'SBS35', 'SBS40', 'SBS91')]

write.table(final_sigs, "final_sigs_CLEANED_UP_2021_09_16.txt", col.names = T, row.names = T, sep = '\t', quote = F)
