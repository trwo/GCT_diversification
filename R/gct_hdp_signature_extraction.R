#HDP v0.1.5 with priors
##This code is based off of scripts previously written by Nicola Roberts, Tim Coorens and Daniel Leongamornlert.

#1. Arrange prior sigs, set up HDP and run posterior sampling chains

##submit as a job to the cluster
#export PATH=/software/R-3.6.1/bin:${PATH}
#export R_HOME=$(R RHOME)
#export R_LIBS_USER="~/custom_libraries/R/x86_64-pc-linux-gnu-library/3.6_hdp"
#for n in $(seq 1 20);
#do
#bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R'select[mem>10000] rusage[mem=10000]' -M10000 -J $n /software/R-3.6.1/bin/Rscript /nfs/users/nfs_t/to3/scripts/hdp_prior_single_chain.R $n
#done 

#/nfs/users/nfs_t/to3/scripts/hdp_prior_single_chain.R
##########
options(stringsAsFactors = F)
library(hdp)
lower_threshold=100

n=as.numeric(commandArgs(T)[1])
mutations=read.table("trinuc_mut_mat.txt")
key_table=read.table("key_table.txt")

ref = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/MutationalPatterns/reformatted_cosmic_sig_v3.txt', header = T, stringsAsFactors = F, row.names = 1, check.names = F)
ref = t(ref)
cosmic.sigs = read.csv('/lustre/scratch119/casm/team294rr/to3/testes/tumour/MutationalPatterns/PCAWG_sigProfiler_SBS_signatures.csv', header = T, stringsAsFactors = F, row.names = 1, check.names = F)
row.names(cosmic.sigs) = row.names(ref)
prior_sigs = as.matrix(cosmic.sigs)

# number of prior signatures to condition on
nps <- ncol(prior_sigs)

#If requiring a minimum number of mutations:
sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table$Sample%in%sample_remove,]

#Hierarchy is set per patient, can change if wanted
freq = table(key_table$Patient)

nps <- ncol(prior_sigs)

luad_prior <- hdp_prior_init(prior_distn = prior_sigs, # matrix of prior sigs
                             prior_pseudoc = rep(1000, nps), # pseudocount weights
                             hh=rep(1, 96), # uniform prior over 96 categories
                             alphaa=c(1, 1), # shape hyperparams for 2 CPs
                             alphab=c(1, 1)) # rate hyperparams for 2 CPs


# make two more CPs available for the data we will add
luad_prior <- hdp_addconparam(luad_prior,
                              alphaa = rep(1,length(freq)+2), # shape hyperparams for x new CPs
                              alphab = rep(1,length(freq)+2)) # rate hyperparams for x new CPs

luad_prior <- hdp_adddp(luad_prior,
                        numdp = nrow(mutations) + 1,
                        ppindex = c(1, rep(1+nps+1:(length(freq)), times=freq)),
                        cpindex = c(3, rep(4:(length(freq)+3), times=freq)))

# assign the data to the relevant DP nodes
luad_prior <- hdp_setdata(luad_prior,
                          dpindex = (1+nps+1)+1:nrow(mutations), 
                          mutations) # mutation counts in all GCTs


hdp_activated <- dp_activate(luad_prior, 
                             dpindex = (1+nps+1)+0:nrow(mutations), 
                             initcc = nps+5,
                             seed = n * 1000)

chain=hdp_posterior(hdp_activated,
                    burnin=40000,
                    n = 200,
                    seed = n * 1000,
                    space = 300,
                    cpiter = 3)

saveRDS(chain,paste0("hdp_prior_chain_",n,".Rdata"))
##########

#2. Review diagnostic plots, extract consensus components and signatures

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/signature_extraction/hdp_0.1.5')

options(stringsAsFactors = F)
library(hdp)

chlist <- vector("list", 20)
for (i in 1:20){
  if(file.exists(paste0("hdp_prior_chain_",i,".Rdata"))){
    chlist[[i]] <- readRDS(paste0("hdp_prior_chain_",i,".Rdata"))
  }
}
if(any(unlist(lapply(chlist,is.null)))) chlist=chlist[-which(unlist(lapply(chlist,is.null)))]

luad_multi <- hdp_multi_chain(chlist)
pdf("QC_plots_chain.pdf") 
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(luad_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(luad_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(luad_multi), plot_data_assigned, bty="L")
dev.off()

#need to remove chains 4, 5, 6, 9, 12, 16, 17

chlist <- vector("list", 20)
for (i in c(1:3, 7, 8, 10, 11, 13:15, 18:20)){
  if(file.exists(paste0("hdp_prior_chain_",i,".Rdata"))){
    chlist[[i]] <- readRDS(paste0("hdp_prior_chain_",i,".Rdata"))
  }
}
if(any(unlist(lapply(chlist,is.null)))) chlist=chlist[-which(unlist(lapply(chlist,is.null)))]

luad_multi <- hdp_multi_chain(chlist)

pdf("QC_plots_chain_after_review.pdf") 
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(luad_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(luad_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(luad_multi), plot_data_assigned, bty="L")
dev.off()

luad_multi <- hdp_extract_components(luad_multi)
saveRDS(luad_multi,"HDP_multi_chain.Rdata")

numcomp(luad_multi) #6
prop.ex(luad_multi) #0.91

pdf("muts_attributed.pdf")
plot_comp_size(luad_multi, bty="L")
dev.off()

# plot components / signatures
# pick your colours
mut_colours=c("dodgerblue","black","red","grey70","olivedrab3","plum2")
# labels along bottom x-axis
trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
# group labels along the top (and controls colour grouping)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))

ref = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/MutationalPatterns/reformatted_cosmic_sig_v3.txt', header = T, stringsAsFactors = F, row.names = 1, check.names = F)
ref = t(ref)
cosmic.sigs = read.csv('/lustre/scratch119/casm/team294rr/to3/testes/tumour/MutationalPatterns/PCAWG_sigProfiler_SBS_signatures.csv', header = T, stringsAsFactors = F, row.names = 1, check.names = F)
row.names(cosmic.sigs) = row.names(ref)
prior_sigs = as.matrix(cosmic.sigs)

# number of prior signatures to condition on
nps <- ncol(prior_sigs)

pdf(paste0("hdp_components_with_priors.pdf"),width=12,height=4)
plot_comp_distn(luad_multi, cat_names=trinuc_context,
                grouping=group_factor, col=mut_colours,
                col_nonsig="grey97", show_group_labels=TRUE)
dev.off()

lower_threshold = 100
mutations=read.table("/lustre/scratch119/casm/team294rr/to3/testes/tumour/HDP/input/trinuc_mut_mat.txt")
key_table=read.table("/lustre/scratch119/casm/team294rr/to3/testes/tumour/HDP/input/key_table.txt")
sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table$Sample%in%sample_remove,]


pdf("signature_attribution.pdf",width=10,height=8)
plot_dp_comp_exposure(luad_multi, dpindices = (1+nps+1)+1:nrow(mutations), incl_nonsig = T,
                      col=c('black', RColorBrewer::brewer.pal(numcomp(luad_multi), "Set1")), cex.names=0.8,
                      ylab_exp = 'Signature exposure', leg.title = 'Signature')
dev.off()

mean_assignment=as.data.frame(comp_dp_distn(luad_multi)$mean)
write.table(mean_assignment,"mean_assignment_prior_hdp.txt")
mean_sigs=as.data.frame(t(comp_categ_distn(luad_multi)$mean))
write.table(mean_sigs,"hdp_prior_sigs.txt")

#3.Deconvolute 'novel' signatures

#Deconvolute HDP signatures into reference signatures and arrive at final set

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/signature_extraction/hdp_0.1.5')

options(stringsAsFactors = F)
library(hdp)
library(RColorBrewer)
library(lsa)
library(lattice)

mut.cols = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)

#Load HDP signatures
mean_sigs = read.table("hdp_prior_sigs.txt", row.names = 1, header = T)
hdp_sigs = mean_sigs
#Load reference signatures, need reformatting
ref = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/MutationalPatterns/reformatted_cosmic_sig_v3.txt', header = T, stringsAsFactors = F, row.names = 1, check.names = F)
ref = t(ref)
cosmic.sigs = read.csv('/lustre/scratch119/casm/team294rr/to3/testes/tumour/MutationalPatterns/PCAWG_sigProfiler_SBS_signatures.csv', header = T, stringsAsFactors = F, row.names = 1, check.names = F)
row.names(cosmic.sigs) = row.names(ref)
prior_sigs = as.matrix(cosmic.sigs)
ref = prior_sigs

#Assess cosine similarities for all unassigned signatures (X0, N1, N2)
cosine_matrix=data.frame(matrix(nrow=ncol(hdp_sigs[, c('X0', 'N1', 'N2')]), ncol=ncol(ref)))
rownames(cosine_matrix)=colnames(hdp_sigs[, c('X0', 'N1', 'N2')])
colnames(cosine_matrix)=colnames(ref)

for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n,m] <- cosine(x=hdp_sigs[,rownames(cosine_matrix)[n]],
                                 y=ref[,colnames(cosine_matrix)[m]])
  }
}

write.table(cosine_matrix, "Cosine_similarities_unassigned_to_priors.txt",sep="\t",quote=F)

#cosine_matrix = read.table('Cosine_similarities_unassigned_to_priors.txt', header = T, sep = '\t')
#rownames(cosine_matrix) = c('X0', 'N1', 'N2')
hdp_sigs = hdp_sigs[, c('X0', 'N1', 'N2')]

#plot output
pdf("cosine_similarities_unassigned_to_priors.pdf", height=5, width=15)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

#select which signatures to decompose into reference signatures
sigs_to_decompose=rowSums(cosine_matrix>0.85)
for(n in 1:nrow(cosine_matrix)){
  print(paste0(rownames(cosine_matrix)[n],": ",paste(colnames(cosine_matrix)[cosine_matrix[n,]>0.7],collapse=",")))
}

#N1 has large error bars and does not look stable across the chains

#"0: SBS3,SBS5,SBS40" stable, flat signature, closest to SBS5, followed by SBS40, then SBS3
#"N1: SBS3,SBS5,SBS39,SBS40" assortment of flat signature, unstable
#"N2: SBS31,SBS35" platinum signature

#First iteration; decomposed hdp sigs into all suspected sigs 
gdsigs=c("SBS3", "SBS5", "SBS31", "SBS35", "SBS39", "SBS40") #0.7 threshold

signatures=t(ref[,gdsigs])
sample_list= c('X0', 'N1', 'N2')
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

#4. Plot initial deconvolution and save results
pdf("Deconvolution_hdp_sigs_after_prior_R1.pdf", height=5, width=10)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot((signature_fraction[nrow(signature_fraction):1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

write.table(signature_fraction, "hdp_unassigned_sigs_broken_down_into_cosmic_sigs.txt", sep="\t", col.names=T, row.names = T, quote=F)
