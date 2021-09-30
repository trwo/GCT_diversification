
#/nfs/users/nfs_t/to3/scripts/hdp_prior_single_chain_updated_20210916.R

options(stringsAsFactors = F)
library(hdp)
lower_threshold=100

n=as.numeric(commandArgs(T)[1])
mutations=read.table("trinuc_mut_mat.txt", header = T)
key_table=read.table("key_table.txt", header = T)

ref = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/signature_extraction/hdp_new_swater/wkdir/final_sigs_2021_09_16.txt', header = T, stringsAsFactors = F, sep = '\t')
ref = ref[, colnames(ref) != 'N1']
prior_sigs = as.matrix(ref)

# number of prior signatures to condition on (8)
nps <- ncol(prior_sigs)

# take a look
#prior_sigs[1:6,1:6]
# columns should sum to one
#colSums(prior_sigs)

# have a look at the first two prior sigs:
#barplot(prior_sigs[,1], las=2, cex.names=0.5)
#barplot(prior_sigs[,'SBS18'], las=2, cex.names=0.5)
#barplot(prior_sigs[,'SBS91'], las=2, cex.names=0.5)

#If requiring a minimum number of mutations:
sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table$Sample%in%sample_remove,]

#Hierarchy is set per patient, can change if wanted
freq = table(key_table$Patient)

#nps <- ncol(prior_sigs)

hdp_prior <- hdp_prior_init(prior_distn = prior_sigs, # matrix of prior sigs
                            prior_pseudoc = rep(1000, nps), # pseudocount weights
                            hh=rep(1, 96), # uniform prior over 96 categories
                            alphaa=c(1, 1), # shape hyperparams for 2 CPs
                            alphab=c(1, 1)) # rate hyperparams for 2 CPs
#hdp_prior
#numdp(hdp_prior)
#pseudoDP(hdp_prior)
#conparam(hdp_prior)
#ppindex(hdp_prior)
#cpindex(hdp_prior)
#dpstate(hdp_prior) # 2 for active node, 1 for frozen, 0 for heldout

# make two more CPs available for the data we will add
hdp_prior <- hdp_addconparam(hdp_prior,
                             alphaa = rep(1,length(freq)+2), # shape hyperparams for x new CPs
                             alphab = rep(1,length(freq)+2)) # rate hyperparams for x new CPs

hdp_prior <- hdp_adddp(hdp_prior,
                       numdp = nrow(mutations) + 1,
                       ppindex = c(1, rep(1+nps+1:(length(freq)), times=freq)),
                       cpindex = c(3, rep(4:(length(freq)+3), times=freq)))

# assign the data to the relevant DP nodes
hdp_prior <- hdp_setdata(hdp_prior,
                         dpindex = (1 + nps + 1) + 1:nrow(mutations), 
                         mutations) # mutation counts in all GCTs


hdp_activated <- dp_activate(hdp_prior, 
                             dpindex = (1+nps+1)+0:nrow(mutations), 
                             initcc = nps+5,
                             seed = n * 1000)

chain=hdp_posterior(hdp_activated,
                    burnin=80000,
                    n=100,
                    seed=n*1000,
                    space=200,
                    cpiter=3)

saveRDS(chain,paste0("hdp_prior_chain_",n,".Rdata"))
