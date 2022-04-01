#WGD estimation analysis
#Code for analysis adapted from work by Moritz Gerstung and Santiago Gonzalez

#1) FUNCTIONS

#Import some functions from PCAWG paper
#source('/nfs/users/nfs_t/to3/scripts/additional_mutationtimeR_functions.R') #source extra, unchanged functions

#edited functions from mutationtimeR
findMainCluster <- function(bb, min.dist = 0.05){ #at least 20 SNVs
  w <- which(bb$n.snv_mnv >= 20 & !is.na(bb$time) & ((bb$time.up - bb$time.lo) < 0.5)) #edited to exclude regions with very wide confidence intervals when estimating WGD time, to3 03/03/21
  s <- seq(0,1,0.001) #increased granularity in view of how early WGD is, to3 03/03/21
  l2 <- pmin(bb$time.lo, bb$time - min.dist)[w]
  u2 <- pmax(bb$time.up, bb$time + min.dist)[w]
  l1 <- (l2 +  bb$time[w])/2
  u1 <- (u2 +  bb$time[w])/2
  wd <- as.numeric(width(bb)[w])
  o <- sapply(s, function(i) sum(wd * ( (l2 <= i & u2 >=i) + (l1 <= i & u1 >= i))))
  s[which.max(o)]
}

#remove the deamination component as we have too few mutations, especially associated with SBS1 to do this and filter low confidence timing segments from analysis
computeWgdParam <- function(vcf, bb, clusters, purity, sex){ #added sex in
  # 1. Find segments compatible with WGD
  min.dist <- 0.05
  m <- findMainCluster(bb) #modifed as above to higher confidence regions on which to time WGD, to3 03/03/21
  l <- pmin(bb$time.lo, bb$time - min.dist)
  u <- pmax(bb$time.up, bb$time + min.dist)
  o <- which(l <= m & u >= m & ((bb$time.up - bb$time.lo) < 0.5)) #only going to count and time over segments with higher confidence, to3 03/03/21
  #note that this line has been left as is in fractionGenomeWgdCompatible as there we only want to know the regions compatable with WGD, not necessarily depend on them all for downstream analyses, to3 03/03/21
  # 2. Find substitutions in compatible segments
  intersect_seg = findOverlaps(vcf, bb[o], ignore.strand = T) #%over% method not working so this is a workaround, to3 03/03/21
  vcf <- vcf[queryHits(intersect_seg),] 
  w <- which(info(vcf)$MajCN==2 & sapply(info(vcf)$CNID, length)==1)
  v <- vcf[w]
  if(nrow(v)< 100) return(NULL) # At least 100 SNVs total over which to estimate the time of acquisition
  seqnames(rowRanges(v)) <- factor(3-info(v)$MinCN, levels=seqlevels(v))
  v <- sort(v)
  
  # 3. Merged CN segments
  b <- GRanges(1:3, IRanges(rep(1,3),rep(max(end(v)),3)), copy_number=4:2, major_cn=2, minor_cn=2:0, clonal_frequency=as.numeric(purity))
  
  # 4. Calculate times
  l <- computeMutCn(v, b, clusters, purity, isWgd=TRUE, n.boot=200, rho=0.01, xmin=3, gender = sex) #recalculates times calculated before on single segments
  b$n.snv_mnv <- l$n <- table(factor(info(v)$MinCN, levels=2:0))
  l$time <- bbToTime(b, l$P)
  
  row.names(l$D) = row.names(info(v)) #get sub annotations back, to3 03/03/21
  l$cn_state_adj_burden = sum(l$D$MutCN == 2, na.rm = T) + (sum(l$D$MutCN == 2, na.rm = T) * ((unname(l$n['0']) +  unname(l$n['1'])) /  (unname(l$n['0']) +  unname(l$n['1']) + unname(l$n['2'])))) #correct for the 2+1 and 2+0 regions where pre-duplication mutations on the minor allele cannot be counted, to3 04/03/21
  l$hq_seg_width = sum(width(bb[o][bb[o]$major_cn == 2])) #over what length of the genome were these substitutions identified? to3 03/03/21
  l$adj_genome_sub_burden = l$cn_state_adj_burden * (genome_length / l$hq_seg_width) #number of mutations estimated to be on both copies, corrected for genome coverage by good confidence WGD segments, to3 04/03/21
  
  return(l)
}

#2) INPUT

#Create lists of pertinent input data, mostly dervied from mutationtimeR output

library(VariantAnnotation)

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/mutationtimeR_new_swater/wkdir')

chrom_sizes <- read.table('/lustre/scratch119/casm/team294rr/rs30/References/GRCh37d5/chrom.sizes.txt', header = F, sep = '\t', stringsAsFactors = F)
genome_length = sum(chrom_sizes$V2) #3095677412

bb_files = list.files(pattern = '_post_mutationtime_segments.rds')

finalBB <- list()
for( f in bb_files){
  g = readRDS(f)
  finalBB[[unlist(strsplit(f, '_post_mutationtime_segments.rds'))]] <- g
}

snv_files = list.files(pattern = '_ssm_post_mutationtime.vcf')

finalSnv = list()
for(f in snv_files){
  vcf = readVcf(f)
  finalSnv[[unlist(strsplit(f, '_ssm_post_mutationtime.vcf'))]] <- vcf
}

cluster_files = list.files(pattern = '_10000iters_3000burnin_bestClusterInfo.txt')

finalClusters <- list()
for(f in cluster_files){
  clusters <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)
  sample = unlist(strsplit(f, '_10000iters_3000burnin_bestClusterInfo.txt'))
  bb = read.table(paste0(sample, '_post_mutationtime_segments.txt'), sep = '\t', header = T, stringsAsFactors = F)
  clusters$proportion = clusters$location * max(bb$clonal_frequency)
  clusters <- clusters[, c('cluster.no', 'proportion', 'no.of.mutations')]
  names(clusters) <- c('cluster', 'proportion', 'n_ssms')
  finalClusters[[sample]] <- clusters
}

purity_files = list.files(pattern = '_rho_and_psi.txt')

finalPurity <- numeric()
for(f in purity_files){
  pur_df = read.table(f, sep = '\t', header = T, stringsAsFactors = F)
  purity = pur_df['FRAC_GENOME', 'rho']
  finalPurity[unlist(strsplit(f, '_rho_and_psi.txt'))] = purity
}

manifest <- read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/supplementary_data/WGS_summary_metadata.txt', sep = '\t', header = T, stringsAsFactors = F)
filepath = '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/mutationtimeR_new_swater/wkdir'
sample.list <- unlist(strsplit(list.files(paste0(filepath), pattern = '_10000iters_3000burnin_bestClusterInfo.txt'), '_10000iters_3000burnin_bestClusterInfo.txt'))
reference_metadata <- manifest[manifest$Sample %in% c(sample.list), c('Sample', 'Sex')]
row.names(reference_metadata) = c(1:length(reference_metadata$Sample))


sex_vec = reference_metadata$Sex
names(sex_vec) = reference_metadata$Sample

#3) WGD TIMING IN MUTATION TIME

#Sanity check for final, eligible samples - are they consistent with WGD?

finalPloidy <- sapply(finalBB, averagePloidy)
names(finalPloidy) <- names(finalBB)

finalHom <- sapply(finalBB, averageHom)
names(finalHom) <- names(finalBB)

isWgd <- .classWgd(finalPloidy, finalHom) #2.9 -2*hom <= ploidy
table(isWgd) #all samples classed as WGD based on ploidy and homozygosity

fracGenomeWgdComp <- t(sapply(finalBB, function(bb) {
  fgw <- try(fractionGenomeWgdCompatible(bb)); 
  if(class(fgw)!='try-error') fgw
  else rep(NA,10)}))
rownames(fracGenomeWgdComp) <- names(finalBB)
fracGenomeWgdComp = data.frame(fracGenomeWgdComp)
fracGenomeWgdComp$patient = substr(row.names(fracGenomeWgdComp), 0, 7)
aggregate(time.wgd ~ patient, fracGenomeWgdComp, median)

#   patient time.wgd
#  PD42036   0.0200
#  PD43296   0.0170
#  PD43298   0.0485
#  PD43299   0.3950
#  PD45543   0.0430
#  PD45544   0.0010
#  PD45545   0.0150
#  PD46269   0.0010
#  PD46271   0.0175
# PD46966   0.0380
# PD46969   0.0010
# PD50654   0.0160
# PD50655   0.8140
# PD50656   0.5680
# PD50657   0.0520
# PD50660   0.4890
# PD50661   0.2290

#4) ESTIMATE PRE-DUPLICATION SUB BURDEN
wgdParam <- mclapply(names(finalSnv)[isWgd], function(ID){
  try(computeWgdParam(finalSnv[[ID]], finalBB[[ID]], clusters=finalClusters[[ID]], purity=finalPurity[ID], sex=sex_vec[ID]))
})
names(wgdParam) <- names(finalSnv)[isWgd]
void <- sapply(wgdParam, function(x) is.null(x) | class(x)=="try-error") #none

saveRDS(wgdParam, 'mutationtimeR_WGD_segments_20210918.rds')

predup_burden_df = data.frame(sample = names(wgdParam))
row.names(predup_burden_df) = predup_burden_df$sample

predup_burden_df$adj_genome_sub_burden = NA
predup_burden_df$hq_seg_width = NA

for(i in 1:nrow(predup_burden_df)){
  predup_burden_df$adj_genome_sub_burden[i] = wgdParam[[predup_burden_df$sample[i]]]$adj_genome_sub_burden
  predup_burden_df$hq_seg_width[i] = wgdParam[[predup_burden_df$sample[i]]]$hq_seg_width
}

predup_burden_df$patient = substr(predup_burden_df$sample, 0, 7)

write.table(fracGenomeWgdComp, 'gct_wgd_estimate_mut_time_20210918.txt', col.names = T, row.names = T, sep = '\t', quote = F)
write.table(predup_burden_df, 'gct_predup_sub_burden_estimate_20210918.txt', col.names = T, row.names = F, sep = '\t', quote = F)


#5) PLOT AGE AGAINST PRE-DUPLICATION SUB BURDEN (FIG. 2A)

library(reshape2)
library(ggplot2)
library(ggpubr)
library(VariantAnnotation)
library(ggbeeswarm)
library(investr)

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/mutationtimeR_new_swater/wkdir')

fracGenomeWgdComp = read.table('gct_wgd_estimate_mut_time_20210918.txt', header = T, sep = '\t', stringsAsFactors = F)
predup_burden_df = read.table('gct_predup_sub_burden_estimate_20210918.txt', header = T, sep = '\t', stringsAsFactors = F)

predup_burden_df_per_patient = aggregate(adj_genome_sub_burden ~ patient, predup_burden_df, mean)

#  patient adj_genome_sub_burden
#  PD42036              2.880271
#  PD43296              8.520985
#  PD43298             17.974153
#  PD43299            399.311492
#  PD45543              0.000000
#  PD45544              0.000000
#  PD45545              8.143693
#  PD46269              2.056881
#  PD46271              5.366702
#  PD46966              6.169503
#  PD46969             11.098779
#  PD50654             29.871324
#  PD50655            788.782149
#  PD50656            312.862374
#  PD50657             10.601263
#  PD50660            399.891822
#  PD50661            359.694023

pt_ages = manifest[!duplicated(manifest$Case.ID), c('Case.ID', 'Age', 'Diagnosis', 'Organ')]

predup_burden_df_per_patient$age = NA
predup_burden_df_per_patient$dx = NA
predup_burden_df_per_patient$site = NA
for(i in 1:nrow(predup_burden_df_per_patient)){
  predup_burden_df_per_patient$age[i] = pt_ages[pt_ages$Case.ID == predup_burden_df_per_patient$patient[i],]$Age
  predup_burden_df_per_patient$dx[i] = pt_ages[pt_ages$Case.ID == predup_burden_df_per_patient$patient[i],]$Diagnosis
  predup_burden_df_per_patient$site[i] = pt_ages[pt_ages$Case.ID == predup_burden_df_per_patient$patient[i],]$Organ
}

median(predup_burden_df_per_patient[predup_burden_df_per_patient$age > 12,]$adj_genome_sub_burden) #5.768102
range(predup_burden_df_per_patient[predup_burden_df_per_patient$age > 12,]$adj_genome_sub_burden) #0.00000 17.97415
median(predup_burden_df_per_patient[predup_burden_df_per_patient$age <= 12,]$adj_genome_sub_burden) #359.694
range(predup_burden_df_per_patient[predup_burden_df_per_patient$age <= 12,]$adj_genome_sub_burden) #10.60126 788.78215

ggplot(predup_burden_df_per_patient) + geom_point(mapping = aes(x = age, y = adj_genome_sub_burden, col = dx, shape = site)) +
  theme_pubr() + ylab('Estimated pre-duplication substitution burden') + xlab('Age (years)') +
  theme(legend.position = 'right')

#all models
asympt.model = nls(adj_genome_sub_burden ~ SSasymp(age, Asym, R0, lrc), predup_burden_df_per_patient) #asymptotic regression
null.model = lm(formula = adj_genome_sub_burden ~ 1, data = predup_burden_df_per_patient) # null model
#predup_burden_df_per_patient$log.adj_genome_sub_burden = log(predup_burden_df_per_patient$adj_genome_sub_burden + 1)
#log.model = lm(formula = log.adj_genome_sub_burden ~ age, data = predup_burden_df_per_patient) # log model
simple.linear.model = lm(formula = adj_genome_sub_burden ~ age, data = predup_burden_df_per_patient) #simple linear

#identify best fit
AIC(asympt.model, null.model, simple.linear.model) #asymptotic regression fits best

summary(asympt.model)

#generate confidence intervals for asymptotic regression model
age_range = seq(0.75, 58, 0.25) #new age range to generate interval over
cimat <- as.data.frame(investr::predFit(asympt.model, data.frame(age = age_range), interval = "confidence", level = 0.95))
cimat$age = age_range

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/figures/age_vs_pre_dup_burden_plot_20210919.pdf', height = 3,  width = 4, useDingbats = F)
ggplot() +
  theme_pubr() + ylab('Estimated pre-duplication substitution burden') + xlab('Age (years)') +
  theme(legend.position = 'none') +
  geom_ribbon(data = cimat, mapping = aes(x = age, ymin = lwr, ymax = upr), color = NA, alpha = 0.1, fill = "grey50") +
  #geom_point(predup_burden_df_per_patient, mapping = aes(x = age, y = adj_genome_sub_burden, fill = dx, shape = site), col = 'black', size = 5) + 
  geom_point(predup_burden_df_per_patient, mapping = aes(x = age, y = adj_genome_sub_burden, fill = dx), shape = 21, col = 'black', size = 5) + 
  stat_smooth(data = predup_burden_df_per_patient, mapping = aes(x = age, y = adj_genome_sub_burden), method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = F, lty = 2, col = 'red') +
  scale_fill_manual(values = c('#9d9d9c', '#ffffff', '#1d1d1b'))
#+ scale_shape_manual(values = c(21, 24, 22))
dev.off()

#6) GCT WGD TIMING VS PCAWG (FIG. 2B)

library(RColorBrewer)
library(ggplot2)

summary_data = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/supplementary_data/WGS_summary_metadata.txt', header = T, sep = '\t', stringsAsFactors = F)
row.names(summary_data) = summary_data$Sample

#read in PCAWG data
pcawg_wgd = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/original_nature_submission/reference_datasets/20190824_timing_gains_adj.txt', header = T, sep = '\t', stringsAsFactors = F)

pcawg_wgd[pcawg_wgd$samplename %in% pcawg_wgd[pcawg_wgd$type == 'WGD',]$samplename,]


sample_hist_key = pcawg_wgd[!duplicated(pcawg_wgd[, c('samplename', 'histology_abbreviation')]), c('samplename', 'histology_abbreviation')]
pcawg_wgd = pcawg_wgd[pcawg_wgd$histology_abbreviation %in% sort(unique(sample_hist_key$histology_abbreviation))[table(sample_hist_key$histology_abbreviation) >= 10], ] #only keep tumour types with at least 10 samples to make proportion values sensible

sample_hist_key = sample_hist_key[sample_hist_key$histology_abbreviation %in% sort(unique(sample_hist_key$histology_abbreviation))[table(sample_hist_key$histology_abbreviation) >= 10], ] #only do it for tumour types with at least 10 tumours represented

sample_hist_key$WGD = NA
sample_hist_key$WGD_timing = NA

for(i in 1:nrow(sample_hist_key)){
  sample_hist_key$WGD[i] = 'No'
  if(sample_hist_key$samplename[i] %in% pcawg_wgd[pcawg_wgd$chr == 'WGD',]$samplename){
    sample_hist_key$WGD[i] = 'Yes'
    sample_hist_key$WGD_timing[i] = pcawg_wgd[pcawg_wgd$samplename == sample_hist_key$samplename[i] & pcawg_wgd$chr == 'WGD',]$time
  } 
}

#read in GCT data
timing_df = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/mutationtimeR_new_swater/wkdir/gct_wgd_estimate_mut_time_20210918.txt', sep = '\t', header = T, stringsAsFactors = F)

median_wgd = aggregate(time.wgd ~ patient, timing_df, median)
median_wgd$age = NA
for(i in 1:nrow(median_wgd)){
  median_wgd$age[i] = unique(summary_data[summary_data$Case.ID == median_wgd$patient[i],]$Age)
}

median_wgd_postpubertal = median_wgd[median_wgd$age > 12,]
median_wgd_postpubertal$histology_abbreviation = 'Postpubertal_GCT'
median_wgd_postpubertal = median_wgd_postpubertal[, c(1,2,4)]

gct_samples = data.frame(samplename = median_wgd_postpubertal$patient, histology_abbreviation = rep('Postpubertal_GCT', nrow(median_wgd_postpubertal)), WGD = rep('Yes', nrow(median_wgd_postpubertal)), WGD_timing = median_wgd_postpubertal$time.wgd)

median(median_wgd_postpubertal[median_wgd_postpubertal$patient %in% c('PD45543', 'PD45544', 'PD45545'),]$time.wgd) #post-pubertal ovary 0.015
median(median_wgd_postpubertal[!median_wgd_postpubertal$patient %in% c('PD45543', 'PD45544', 'PD45545'),]$time.wgd) #post-pubertal testis 0.0175

#combined data gct and pcawg data
sample_hist_key = rbind(sample_hist_key, gct_samples)

sample_hist_key$WGD_timing_plot = round(sample_hist_key$WGD_timing, digits = 2)
sample_hist_key[is.na(sample_hist_key$WGD_timing_plot),]$WGD_timing_plot = 'None'
sample_hist_key$WGD_timing_plot = factor(sample_hist_key$WGD_timing_plot, levels = rev(sort(unique(sample_hist_key$WGD_timing_plot))))

#rank by proportion of samples that undergo WGD and then arrange by the estimated time at which it occurs
ranked_histo = as.data.frame(unique(sample_hist_key$histology_abbreviation))
names(ranked_histo) = 'Dx'

ranked_histo$noWGD_count = NA
ranked_histo$total = NA

for(i in 1:nrow(ranked_histo)){
  ranked_histo$noWGD_count[i] = nrow(sample_hist_key[sample_hist_key$histology_abbreviation == ranked_histo$Dx[i] & sample_hist_key$WGD == 'No',])
  ranked_histo$total[i] = nrow(sample_hist_key[sample_hist_key$histology_abbreviation == ranked_histo$Dx[i],])
}

ranked_histo$prop_noWGD = ranked_histo$noWGD_count / ranked_histo$total
sample_hist_key$histology_abbreviation = factor(sample_hist_key$histology_abbreviation, levels = ranked_histo[order(ranked_histo$prop_noWGD, decreasing = F),]$Dx)

#g <- colorRampPalette(RColorBrewer::brewer.pal(4,"Set1")[c(3,2,4)])(101)
g <- colorRampPalette(c('black', 'grey95'))(101)
names(g) = seq(0.0, 1, 0.01)
h = c(g, 'white')
names(h)[102] = 'None'

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/figures/wgd_gct_v_pcawg_plot_20210919.pdf', height = 8, width = 6)
ggplot(sample_hist_key, aes(y = histology_abbreviation, fill = WGD_timing_plot)) +
  geom_bar(position="fill") + 
  scale_fill_manual(values = h) + 
  coord_cartesian(expand = F) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('Proportion undergone WGD')  + ylab('')
dev.off()

#7) ESTIMATE NUMBER OF CELL DIVISIONS PRIOR TO WGD

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/mutationtimeR_new_swater/wkdir')

predup_burden_df = read.table('gct_predup_sub_burden_estimate_20210918.txt', header = T, sep = '\t', stringsAsFactors = F)

predup_burden_df_per_patient = aggregate(adj_genome_sub_burden ~ patient, predup_burden_df, mean)
pt_ages = manifest[!duplicated(manifest$Case.ID), c('Case.ID', 'Age', 'Diagnosis', 'Organ')]

predup_burden_df_per_patient$age = NA
predup_burden_df_per_patient$dx = NA
predup_burden_df_per_patient$site = NA
for(i in 1:nrow(predup_burden_df_per_patient)){
  predup_burden_df_per_patient$age[i] = pt_ages[pt_ages$Case.ID == predup_burden_df_per_patient$patient[i],]$Age
  predup_burden_df_per_patient$dx[i] = pt_ages[pt_ages$Case.ID == predup_burden_df_per_patient$patient[i],]$Diagnosis
  predup_burden_df_per_patient$site[i] = pt_ages[pt_ages$Case.ID == predup_burden_df_per_patient$patient[i],]$Organ
}

#Rahbari et al 2016 report 0.5-0.7 mutations per haploid genome per cell division. We can, therefore, extrapolate to 1-1.4 mutations per diploid genome per cell division (the assumed diploid PGC precursor to postpubertal GCTs)
#this is interpreted as 1.2 +/- 0.2
#assuming no uncertainty for the pre-duplication substitution burden
#Number of cell division prior to WGD (f) = Pre-WGD sub burden (A) / Estimated subs per diploid genome per cell division post-PGC specification (B)
#sig_f = f * (sig_B / B) where B = 1.2 and sig_B = 0.2

predup_burden_df_per_patient$estimated_postPGC_cell_divisions = predup_burden_df_per_patient$adj_genome_sub_burden / 1.2 #estimated subs per diploid genome per cell division post PGC specification

#take estimate of postPGC mutation rate per cell division in testis and ovary from Rahbari et al. 2016
predup_burden_df_per_patient$estimated_postPGC_cell_divisions_err = predup_burden_df_per_patient$estimated_postPGC_cell_divisions * (0.2 / 1.2)


median(predup_burden_df_per_patient[predup_burden_df_per_patient$age > 12,]$estimated_postPGC_cell_divisions) #4.806752
range(predup_burden_df_per_patient[predup_burden_df_per_patient$age > 12,]$estimated_postPGC_cell_divisions) #0.00000 14.97846
median(predup_burden_df_per_patient[predup_burden_df_per_patient$age > 12,]$estimated_postPGC_cell_divisions_err) #0.8011253