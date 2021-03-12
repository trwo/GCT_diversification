#WGD estimation analysis
#Code for analysis adpated from work by Moritz Gerstung and Santiago Gonzalez

#1. Functions

#Import some functions from PCAWG paper
#source('/nfs/users/nfs_t/to3/scripts/additional_mutationtimeR_functions.R') #source extra, unchanged functions

#edited functions from mutationtimeR
findMainCluster <- function(bb, min.dist = 0.05){ #at least 20 SNVs
  w <- which(bb$n.snv_mnv >= 20 & !is.na(bb$time) & ((bb$time.up - bb$time.lo) < 0.5)) #edited to exclude regions with very wide confidence intervals when, to3 03/03/21 estimating WGD time
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
  #note that this line has been left as is in fractionGenomeWgdCompatible as there we only want to know the regions compatable with WGD, not necessarily depend on them all for downstream analyses, , to3 03/03/21
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

#2. Input

#Create lists of pertinent input data, mostly dervied from mutationtimeR output
setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/mutationtimeR')

chrom_sizes <- read.table('/lustre/scratch119/casm/team294rr/rs30/Testicular_project/Final_Analyses/metadata/male.hg19.chrom.sizes.tsv', header = F, sep = '\t', stringsAsFactors = F)
genome_length = sum(chrom_sizes[1:24,]$V2) #3095677412

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

purity_files = list.files(pattern = '.battenberg.rho_and_psi_updated.txt')

finalPurity <- numeric()
for(f in purity_files){
  pur_df = read.table(f, sep = '\t', header = T, stringsAsFactors = F)
  purity = pur_df['FRAC_GENOME', 'rho']
  finalPurity[unlist(strsplit(f, '.battenberg.rho_and_psi_updated.txt'))] = purity
}

sex_vec = ifelse(grepl('PD455', purity_files), 'female', 'male')
names(sex_vec) = unlist(strsplit(purity_files, '.battenberg.rho_and_psi_updated.txt'))

#3. WGD timing in mutation time

#Sanity check for final, eligible samples - are they consistent with WGD?

finalPloidy <- sapply(finalBB, averagePloidy)
names(finalPloidy) <- names(finalBB)

finalHom <- sapply(finalBB, averageHom)
names(finalHom) <- names(finalBB)

isWgd <- .classWgd(finalPloidy, finalHom)
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
#  PD42036   0.0220
#  PD43296   0.0155
#  PD43298   0.0570
#  PD43299   0.3935
#  PD45543   0.0455
#  PD45544   0.0290
#  PD45545   0.0195
#  PD46269   0.0100
#  PD46271   0.0175
#  PD46966   0.0380
#  PD46969   0.0065

#4. Estimate pre-duplication substitution burden
wgdParam <- mclapply(names(finalSnv)[isWgd], function(ID){
  try(computeWgdParam(finalSnv[[ID]], finalBB[[ID]], clusters=finalClusters[[ID]], purity=finalPurity[ID], sex=sex_vec[ID]))
})
names(wgdParam) <- names(finalSnv)[isWgd]
void <- sapply(wgdParam, function(x) is.null(x) | class(x)=="try-error") #none

saveRDS(wgdParam, 'mutationtimeR_WGD_segments_20210303.rds')

predup_burden_df = data.frame(sample = names(wgdParam))
row.names(predup_burden_df) = predup_burden_df$sample

predup_burden_df$adj_genome_sub_burden = NA
predup_burden_df$hq_seg_width = NA

for(i in 1:nrow(predup_burden_df)){
  predup_burden_df$adj_genome_sub_burden[i] = wgdParam[[predup_burden_df$sample[i]]]$adj_genome_sub_burden
  predup_burden_df$hq_seg_width[i] = wgdParam[[predup_burden_df$sample[i]]]$hq_seg_width
}

predup_burden_df$patient = substr(predup_burden_df$sample, 0, 7)
predup_burden_df_per_patient = aggregate(adj_genome_sub_burden ~ patient, predup_burden_df, median)

#   patient adj_genome_sub_burden
#  PD42036              0.000000
#  PD43296              7.151071
#  PD43298             13.736280
#  PD43299            419.612430
#  PD45543              0.000000
#  PD45544              0.000000
#  PD45545              7.617121
#  PD46269              0.000000
#  PD46271              5.903524
#  PD46966              7.210443
#  PD46969             15.341817

median(predup_burden_df_per_patient[predup_burden_df_per_patient$patient != 'PD43299',]$adj_genome_sub_burden) #6.527298

write.table(fracGenomeWgdComp, 'gct_wgd_estimate_mut_time_20210304.txt', col.names = T, row.names = T, sep = '\t', quote = F)
write.table(predup_burden_df, 'gct_predup_sub_burden_estimate_20210304.txt', col.names = T, row.names = F, sep = '\t', quote = F)

fracGenomeWgdComp = read.table('gct_wgd_estimate_mut_time_20210304.txt', header = T, sep = '\t', stringsAsFactors = F)

library(reshape2)
library(ggplot2)
library(ggpubr)
library(VariantAnnotation)
library(ggbeeswarm)

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/wgd_timing_per_microbiopsy_20210309.pdf', useDingbats = F, height = 6, width = 6)
ggplot(fracGenomeWgdComp) +
  geom_quasirandom(mapping = aes(x = patient, y = time.wgd), col = '#999999') +
  geom_pointrange(mapping = aes(x = patient, y = time.wgd),
                  stat = "summary", shape=19, 
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill="black") +
  theme_pubr() + ylab('Estimated WGD timing') + xlab('') +
  scale_x_discrete(expand = c(0,0.5)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0.01)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.9))
dev.off()

table(fracGenomeWgdComp$patient)
#PD42036 PD43296 PD43298 PD43299 PD45543 PD45544 PD45545 PD46269 PD46271 PD46966 PD46969 
#      5      22      20       4       2       5       2       3       2       3       4 

#5. Plot pre-WGD sub burden against truncal mutation burden

library(reshape2)
library(ggplot2)
library(ggpubr)
library(VariantAnnotation)
library(ggbeeswarm)

sub_prop = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/supplementary_data/gct_trunk_v_overall_subs.txt', header = T, sep = '\t', stringsAsFactors = F, row.names = 1)
sub_prop$patient = row.names(sub_prop)

sub_prop.m = reshape2::melt(sub_prop[, c('Truncal.Overall', 'X1...Truncal.Overall.', 'patient')], id.vars = 'patient')
sub_prop.m$variable = factor(sub_prop.m$variable, levels = c('X1...Truncal.Overall.', 'Truncal.Overall'))
sub_prop.m$patient = factor(sub_prop.m$patient, levels = c("PD46969", "PD46966", "PD46271", "PD46270", "PD46269", "PD45545", "PD45544", "PD45543", "PD43298", "PD43296", "PD42036", "PD43299"))

timing_tab = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/mutationtimeR/gct_predup_sub_burden_estimate_20210304.txt', header = T, sep = '\t', stringsAsFactors = F)

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/pre_dup_burden_per_microbiopsy_20210304.pdf', useDingbats = F, height = 6, width = 6)
ggplot(timing_tab) +
  geom_quasirandom(mapping = aes(x = patient, y = adj_genome_sub_burden), col = '#999999') +
  geom_pointrange(mapping = aes(x = patient, y = adj_genome_sub_burden),
                  stat = "summary", shape=19, 
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill="black") +
  theme_pubr() + ylab('Estimated pre-duplication substution burden') + xlab('') +
  scale_x_discrete(expand = c(0,0.5)) +
  scale_y_continuous(limits = c(0, 500), expand = c(0, 0.1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.9))
dev.off()

final_wgd_plotting_agg = aggregate(adj_genome_sub_burden ~ patient, timing_tab, median)

final_wgd_plotting_agg$trunc_mt = NA
for(i in 1:nrow(final_wgd_plotting_agg)){
  final_wgd_plotting_agg$trunc_mt[i] = sub_prop[final_wgd_plotting_agg$patient[i], ]$Invasive.truncal.sum
}

median(final_wgd_plotting_agg[final_wgd_plotting_agg$patient != 'PD43299',]$adj_genome_sub_burden) #6.527298
final_wgd_plotting_agg$pre_WGD_burden_log = log(final_wgd_plotting_agg$adj_genome_sub_burden + 1) #to keep all values over 0
final_wgd_plotting_agg$trunc_mt_log = log(final_wgd_plotting_agg$trunc_mt + 1)

final_wgd_plotting_agg.m = reshape2::melt(final_wgd_plotting_agg[, c('patient', 'pre_WGD_burden_log', 'trunc_mt_log')], id.vars = 'patient')

final_wgd_plotting_agg.m$Patient = factor(final_wgd_plotting_agg.m$patient, levels = c("PD46969", "PD46966", "PD46271", "PD46270", "PD46269", "PD45545", "PD45544", "PD45543", "PD43298", "PD43296", "PD42036", "PD43299"))
final_wgd_plotting_agg.m$variable = factor(final_wgd_plotting_agg.m$variable, levels = c('trunc_mt_log', 'pre_WGD_burden_log'))

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/gct_sub_prop_predup_v_trunc_mt_20210304.pdf', height = 4.5, width = 3, useDingbats = F)
ggplot(final_wgd_plotting_agg.m) +
  geom_bar(mapping = aes(x = value, y = Patient, fill = variable), position="dodge", stat="identity") + scale_fill_manual(values = c('grey80', 'grey40')) + coord_cartesian(expand = F, xlim = c(log(0 + 1), log(10000 + 1))) + ylab('') + xlab('Number of substitutions') + theme_pubr(legend = 'none') + scale_x_continuous(breaks = c(log(0 + 1), log(10 + 1), log(100 + 1), log(1000 + 1), log(10000 + 1)), labels = c(0, 10, 100, 1000, 10000)) + geom_vline(xintercept = median(final_wgd_plotting_agg.m[final_wgd_plotting_agg.m$variable == 'pre_WGD_burden_log' & final_wgd_plotting_agg.m$patient != 'PD43299', ]$value), lty = 2) #2.015091, actual value underlying this is 6.527298
dev.off()

#6. Estimate the number of cell divisions prior to WGD

final_wgd_plotting_agg$estimated_postPGC_cell_divisions = final_wgd_plotting_agg$adj_genome_sub_burden / 0.6 #middle of range for estimated subs per cell division

#take estimate of postPGC mutation rate per cell division in testis and ovary from Rahbari et al. 2016
final_wgd_plotting_agg$estimated_postPGC_cell_divisions_upper = final_wgd_plotting_agg$adj_genome_sub_burden / 0.5
final_wgd_plotting_agg$estimated_postPGC_cell_divisions_lower = final_wgd_plotting_agg$adj_genome_sub_burden / 0.7

median(final_wgd_plotting_agg[final_wgd_plotting_agg$patient != 'PD43299',]$estimated_postPGC_cell_divisions) #10.87883
median(final_wgd_plotting_agg[final_wgd_plotting_agg$patient != 'PD43299',]$estimated_postPGC_cell_divisions_lower) #9.324711
median(final_wgd_plotting_agg[final_wgd_plotting_agg$patient != 'PD43299',]$estimated_postPGC_cell_divisions_upper) #13.0546

min(final_wgd_plotting_agg[final_wgd_plotting_agg$patient != 'PD43299',]$estimated_postPGC_cell_divisions_lower) #0
max(final_wgd_plotting_agg[final_wgd_plotting_agg$patient != 'PD43299',]$estimated_postPGC_cell_divisions_upper) #30.68363

#7. Compare WGD timing to PCAWG, data kindly provided by Peter van Loo

summary_data = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/supplementary_data/GCT_DNA_supplementary_data.txt', header = T, sep = '\t', stringsAsFactors = F)
row.names(summary_data) = summary_data$Sample

#read in PCAWG data
pcawg_wgd = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/reference_datasets/20190824_timing_gains_adj.txt', header = T, sep = '\t', stringsAsFactors = F)

pcawg_wgd[pcawg_wgd$samplename %in% pcawg_wgd[pcawg_wgd$type == 'WGD',]$samplename,]


sample_hist_key = pcawg_wgd[!duplicated(pcawg_wgd[, c('samplename', 'histology_abbreviation')]), c('samplename', 'histology_abbreviation')]
pcawg_wgd = pcawg_wgd[pcawg_wgd$histology_abbreviation %in% sort(unique(sample_hist_key$histology_abbreviation))[table(sample_hist_key$histology_abbreviation) >= 10], ] #only keep tumour types with at least 10 samples to make proportion values sensible

sample_hist_key = sample_hist_key[sample_hist_key$histology_abbreviation %in% sort(unique(sample_hist_key$histology_abbreviation))[table(sample_hist_key$histology_abbreviation) >= 10], ]

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
timing_df = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/mutationtimeR/gct_wgd_estimate_mut_time_20210304.txt', sep = '\t', header = T, stringsAsFactors = F)

median_wgd = aggregate(time.wgd ~ patient, timing_df, median)
median_wgd_postpubertal = median_wgd[median_wgd$patient != 'PD43299',]
median_wgd_postpubertal$histology_abbreviation = 'Postpubertal_GCT'

gct_samples = data.frame(samplename = median_wgd_postpubertal$patient, histology_abbreviation = rep('Postpubertal_GCT', nrow(median_wgd_postpubertal)), WGD = rep('Yes', nrow(median_wgd_postpubertal)), WGD_timing = median_wgd_postpubertal$time.wgd)

median(median_wgd_postpubertal[median_wgd_postpubertal$patient %in% c('PD45543', 'PD45544', 'PD45545'),]$time.wgd) #post-pubertal ovary 0.029
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
sample_hist_key$histology_abbreviation = factor(sample_hist_key$histology_abbreviation, levels = ranked_histo[order(ranked_histo$prop_noWGD, decreasing = T),]$Dx)

g <- colorRampPalette(c('grey80', 'black'))(101)
names(g) = seq(0.0, 1, 0.01)
h = c(g, 'white')
names(h)[102] = 'None'

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/wgd_gct_v_pcawg_plot_20210304.pdf', height = 6, width = 8)
ggplot(sample_hist_key, aes(x = histology_abbreviation, fill = WGD_timing_plot)) +
  geom_bar(position="fill") + 
  scale_fill_manual(values = h) + 
  coord_cartesian(expand = F) +
  theme_classic() +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2)) +
  xlab('')  + ylab('Proportion undergone WGD')
dev.off()