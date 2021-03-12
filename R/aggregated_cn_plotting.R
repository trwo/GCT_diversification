#Aggregated copy number plotting - study cohort vs TCGA

#Study cohort
##Turn our battenberg files into a long segment format, credit to Stefan Dentro for the conversion script
##cd lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/final_variant_calls/CNA/
  
##for f in $(ls *.battenberg.subclones.txt);
##do
##filename=$(basename -- "$f")
##sampleID="${filename%.battenberg.subclones.txt}"
##python /nfs/users/nfs_t/to3/scripts/copynumber.py -i $sampleID.battenberg.subclones.txt -c $sampleID.battenberg.rho_and_psi_updated.txt -o . -r
##done

##submit as job to the cluster
##cd /lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/combined_CN_RNA_exp
##bsub -G team294-vwork -o $PWD/log.%J -e $PWD/err.%J -q yesterday -R'select[mem>30000] rusage[mem=30000]' -M30000 -J 'gct_cn' /software/R-3.6.1/bin/Rscript lcm_gct_WG_CN.R

##########
#lcm_gct_WG_CN.R
#in-house LCM data
setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/combined_CN_RNA_exp/')

library(plyr)

#locate copy number files and the eligible samples
project_dir <- "/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/final_variant_calls/CNA/"
dna_cuts = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/supplementary_data/GCT_DNA_supplementary_data.txt', header = T, sep = '\t', stringsAsFactors = F)
pure_dna_cuts = dna_cuts[dna_cuts$purity > 0.4 & !is.na(dna_cuts$purity) & dna_cuts$Histology != 'GCNIS',] #pure, invasive only
samples <- pure_dna_cuts$Sample #all eligible samples from all donors
donors <- unique(substr(samples,1,7)) #all donors

#set bin size
bin_size <- 10000
filter_size = 10000 #remove CN changes smaller than this

#chromosome lengths, each will be rounded up to the nearest multiple of the bin size
chrom_sizes <- read.table('/lustre/scratch119/casm/team294rr/rs30/Testicular_project/Final_Analyses/metadata/male.hg19.chrom.sizes.tsv', header = F, sep = '\t', stringsAsFactors = F)
chrom_sizes[,2] = round_any(chrom_sizes[,2], bin_size, ceiling) #chunk each chromosome into 10kb bins in order to accurate characterise the CN at the beginning and ends of chromosomes, has to be ceiling to capture the ends
chrom_sizes[,3] <- cumsum(as.numeric(chrom_sizes[,2]))
chrom_sizes = chrom_sizes[!chrom_sizes$V1 %in% c('chrM', 'chrX', 'chrY'),]
chrom_sizes$V4 = c(0, chrom_sizes$V3[1:(length(chrom_sizes$V3) - 1)])
chrom_sizes$V1 = unlist(lapply(strsplit(chrom_sizes$V1, 'chr'), '[[', 2))
names(chrom_sizes) = c('chr', 'chr_length', 'cumsum_length', 'prior_wg_length')

#create array with per donor and per bin
cn_arr <- array(data = NA, dim = c((chrom_sizes[22, 'cumsum_length'] / bin_size), length(donors)))
colnames(cn_arr) = donors
cn_df = data.frame(cn_arr)

for(i in 1:length(donors)) {
  print(donors[i])
  donor_samples = samples[grepl(donors[i], samples)]
  donor_sample_arr = array(data = NA, dim = c((chrom_sizes[22, 'cumsum_length'] / bin_size), length(donor_samples))) #need copy number per site per sample first before averaging across donor and ultimately across all tumours
  colnames(donor_sample_arr) = donor_samples
  donor_sample_df = data.frame(donor_sample_arr)
  for(j in 1:length(donor_samples)){
    if(file.exists(paste0(project_dir,"/",donor_samples[j],"_segments.txt"))){
      print(donor_samples[j])
      segments <- read.table(paste0(project_dir,"/",donor_samples[j],"_segments.txt"), sep = "\t", stringsAsFactors = F, header = T)
      segments = segments[segments$chromosome %in% c(1:22),] #avoid sex chromosomes due to sampling of men and women
      segments = segments[order(segments$clonal_frequency, decreasing = T),] #order by cn state with greatest clonal frequency first
      segments = segments[!duplicated(segments[,1:3]),] #only keep the integer value for the copy number state with the highest clonal frequency
      segments = segments[, c('chromosome', 'start', 'end', 'copy_number')]
      names(segments) = c('chr', 'startpos', 'endpos', 'ntot')
      
      
      segments = segments[(segments$endpos - segments$startpos) > filter_size,] #minimum copy number segment size
      segments$startpos = round_any(segments$startpos, bin_size, round) #round start and end positions to nearest 10kb
      segments$endpos = round_any(segments$endpos, bin_size, round)
      segments$startcumpos <- segments$startpos #treating genome as one continuous length, generate cumsum postions
      segments$endcumpos <- segments$endpos
      
      
      for(k in 1:nrow(segments)) { #annotate segments as cumulative positions
        segments$startcumpos[k] = segments$startpos[k] + chrom_sizes$prior_wg_length[which(chrom_sizes$chr == segments$chr[k])]
        segments$endcumpos[k] <- segments$endpos[k] + chrom_sizes$prior_wg_length[which(chrom_sizes$chr == segments$chr[k])]
      }
      
      sample_vec = c()
      for(m in 1:nrow(donor_sample_df)){ #get copy number estimates per bin
        if(nrow(segments[segments$startcumpos <= (bin_size * (m - 1)) & segments$endcumpos >= (bin_size * m),]) > 0) sample_vec = c(sample_vec, segments[segments$startcumpos <= (bin_size * (m - 1)) & segments$endcumpos >= (bin_size * m),]$ntot)
        else {
          sample_vec = c(sample_vec, NA)
        }
      }
      donor_sample_df[, donor_samples[j]] = sample_vec
    }
  }
  cn_df[, donors[i]] = apply(donor_sample_df, 1, median, na.rm = T) #median per tumour for that bin
}

write.table(cn_df, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/combined_CN_RNA_exp/lcm_gct_cn_10000_bp_bin_per_GCT_WG_median.txt', col.names = T, row.names = F, sep = '\t', quote = F) 

##########

##read in output and add further stats
setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/combined_CN_RNA_exp/')

cn_df = read.table('lcm_gct_cn_10000_bp_bin_per_GCT_WG_median.txt', header = T, sep = '\t', stringsAsFactors = F)
cn_df$start = c(1:nrow(cn_df))*10000 - 9999
cn_df$end = c(1:nrow(cn_df))*10000

##adjust for ploidy, to express copy number relative to the overall ploidy 
dna_cuts = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/supplementary_data/GCT_DNA_supplementary_data.txt', header = T, sep = '\t', stringsAsFactors = F)
pure_dna_cuts = dna_cuts[dna_cuts$purity > 0.4 & !is.na(dna_cuts$purity) & dna_cuts$Histology != 'GCNIS',] #pure, invasive only
lcm_ploidy_table = aggregate(ploidy ~ Case_ID, pure_dna_cuts, median)
row.names(lcm_ploidy_table) = lcm_ploidy_table$Case_ID
cn_df[,1:12] = sweep(cn_df[,1:12], 2, lcm_ploidy_table[colnames(cn_df[,1:12]),]$ploidy, FUN = '-', check.margin = T) #adjust the copy number calls for each tumour by their median ploidy

cn_df$all_samples = rowMeans(cn_df[,1:12], na.rm = T)
cn_df$postp_samples = rowMeans(cn_df[,c(1:3, 5:12)], na.rm = T) 
cn_df$uq_all_samples = apply(cn_df[,c(1:12)], MARGIN = 1, FUN = function(z) {quantile(z, 0.75, na.rm = T)})
cn_df$lq_all_samples = apply(cn_df[,c(1:12)], MARGIN = 1, FUN = function(z) {quantile(z, 0.25, na.rm = T)})
cn_df$uq_postp_samples = apply(cn_df[,c(1:3, 5:12)], MARGIN = 1, FUN = function(z) {quantile(z, 0.75, na.rm = T)})
cn_df$lq_postp_samples = apply(cn_df[,c(1:3, 5:12)], MARGIN = 1, FUN = function(z) {quantile(z, 0.25, na.rm = T)})

write.table(cn_df, 'lcm_gct_cn_10000_bp_bin_per_GCT_WG_median_extra_cols.txt', col.names = T, row.names = F, sep = '\t', quote = F)

#TCGA data
##Use ASCAT values, most recent run provided by Peter van Loo, limited to tumours with >40% purity
##cd /lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/combined_CN_RNA_exp
##bsub -G team294-vwork -o $PWD/log.%J -e $PWD/err.%J -q yesterday -R'select[mem>30000] rusage[mem=30000]' -M30000 -J 'tcga_cn' /software/R-3.6.1/bin/Rscript tcga_gct_WG_CN.R

##########
#tcga_gct_WG_CN.R
#TCGA reference data
setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/combined_CN_RNA_exp/')

library(plyr)

#get metadata, including ploidy values for all eligible tumours
tcga_ploidy = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/reference_datasets/tcga_tgct/summary.ascatTCGA.penalty70.hg19.tsv', header = T, sep = '\t', stringsAsFactors = F)
tcga_ploidy = tcga_ploidy[(tcga_ploidy$cancer_type == 'TGCT') & (tcga_ploidy$pass == TRUE) & (tcga_ploidy$rep == TRUE) & (tcga_ploidy$purity > 0.4), ] #103 samples

#get copy segments for all tumours in a single dataframe
tcga_segment_files = list.files(path = '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/reference_datasets/tcga_tgct', pattern = 'segments.txt')
paste0(tcga_ploidy$name, '.segments.txt') %in% tcga_segment_files #all present in segment files

tcga_ascat_segments = c()
for(i in tcga_ploidy$name){
  tumour_file = read.table(paste0('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/reference_datasets/tcga_tgct/', i, '.segments.txt'), header = T, sep = '\t', stringsAsFactors = F)
  tcga_ascat_segments = rbind(tcga_ascat_segments, tumour_file)
}

names(tcga_ascat_segments) = c('tumour', 'chr', 'start', 'end', 'majA', 'minA')
tcga_ascat_segments$ntot = tcga_ascat_segments$majA + tcga_ascat_segments$minA

samples <- unique(tcga_ascat_segments$tumour)

#set bin size
bin_size <- 10000
filter_size = 10000 #remove CN changes smaller than this

#chromosome lengths, each will be rounded up to the nearest multiple of 10,000
chrom_sizes <- read.table('/lustre/scratch119/casm/team294rr/rs30/Testicular_project/Final_Analyses/metadata/male.hg19.chrom.sizes.tsv', header = F, sep = '\t', stringsAsFactors = F)
chrom_sizes[,2] = round_any(chrom_sizes[,2], bin_size, ceiling) #chunk each chromosome into bins in order to accurate characterise the CN at the beginning and ends of chromosomes, has to be ceiling to capture the ends
chrom_sizes[,3] <- cumsum(as.numeric(chrom_sizes[,2]))
chrom_sizes = chrom_sizes[!chrom_sizes$V1 %in% c('chrM', 'chrX', 'chrY'),]
chrom_sizes$V4 = c(0, chrom_sizes$V3[1:(length(chrom_sizes$V3) - 1)])
chrom_sizes$V1 = unlist(lapply(strsplit(chrom_sizes$V1, 'chr'), '[[', 2))
names(chrom_sizes) = c('chr', 'chr_length', 'cumsum_length', 'prior_wg_length')

tcga_ascat_segments = tcga_ascat_segments[(tcga_ascat_segments$end - tcga_ascat_segments$start) > filter_size,]
tcga_ascat_segments$mod_start = round_any(tcga_ascat_segments$start, bin_size, round) #fit genome to bins
tcga_ascat_segments$mod_end = round_any(tcga_ascat_segments$end, bin_size, round)
tcga_ascat_segments = tcga_ascat_segments[tcga_ascat_segments$chr %in% c(1:22),]

tcga_ascat_segments$startcumpos = tcga_ascat_segments$mod_start
tcga_ascat_segments$endcumpos = tcga_ascat_segments$mod_end

for(j in 1:nrow(tcga_ascat_segments)) { #annotate segments as cumulative positions
  tcga_ascat_segments$startcumpos[j] = tcga_ascat_segments$mod_start[j] + chrom_sizes$prior_wg_length[which(chrom_sizes$chr == tcga_ascat_segments$chr[j])]
  tcga_ascat_segments$endcumpos[j] <- tcga_ascat_segments$mod_end[j] + chrom_sizes$prior_wg_length[which(chrom_sizes$chr == tcga_ascat_segments$chr[j])]
}

#create array per sample and per 10kb bin cells
tcga_cn_arr <- array(data = NA, dim = c((chrom_sizes[22, 'cumsum_length'] / bin_size), length(samples)))
colnames(tcga_cn_arr) = samples
tcga_cn_df = data.frame(tcga_cn_arr)
colnames(tcga_cn_df) = gsub('\\.', '-', colnames(tcga_cn_df))

for(k in 1:length(samples)){
  print(samples[k])
  for(m in 1:nrow(tcga_cn_df)){ #get copy number estimates per 10000bp bin
    if(nrow(tcga_ascat_segments[tcga_ascat_segments$startcumpos <= (10000 * (m - 1)) & tcga_ascat_segments$endcumpos >= (10000 * m) & tcga_ascat_segments$tumour == samples[k],]) > 0){
      tcga_cn_df[m, samples[k]] = tcga_ascat_segments[tcga_ascat_segments$startcumpos <= (10000 * (m - 1)) & tcga_ascat_segments$endcumpos >= (10000 * m) & tcga_ascat_segments$tumour == samples[k],]$ntot
    }
  }
}

write.table(tcga_cn_df, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/combined_CN_RNA_exp/tcga_cn_10000_bp_GCT_WG_median.txt', col.names = T, row.names = F, sep = '\t', quote = F) 
##########

##read in output and add further stats
tcga_cn_df = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/combined_CN_RNA_exp/tcga_cn_10000_bp_GCT_WG_median.txt', header = T, sep = '\t', stringsAsFactors = F, check.names = F)

#read in ploidy values
#identify samples with RNA seq data
tcga_ploidy = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/reference_datasets/tcga_tgct/summary.ascatTCGA.penalty70.hg19.tsv', header = T, sep = '\t', stringsAsFactors = F)
tcga_ploidy = tcga_ploidy[(tcga_ploidy$cancer_type == 'TGCT') & (tcga_ploidy$pass == TRUE) & (tcga_ploidy$rep == TRUE) & (tcga_ploidy$purity > 0.4), ] #103 samples
row.names(tcga_ploidy) = tcga_ploidy$name

#substract the ploidy from the counts for each tumour
tcga_cn_df = sweep(tcga_cn_df, 2, tcga_ploidy[colnames(tcga_cn_df),]$ploidy, FUN = '-', check.margin = T)
#tcga_cn_df$all_samples = apply(tcga_cn_df[, 1:103], 1, median, na.rm = T)
tcga_cn_df$all_samples = rowMeans(tcga_cn_df[, 1:103], na.rm = T)
tcga_cn_df$start = c(1:nrow(tcga_cn_df))*10000 - 9999
tcga_cn_df$end = c(1:nrow(tcga_cn_df))*10000
tcga_cn_df$uq = apply(tcga_cn_df[,1:103], MARGIN = 1, FUN = function(z) {quantile(z, 0.75, na.rm = T)})
tcga_cn_df$lq = apply(tcga_cn_df[,1:103], MARGIN = 1, FUN = function(z) {quantile(z, 0.25, na.rm = T)})

write.table(tcga_cn_df, 'tcga_cn_10000_bp_GCT_WG_median_extra_cols.txt', col.names = T, row.names = F, sep = '\t', quote = F)

##########
###PLOT###
##########

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/combined_CN_RNA_exp/')

library(ggplot2)
library(ggpubr)
library(plyr)

cn_df = read.table('lcm_gct_cn_10000_bp_bin_per_GCT_WG_median_extra_cols.txt', header = T, sep = '\t', stringsAsFactors = F, check.names = F)
tcga_cn_df = read.table('tcga_cn_10000_bp_GCT_WG_median_extra_cols.txt', header = T, sep = '\t', stringsAsFactors = F, check.names = F)

chrom_sizes <- read.table('/lustre/scratch119/casm/team294rr/rs30/Testicular_project/Final_Analyses/metadata/male.hg19.chrom.sizes.tsv', header = F, sep = '\t', stringsAsFactors = F)
chrom_sizes[,2] = round_any(chrom_sizes[,2], bin_size, ceiling) #chunk each chromosome into 10kb bins in order to accurate characterise the CN at the beginning and ends of chromosomes, has to be ceiling to capture the ends
chrom_sizes[,3] <- cumsum(as.numeric(chrom_sizes[,2]))
chrom_sizes = chrom_sizes[!chrom_sizes$V1 %in% c('chrM', 'chrX', 'chrY'),]
chrom_sizes$V4 = c(0, chrom_sizes$V3[1:(length(chrom_sizes$V3) - 1)])
chrom_sizes$V1 = unlist(lapply(strsplit(chrom_sizes$V1, 'chr'), '[[', 2))
names(chrom_sizes) = c('chr', 'chr_length', 'cumsum_length', 'prior_wg_length')

#get genes loci, examine overall hotspots
data = read.table("/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/featureCounts/lcm_rna_featurecount_data_preQC_20201013.txt", header = T, sep = '\t', stringsAsFactors = F, check.names = T)
row.names(data) = data$Geneid
data$Start = as.numeric(unlist(lapply(strsplit(data$Start, ';'), '[[', 1)))
data$End = as.numeric(unlist(lapply(strsplit(data$End, ';'), '[[', 1)))
data$startcumpos = data$Start
data$endcumpos = data$End
data = data[data$Chr %in% c(1:22),] #only autosomal chromosomes

for(k in 1:nrow(data)) { #annotate segments as cumulative positions
  if(data$Chr[k] != "1"){
    data$startcumpos[k] = data$Start[k] + chrom_sizes$prior_wg_length[which(chrom_sizes$chr == data$Chr[k])]
    data$endcumpos[k] = data$End[k] + chrom_sizes$prior_wg_length[which(chrom_sizes$chr == data$Chr[k])]
  }
}
data = data[, c('Geneid', 'Chr', 'Start', 'End', 'startcumpos', 'endcumpos')]
genes_to_plot = data[c('KRAS', 'KIT', 'CCND2'),]
genes_to_plot$plot_pos = NA
for(i in 1:nrow(genes_to_plot)){
  genes_to_plot$plot_pos[i] = cn_df[cn_df$start < genes_to_plot$startcumpos[i] & cn_df$end > genes_to_plot$startcumpos[i], ]$postp_samples
}

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/lcm_gct_v_tcga_cn_20210209.pdf', height = 3, width = 12, useDingbats = F)
ggplot() + 
  scale_x_continuous(expand = c(0,0), breaks = c(chrom_sizes$cumsum_length - (chrom_sizes$chr_length / 2)), labels = c(1:22), limits = c(0, 2881130000)) +
  geom_ribbon(tcga_cn_df, mapping = aes(x = start, ymin = lq, ymax = uq), alpha = 0.4, colour = 'transparent', fill = '#999999') +
  geom_point(cn_df, mapping = aes(x = start, y = postp_samples), size = 0.5) + 
  geom_vline(xintercept = chrom_sizes$prior_wg_length) + labs(y = 'copy number - ploidy', x = 'chr') + theme_pubr() + theme(panel.border = element_rect(fill = NA)) + coord_cartesian(expand = F) + geom_hline(yintercept = 0, col = '#999999', lty = 2) + ylim(-2, 8) + geom_text(genes_to_plot, mapping = aes(x = startcumpos, y = plot_pos + 0.5), label = c('KRAS', 'KIT', 'CCND2'))
dev.off()
