
setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/original_nature_submission/analyses/combined_CN_RNA_exp/')

library(plyr)

#get metadata, including ploidy values for all eligible tumours
tcga_ploidy = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/reference_datasets/tcga_tgct/summary.ascatTCGA.penalty70.hg19.tsv', header = T, sep = '\t', stringsAsFactors = F)
tcga_ploidy = tcga_ploidy[(tcga_ploidy$cancer_type == 'TGCT') & (tcga_ploidy$pass == TRUE) & (tcga_ploidy$rep == TRUE) & (tcga_ploidy$purity > 0.4), ] #103 samples

#get copy segments for all tumours in a single dataframe
tcga_segment_files = list.files(path = '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/original_nature_submission/reference_datasets/tcga_tgct', pattern = 'segments.txt')
paste0(tcga_ploidy$name, '.segments.txt') %in% tcga_segment_files #all present in segment files

tcga_ascat_segments = c()
for(i in tcga_ploidy$name){
  tumour_file = read.table(paste0('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/original_nature_submission/reference_datasets/tcga_tgct/', i, '.segments.txt'), header = T, sep = '\t', stringsAsFactors = F)
  tcga_ascat_segments = rbind(tcga_ascat_segments, tumour_file)
}

names(tcga_ascat_segments) = c('tumour', 'chr', 'start', 'end', 'majA', 'minA')
tcga_ascat_segments$ntot = tcga_ascat_segments$majA + tcga_ascat_segments$minA

samples <- unique(tcga_ascat_segments$tumour)

#set bin size
bin_size <- 10000
filter_size = 10000 #remove CN changes smaller than this
adjust = bin_size - 1 #for arranging segments

#chromosome lengths, each will be rounded up to the nearest multiple of 10,000
chrom_sizes <- read.table('/lustre/scratch119/casm/team294rr/rs30/Testicular_project/Final_Analyses/metadata/male.hg19.chrom.sizes.tsv', header = F, sep = '\t', stringsAsFactors = F)
chrom_sizes[,2] = round_any(chrom_sizes[,2], bin_size, ceiling) #chunk each chromosome into 10kb bins in order to accurate characterise the CN at the beginning and ends of chromosomes, has to be ceiling to capture the ends
chrom_sizes[,3] <- cumsum(as.numeric(chrom_sizes[,2]))
chrom_sizes = chrom_sizes[!chrom_sizes$V1 %in% c('chrM', 'chrX', 'chrY'),]
chrom_sizes$V4 = c(0, chrom_sizes$V3[1:(length(chrom_sizes$V3) - 1)])
chrom_sizes$V1 = unlist(lapply(strsplit(chrom_sizes$V1, 'chr'), '[[', 2))
names(chrom_sizes) = c('chr', 'chr_length', 'cumsum_length', 'prior_wg_length')

tcga_ascat_segments = tcga_ascat_segments[(tcga_ascat_segments$end - tcga_ascat_segments$start) > filter_size,]
tcga_ascat_segments$mod_start = round_any(tcga_ascat_segments$start, bin_size, round) #fit genome to bins of 10,000bp
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

write.table(tcga_cn_df, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/original_nature_submission/analyses/combined_CN_RNA_exp/tcga_cn_10000_bp_GCT_WG_median.txt', col.names = T, row.names = F, sep = '\t', quote = F) 