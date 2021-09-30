setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/combined_CN_RNA_exp/')

library(plyr)

#locate copy number files and the eligible samples
project_dir <- "/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/variant_calls/battenberg/00_final_calls/long_format_bb"
dna_cuts = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/supplementary_data/WGS_summary_metadata.txt', header = T, sep = '\t', stringsAsFactors = F)
pure_dna_cuts = dna_cuts[dna_cuts$Purity > 0.4 & dna_cuts$Histology != 'GCNIS' & dna_cuts$Experimental_arm == 'LCM',] #pure, invasive only
samples <- pure_dna_cuts$Sample #all samples from all donors
donors <- unique(substr(samples,1,7)) #all donors

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

#create array with per donor and per 10kb bin cells
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
    if(file.exists(paste0(project_dir,"/",donor_samples[j],"_subclones_segments.txt"))){
      print(donor_samples[j])
      segments <- read.table(paste0(project_dir,"/",donor_samples[j],"_subclones_segments.txt"), sep = "\t", stringsAsFactors = F, header = T)
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
      for(m in 1:nrow(donor_sample_df)){ #get copy number estimates per 10000bp bin
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

write.table(cn_df, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/combined_CN_RNA_exp/lcm_gct_cn_10000_bp_bin_per_GCT_WG_median.txt', col.names = T, row.names = F, sep = '\t', quote = F) 