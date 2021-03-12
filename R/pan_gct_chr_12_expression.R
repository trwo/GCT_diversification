#pan-invasive GCT chromosome 12 expression and downstream analysis

#1. perform pan-invasive GCT DE analysis

library(limma)
library(edgeR)

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/pre_processing')

#read in raw count data
raw_data = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/filtering/germ_cell_lcm_rna_preQC.txt', header = T, sep = '\t', stringsAsFactors = F) #raw count data
metadata = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/supplementary_data/GCT_RNA_metadata.txt', header = T, sep = '\t', stringsAsFactors = F) #read in metadata for filtered, final cuts
row.names(metadata) = metadata$Sample_ID
gct_mat = raw_data[, 7:ncol(raw_data)]
row.names(gct_mat) = raw_data$Geneid

#remove cuts that failed QC or are of in situ disease - we want to compare invasive GCT tissues vs normal seminiferous tubules only
gct_mat = gct_mat[, colnames(gct_mat) %in% metadata[metadata$Updated.Description != 'GCNIS',]$Sample_ID] #382 samples after removing GCNIS

#define factors that will be input into the model
histology = factor(metadata[colnames(gct_mat), ]$Updated.Description)
tumour = factor(metadata[colnames(gct_mat), ]$Normal.Tumour)
metadata = metadata[colnames(gct_mat),]

#preprocessing
y <- DGEList(counts = gct_mat, group = tumour)
keep <- filterByExpr(y, group = tumour)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y) #12431 genes kept
y <- estimateCommonDisp(y) #to generate pseudocounts, adjusted for library size

#visualise if microbiopsies from the same histology cluster and if there is patient effect
plotMDS(y, labels = tumour, col = as.numeric(histology))

#specify model to be fitted, no intercept to be fitted, each normal-tumour status will have its own group mean
design = model.matrix(~ 0 + tumour)
colnames(design) = levels(tumour)

#run double voom to optimise regression weights, adjust for histology as 'block' to find underlying generic tumour signal
v <- voom(y, design)
corfit <- duplicateCorrelation(v, design, block = histology)
v <- voom(y, design, block = histology, correlation = corfit$consensus)
corfit <- duplicateCorrelation(v, design, block = histology)
fit <- lmFit(v, design, block = histology, correlation = corfit$consensus)
corfit$consensus #0.1980107

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/limma_voom_de')

saveRDS(v, 'LCM_invasive_GCT_v_testis_DE_voom_object_20210305.RDS')
saveRDS(fit, 'LCM_invasive_GCT_v_testis_DE_model_fit_20210305.RDS')

v = readRDS('LCM_invasive_GCT_v_testis_DE_voom_object_20210305.RDS')
fit = readRDS('LCM_invasive_GCT_v_testis_DE_model_fit_20210305.RDS')

cont_matrix = makeContrasts(Tumour - Normal, levels = colnames(coef(fit)))
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)
top.table <- topTable(fit2, number = nrow(fit$coefficients), sort.by = "P")

raw_data = read.table("/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/filtering/germ_cell_lcm_rna_preQC.txt", header = T, sep = '\t', stringsAsFactors = F, check.names = T)
row.names(raw_data) = raw_data$Geneid
raw_data$Start = unlist(lapply(strsplit(raw_data$Start, ';'), '[[', 1))

top.table$chr = raw_data[row.names(top.table),]$Chr
top.table$start = raw_data[row.names(top.table),]$Start
top.table$start = as.numeric(top.table$start)

write.table(top.table, 'invasive_gct_v_normtestis_voom_de_20210304.txt', col.names = T, sep = '\t', row.names = T, quote = F)

#2. Plot rolling average logFC in gene expression across chromosome 12, with the (copy number - ploidy) overlaid

#fetch cn data from Fig. 1
cn_df = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/combined_CN_RNA_exp/lcm_gct_cn_10000_bp_bin_per_GCT_WG_median_extra_cols.txt', header = T, sep = '\t', stringsAsFactors = F)

#limit to chromosome 12 (cumulative coordinates from 10kb bins in file, hence odd coordinate values)
lcm_chr12_df = cn_df[(cn_df$start > 1950970000) & (cn_df$end <= 2084830000), ]
lcm_chr12_df$chr = 'chr12'
lcm_chr12_df$start = lcm_chr12_df$start - 1950970000
lcm_chr12_df$end = lcm_chr12_df$end - 1950970000
lcm_chr12_df = lcm_chr12_df[!is.na(lcm_chr12_df$all_samples),]

#fetch logFC data from invasive v normal comparison
top.table.lcm = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/limma_voom_de/invasive_gct_v_normtestis_voom_de_20210304.txt', header = T, sep = '\t', stringsAsFactors = F)

#limit to chromsome 12 data for plotting
lcm_chr12_exp_data = top.table.lcm[top.table.lcm$chr == 12,]
lcm_chr12_exp_data$chr = paste0('chr', lcm_chr12_exp_data$chr)
lcm_chr12_exp_data = lcm_chr12_exp_data[order(lcm_chr12_exp_data$start),]

library(karyoploteR)
library(zoo)

#generic plotting settings
chromosomes = 'chr12'
dna_r0 <- 0.9
dna_r1 <- 0
rna_r0 <- 0
rna_r1 <- 0.9
karyoplotr_plottype = 1

#in-house rna.data, creating rolling mean
window_size = 50 #must be an even number
y_data = rollmean(lcm_chr12_exp_data$logFC, k = window_size, align = 'center')
x_data = lcm_chr12_exp_data$start[(window_size / 2):(length(lcm_chr12_exp_data$start) - (window_size / 2))]
pos_y_data = ifelse(y_data > 0, y_data, 0) 
neg_y_data = ifelse(y_data < 0, y_data, 0)

#plot in-house data
title = 'In-house'

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/gct_cn_v_exp_karyoplotter_20210305.pdf', height = 6, width = 6, useDingbats = F)
kp <- karyoploteR::plotKaryotype(genome = "hg19", chromosomes = 'chr12', plot.type = karyoplotr_plottype, cex = 0.85, main = gsub(".*__", "", title))
karyoploteR::kpDataBackground(kp, r0 = dna_r0, r1 = dna_r1, col = 'grey90')
karyoploteR::kpAxis(karyoplot = kp, col = "black", ymin = -1.5, ymax = 1.5, r0 = rna_r0, r1 = rna_r1, cex = 1, tick.pos = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), side = 2)
karyoploteR::kpBars(kp, chr = 'chr12', x0 = x_data - 100000, x1 = x_data + 100000, y0 = 0, y1 = pos_y_data, col = '#95D840', data.panel = 1, ymin = -1.5, ymax = 1.5, r0 = rna_r0, r1 = rna_r1, clipping = T, border = 'NA')
karyoploteR::kpBars(kp, chr = 'chr12', x0 = x_data - 100000, x1 = x_data + 100000, y0 = 0, y1 = neg_y_data, col = '#A020F0', data.panel = 1, ymin = -1.5, ymax = 1.5, r0 = rna_r0, r1 = rna_r1, clipping = T, border = 'NA')
karyoploteR::kpAxis(karyoplot = kp, col = "black", ymin = 6, ymax = -3, r1 = dna_r1, r0 = dna_r0, cex = 1, tick.pos = c(-3, 0, 3, 6), side = 1)
karyoploteR::kpPoints(kp, chr = lcm_chr12_df$chr, x = lcm_chr12_df$start, y = lcm_chr12_df$all_samples, col = "black", ymin = 6, ymax = -3, r0 = dna_r0, r1 = dna_r1, clipping = T) #reintroduced prepubertal samples
dev.off()

#3. Is the 12p overexpression due to chance, compared to regions with copy number states nearer baseline ploidy?

#fetch cn data from Fig. 1
cn_df = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/combined_CN_RNA_exp/lcm_gct_cn_10000_bp_bin_per_GCT_WG_median_extra_cols.txt', header = T, sep = '\t', stringsAsFactors = F)

#limit to chromosome 12 (cumulative coordinates from 10kb bins in file, hence strange-looking coordinates)
lcm_chr12_df = cn_df[(cn_df$start > 1950970000) & (cn_df$end <= 2084830000), ]
lcm_chr12_df$chr = 'chr12'
lcm_chr12_df$start = lcm_chr12_df$start - 1950970000
lcm_chr12_df$end = lcm_chr12_df$end - 1950970000
lcm_chr12_df = lcm_chr12_df[!is.na(lcm_chr12_df$all_samples),]

#fetch logFC data from invasive v normal comparison

top.table.lcm = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/limma_voom_de/invasive_gct_v_normtestis_voom_de_20210304.txt', header = T, sep = '\t', stringsAsFactors = F)
top.table.lcm = top.table.lcm[top.table.lcm$chr %in% c(1:22),] #only performed copy number plotting across autosomal genome for comparison

#limit to chromsome 12 data for plotting
lcm_chr12_exp_data = top.table.lcm[top.table.lcm$chr == 12,]
lcm_chr12_exp_data$chr = paste0('chr', lcm_chr12_exp_data$chr)
lcm_chr12_exp_data = lcm_chr12_exp_data[order(lcm_chr12_exp_data$start),]

#get regions with a CN-ploidy within 0.5 of the overall ploidy located outside of chr12
av_cn = cn_df[!(is.na(cn_df$all_samples)) & (cn_df$all_samples > -0.5) & (cn_df$all_samples < 0.5) & ((cn_df$start < 1950970000) | (cn_df$start > 2084830000)),]

#provided cumulative loci for genes
top.table.lcm$mod_start = top.table.lcm$start
for(k in 1:nrow(top.table.lcm)) {
  if(top.table.lcm$chr[k] != "1"){
    top.table.lcm$mod_start[k] = top.table.lcm$start[k] + chrom_sizes$prior_wg_length[which(chrom_sizes$chr == top.table.lcm$chr[k])]
  }
}

top.table.lcm$av_cn_region = NA

for(i in 1:nrow(top.table.lcm)){
  if(nrow(av_cn[(av_cn$start < top.table.lcm$mod_start[i]) & (av_cn$end > top.table.lcm$mod_start[i]), ]) > 0){
    top.table.lcm$av_cn_region[i] = 'Yes'
  }
  if(nrow(av_cn[(av_cn$start < top.table.lcm$mod_start[i]) & (av_cn$end > top.table.lcm$mod_start[i]), ]) == 0){
    top.table.lcm$av_cn_region[i] = 'No'
  }
}

n_12p_genes = nrow(top.table.lcm[top.table.lcm$chr == 12 & top.table.lcm$start < 35800000,]) #226 used
mean_12p = mean(top.table.lcm[top.table.lcm$chr == 12 & top.table.lcm$start < 35800000,]$logFC) # 0.686984

nRand = 100000
randScores=c()
set.seed(42)
av_cn_exp_regions = top.table.lcm[top.table.lcm$av_cn_region == 'Yes',]

#average across 226 genes from other genomic regions near baseline ploidy
for(i in 1:nRand){
  randScores=c(randScores, mean(sample(av_cn_exp_regions$logFC, n_12p_genes, replace = F)))
}

pVal = sum(randScores >= mean_12p)/nRand #0 random samples have a higher mean logFC than that seen across 12p

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/rand_sampling_gene_exp_v_12p_20210305.pdf', height = 3, width = 4)
ggplot() + geom_density(mapping = aes(x = randScores), fill = '#999999') + xlim(-2, 1) + theme_pubr() + geom_vline(xintercept = mean_12p) + labs(x = 'log2FC of gene expression in regions with baseline ploidy') + coord_cartesian(expand = F)
dev.off()
