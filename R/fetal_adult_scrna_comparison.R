#Comparing expression of fetal and lineage-specific genes between GCTs and reference single cell data
##Using TPM values

#1. Readd in libraries
library(ggplot2)
library(reshape2)
library(cowplot)
library(Seurat)
library(ggbeeswarm)
library(dplyr)
library(ggpubr)
library(gridExtra)

#2. Assemble reference data set across genes that intersect those used in our cohort
setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/reference_datasets')

#read in scRNA seurat objects containing the count data
testis_ref = readRDS('postnatal_testis_final.rds') #neonatal and adult testis - Sohni et al.
fetal_brain = readRDS('pfc_final.rds') #fetal brain tissues - Zhong et al.
adult_brain = readRDS('adult_brain_li.rds') #adult brain - Li et al.
fetal_muscle = readRDS('fetal_muscle_cao.rds') #fetal cardiac, smooth and skeletal muscle - Cao et al.
adult_smooth_musc = readRDS('adult_smooth_muscle.rds') #adult cardiac, smooth and skeletal muscle - Litviňuková et al.

#read in our TPM data
gct_tpm = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/pre_processing/tmm_norm_gct_tpm_20210304.txt', header = T, sep = '\t', stringsAsFactors = F)

#merge reference scRNA data and restrict to genes also found in our data
ref_scrna = merge(testis_ref, y = c(fetal_brain, adult_brain, fetal_muscle, adult_smooth_musc), add.cell.ids = c("testis", "fetal_brain", "adult_brain","fetal_smooth_muscle", "adult_smooth_muscle"), project = 'ref_scrna_data') #49128 features across 32036 samples within 1 assay 

genes_to_keep = row.names(ref_scrna)[row.names(ref_scrna) %in% row.names(pseudo_counts)]

ref_scrna_final = CreateSeuratObject(counts = ref_scrna@assays$RNA@counts[genes_to_keep,], project = 'final_ref_scrna', assay = 'RNA', meta.data = ref_scrna@meta.data) #27569 features across 46195 samples within 1 assay

ref_scrna_final = subset(ref_scrna_final, subset = cell_types != 'Microglia' & cell_types != 'Pre-Spermatogonia') #27168 features across 32036 samples within 1 assay

ref_scrna_final = NormalizeData(object = ref_scrna_final, normalization.method = "LogNormalize")
all.genes <- rownames(ref_scrna_final)
ref_scrna_final <- ScaleData(ref_scrna_final, features = all.genes)
ref_scrna_final@meta.data$cell_types = factor(ref_scrna_final@meta.data$cell_types)

levels(ref_scrna_final@meta.data$cell_types) = c("Adult astrocytes", "Fetal astrocytes", "Fetal excitatory neurons", "Adult excitatory neurons", "Fetal cardiac smooth muscle cells", "Adult interneurons", "Fetal interneurons", "Fetal other smooth muscle cells", "Fetal neuronal stem cells", "Adult oligodendrocytes", "Adult oligodendrocyte progenitor cells", "Fetal oligodendrocyte progenitor cells", "PGC-like", "Adult other smooth muscle cells", "Spermatid", "Spermatocytes", "Spermatogonia")

saveRDS(ref_scrna_final, 'reference_scrna.rds')

#ref_scrna_final = readRDS('reference_scrna.rds')

ref_scrna_final@meta.data$cell_types = as.factor(ref_scrna_final@meta.data$cell_types)
ref_scrna_final@meta.data$cell_types = factor(ref_scrna_final@meta.data$cell_types, levels = rev(c("PGC-like", "Spermatogonia", "Spermatocytes", "Spermatid", "Fetal neuronal stem cells", "Fetal oligodendrocyte progenitor cells", "Adult oligodendrocyte progenitor cells", "Adult oligodendrocytes", "Fetal astrocytes", "Adult astrocytes", "Fetal interneurons", "Adult interneurons", "Fetal excitatory neurons", "Adult excitatory neurons", "Fetal cardiac smooth muscle cells", "Fetal other smooth muscle cells", "Adult other smooth muscle cells")))

##limit gct tpm counts to genes shared between datasets
gct_tpm_filtered = gct_tpm[row.names(ref_scrna_final),]
gct_log2tpm <- log2(as.matrix(gct_tpm_filtered) + 1)
rna_cuts = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/supplementary_data/GCT_RNA_metadata.txt', header = T, sep = '\t', stringsAsFactors = F)
row.names(rna_cuts) = rna_cuts$Sample_ID

#2. Convert reference data to tpm values
ref_scrna_final_df = as.data.frame(ref_scrna_final@assays$RNA@counts)
ref_scrna_final_tpm <- t( t(ref_scrna_final_df) * 1e6 / colSums(ref_scrna_final_df)) #tpm conversion, 
ref_scrna_final_log2tpm <- log2(as.matrix(ref_scrna_final_tpm) + 1)
write.table(ref_scrna_final_log2tpm, 'reference_scrna_log2tpm.txt', col.names = T, row.names = T, sep = '\t', quote = F)

ref_scrna_final@meta.data$sample = row.names(ref_scrna_final@meta.data)
ref_scrna_final@meta.data$source = 'scrna'
ref_scrna_final_meta = ref_scrna_final@meta.data[, c('sample', 'cell_types', 'source')]

rna_cuts$cell_types = rna_cuts$Updated.Description
rna_cuts$sample = row.names(rna_cuts)
rna_cuts$source = 'lcm'
gct_final_meta = rna_cuts[, c('sample', 'cell_types', 'source')]

#3. Combined tpm values from our data and the reference data
all_data_final_log2tpm = cbind(ref_scrna_final_log2tpm, gct_log2tpm)
all_data_final_meta = rbind(ref_scrna_final_meta, gct_final_meta)

colnames(all_data_final_log2tpm) = gsub('\\.', '_', colnames(all_data_final_log2tpm))
row.names(all_data_final_meta) = gsub('-', '_', row.names(all_data_final_meta))
row.names(all_data_final_meta) = gsub('\\.', '_', row.names(all_data_final_meta))

#germ cell comparison
gc_exp = all_data_final_log2tpm[,row.names(all_data_final_meta[all_data_final_meta$cell_types %in% c('PGC-like', 'Spermatogonia', 'Spermatocytes', 'Spermatid', 'GCNIS', 'Seminoma', 'Dysgerminoma', 'Seminiferous_tubule'), ])]
gc_exp.t = as.data.frame(t(as.data.frame(gc_exp)))
gc_exp.t$cell_types = all_data_final_meta[row.names(gc_exp.t),]$cell_types
gc_exp.t$sample = row.names(gc_exp.t)
gc_exp.t.m = reshape2::melt(gc_exp.t[, c('cell_types', 'sample', 'POU5F1', 'NANOG', 'SOX17', 'DDX4', 'TNP1')], id.vars = c('cell_types', 'sample'))

gc_exp.t.m$gene_use = NA
for(i in 1:nrow(gc_exp.t.m)){
  if(gc_exp.t.m$variable[i] %in% c('POU5F1', 'NANOG', 'SOX17')) gc_exp.t.m$gene_use[i] = 'fetal'
  if(gc_exp.t.m$variable[i] == 'DDX4') gc_exp.t.m$gene_use[i] = 'gen_gc'
  if(gc_exp.t.m$variable[i] == 'TNP1') gc_exp.t.m$gene_use[i] = 'spermatid'
}

gc_exp.t.m$gene_use = factor(gc_exp.t.m$gene_use, levels = c('fetal', 'gen_gc', 'spermatid'))
gc_exp.t.m$cell_types = factor(gc_exp.t.m$cell_types, levels = c('PGC-like', 'Spermatogonia', 'Spermatocytes', 'Spermatid', 'GCNIS', 'Seminoma', 'Dysgerminoma', 'Seminiferous_tubule'))

aggregate(value ~ cell_types + variable, gc_exp.t.m, length) #Dysgerminoma will need individual data points for DDX4 and TNP1 plotting

p1 = ggplot(gc_exp.t.m, mapping = aes(x = gene_use, y = value)) + 
  geom_quasirandom(dodge.width=.8, cex=1, aes(x = gene_use, y = value), shape = 21, fill = 'transparent', col = 'transparent') + 
  geom_beeswarm(data = gc_exp.t.m[gc_exp.t.m$cell_types == 'Dysgerminoma' & gc_exp.t.m$variable %in% c('DDX4', 'TNP1'), ], mapping = aes(x = gene_use, y = value), dodge.width=1, cex=1, shape = 21, fill = 'black', col = 'black') + 
  facet_wrap(. ~ cell_types, strip.position = 'bottom', ncol = 8) +
  theme_pubr() + 
  theme(axis.text.x = element_blank(), 
        panel.spacing = unit(0, "mm"), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        legend.position = 'none', 
        axis.ticks.x = element_blank()) +
  labs(y = 'Log2(TPM) expression', x = '') + 
  geom_pointrange(mapping = aes(col = gene_use), stat = "summary", shape=19, fun.min = function(z) {quantile(z,0.25)}, 
                  fun.max = function(z) {quantile(z,0.75)}, 
                  fun = median, size = 0.8, fatten = 5, linetype = 2) + 
  scale_color_manual(values = c('#FF6C5C', '#58CCED', '#072F5F'))

#neural comparison
ne_exp = all_data_final_log2tpm[,row.names(all_data_final_meta[all_data_final_meta$cell_types %in% c('Neuroepithelium', 'Fetal interneurons', 'Fetal neuronal stem cells', 'Fetal oligodendrocyte progenitor cells', 'Fetal astrocytes', 'Fetal excitatory neurons', 'Adult interneurons', 'Adult oligodendrocyte progenitor cells', 'Adult astrocytes', 'Adult excitatory neurons', 'Adult oligodendrocytes'), ])]
ne_exp.t = as.data.frame(t(as.data.frame(ne_exp)))
ne_exp.t$cell_types = all_data_final_meta[row.names(ne_exp.t),]$cell_types
ne_exp.t$sample = row.names(ne_exp.t)
ne_exp.t.m = reshape2::melt(ne_exp.t[, c('cell_types', 'sample', 'PAX6', 'SOX2', 'NES', 'MBP', 'SLC1A2', 'RBFOX3', 'GAD1')], id.vars = c('cell_types', 'sample'))

ne_exp.t.m$gene_use = NA
for(i in 1:nrow(ne_exp.t.m)){
  if(ne_exp.t.m$variable[i] %in% c('PAX6', 'SOX2', 'NES')) ne_exp.t.m$gene_use[i] = 'fetal'
  if(ne_exp.t.m$variable[i] == 'MBP') ne_exp.t.m$gene_use[i] = 'oligo'
  if(ne_exp.t.m$variable[i] == 'SLC1A2') ne_exp.t.m$gene_use[i] = 'astro'
  if(ne_exp.t.m$variable[i] == 'RBFOX3') ne_exp.t.m$gene_use[i] = 'neun'
  if(ne_exp.t.m$variable[i] == 'GAD1') ne_exp.t.m$gene_use[i] = 'inter'
}

ne_exp.t.m$gene_use = factor(ne_exp.t.m$gene_use, levels = c('fetal', 'oligo', 'astro', 'neun', 'inter'))
ne_exp.t.m$cell_types = factor(ne_exp.t.m$cell_types, levels = c('Fetal neuronal stem cells', 'Fetal oligodendrocyte progenitor cells', 'Fetal astrocytes', 'Fetal excitatory neurons', 'Fetal interneurons', 'Adult oligodendrocyte progenitor cells', 'Adult oligodendrocytes', 'Adult astrocytes', 'Adult excitatory neurons', 'Adult interneurons', 'Neuroepithelium'))

aggregate(value ~ cell_types + variable, ne_exp.t.m, length)

p2 = ggplot(ne_exp.t.m, mapping = aes(x = gene_use, y = value)) + 
  geom_quasirandom(dodge.width=.8, cex=1, aes(x = gene_use, y = value), shape = 21, fill = 'transparent', col = 'transparent') + 
  facet_wrap(. ~ cell_types, strip.position = 'bottom', ncol = 11) +
  theme_pubr() + 
  theme(axis.text.x = element_blank(), 
        panel.spacing = unit(0, "mm"), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        legend.position = 'none', 
        axis.ticks.x = element_blank()) +
  labs(y = 'Log2(TPM) expression', x = '') + 
  geom_pointrange(mapping = aes(col = gene_use), stat = "summary", shape=19, fun.min = function(z) {quantile(z,0.25)}, 
                  fun.max = function(z) {quantile(z,0.75)}, 
                  fun = median, size = 0.8, fatten = 5, linetype = 2) + 
  scale_color_manual(values = c('#FF6C5C', '#58CCED', '#3895D3', '#1261A0', '#072F5F'))

#smooth muscle comparison
sm_exp = all_data_final_log2tpm[,row.names(all_data_final_meta[all_data_final_meta$cell_types %in% c("Fetal cardiac smooth muscle cells", "Fetal other smooth muscle cells", "Adult other smooth muscle cells", "Smooth_muscle_teratoma"), ])]
sm_exp.t = as.data.frame(t(as.data.frame(sm_exp)))
sm_exp.t$cell_types = all_data_final_meta[row.names(sm_exp.t),]$cell_types
sm_exp.t$sample = row.names(sm_exp.t)
sm_exp.t.m = reshape2::melt(sm_exp.t[, c('cell_types', 'sample', 'IGF2', 'TAGLN', 'ACTA2', 'MYH11')], id.vars = c('cell_types', 'sample'))

sm_exp.t.m$gene_use = NA
for(i in 1:nrow(sm_exp.t.m)){
  if(sm_exp.t.m$variable[i] %in% c('TAGLN', 'ACTA2', 'MYH11')) sm_exp.t.m$gene_use[i] = 'adult'
  if(sm_exp.t.m$variable[i] == 'IGF2') sm_exp.t.m$gene_use[i] = 'fetal'
}

sm_exp.t.m$gene_use = factor(sm_exp.t.m$gene_use, levels = c('fetal', 'adult'))
sm_exp.t.m$cell_types = factor(sm_exp.t.m$cell_types, levels = c("Fetal cardiac smooth muscle cells", "Fetal other smooth muscle cells", "Adult other smooth muscle cells", "Smooth_muscle_teratoma"))

aggregate(value ~ cell_types + variable, sm_exp.t.m, length) #add data points for GCT data for IGF2

p3 = ggplot(sm_exp.t.m, mapping = aes(x = gene_use, y = value)) + 
  geom_quasirandom(dodge.width=.8, cex=1, aes(x = gene_use, y = value), shape = 21, fill = 'transparent', col = 'transparent') + 
  geom_beeswarm(data = sm_exp.t.m[sm_exp.t.m$cell_types == 'Smooth_muscle_teratoma' & sm_exp.t.m$variable == 'IGF2', ], mapping = aes(x = gene_use, y = value), dodge.width=1, cex=1, shape = 21, fill = 'black', col = 'black') + 
  facet_wrap(. ~ cell_types, strip.position = 'bottom', ncol = 8) +
  theme_pubr() + 
  theme(axis.text.x = element_blank(), 
        panel.spacing = unit(0, "mm"), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        legend.position = 'none', 
        axis.ticks.x = element_blank()) +
  labs(y = 'Log2(TPM) expression', x = '') + 
  geom_pointrange(mapping = aes(col = gene_use), stat = "summary", shape=19, fun.min = function(z) {quantile(z,0.25)}, 
                  fun.max = function(z) {quantile(z,0.75)}, 
                  fun = median, size = 0.8, fatten = 5, linetype = 2) + 
  scale_color_manual(values = c('#FF6C5C', '#58CCED'))

#4. Plot

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/log2tpm_expression_gct_vs_ref_tissue20210305.pdf', width = 14, height = 9, useDingbats = F)
ggdraw() +
  draw_plot(p1, 0, 2/3, 8/11, 1/3) +
  draw_plot(p2, 0, 1/3, 1, 1/3) +
  draw_plot(p3, 0, 0, 4.25/11, 1/3) +
  draw_plot_label(c("A", "B", "C"), c(0, 0, 0), c(1, 2/3, 1/3), size = 15)
dev.off()