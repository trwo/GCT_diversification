#Cytoband-level gene expression by GCT tissue

library(msigdbr)
library(limma)
library(scales)
library(ggpubr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/limma_voom_de')

gct_v = readRDS('limma_voom_de_LCM_GCT_DE_voom_object_20210304.RDS') #original limma-voom object for DE
gct_fit = readRDS('limma_voom_de_LCM_GCT_DE_model_fit_20210304.RDS') #original limma-voom object for DE

contr.matrix <- makeContrasts(
  fat_v_st = Adipose_teratoma - Seminiferous_tubule,
  dys_v_st = Dysgerminoma - Seminiferous_tubule,
  ec_v_st = Embryonal_carcinoma - Seminiferous_tubule,
  gcnis_v_st = GCNIS - Seminiferous_tubule,
  cartilage_v_st = Immature_cartilage_teratoma - Seminiferous_tubule, 
  glandular_v_st = Mature_glandular_teratoma - Seminiferous_tubule, 
  ne_v_st = Neuroepithelium - Seminiferous_tubule, 
  other_ter_v_st = Other_epithelial_teratoma - Seminiferous_tubule, 
  muscle_v_st = Smooth_muscle_teratoma - Seminiferous_tubule, 
  sct_v_st = Syncytiotrophoblasts - Seminiferous_tubule, 
  yst_v_st = Yolk_sac_tumour - Seminiferous_tubule, 
  sem_v_st = Seminoma - Seminiferous_tubule,
  str_v_st = Stroma - Seminiferous_tubule,
  levels = colnames(gct_fit$design))

all_hum_gene_sets = msigdbr(species = "Homo sapiens") #get genes per cytoband

c1_gene_set = list()

for(i in unique(all_hum_gene_sets[all_hum_gene_sets$gs_cat == 'C1', ]$gs_name)){
  c1_gene_set[[i]] = all_hum_gene_sets[all_hum_gene_sets$gs_cat == 'C1' & all_hum_gene_sets$gs_name == i, ]$human_gene_symbol
}
c1_idx <- ids2indices(c1_gene_set,id=rownames(gct_v))

all_histo_regions_enrich = c()
for(i in 1:ncol(contr.matrix)){
  temp = camera(gct_v, c1_idx, gct_fit$design, contrast = contr.matrix[,i])
  tissue = names(contr.matrix[,i][contr.matrix[,i] == 1])
  temp$tissue = tissue
  temp$cytoband = row.names(temp)
  row.names(temp) = NULL
  all_histo_regions_enrich = rbind(all_histo_regions_enrich, temp)
}

all_histo_regions_enrich$chr = substr(unlist(lapply(strsplit(all_histo_regions_enrich$cytoband, 'p|q'), '[[', 1)), 4, 5)
all_histo_regions_enrich$chr = factor(all_histo_regions_enrich$chr, levels = c(1:22, 'X', 'Y'))
all_histo_regions_enrich$arm = substr(unlist(lapply(strsplit(all_histo_regions_enrich$cytoband, paste(unique(paste0('chr', all_histo_regions_enrich$chr)), collapse = '|')), '[[', 2)), 0, 1)
all_histo_regions_enrich$cytoband_pos = substr(unlist(lapply(strsplit(all_histo_regions_enrich$cytoband, 'p|q'), '[[', 2)), 0, 2)

all_histo_regions_enrich = all_histo_regions_enrich[order(all_histo_regions_enrich$tissue, all_histo_regions_enrich$chr, all_histo_regions_enrich$arm, all_histo_regions_enrich$cytoband_pos),]
all_histo_regions_enrich$cytoband = factor(all_histo_regions_enrich$cytoband, levels = unique(all_histo_regions_enrich$cytoband))
all_histo_regions_enrich$sig_variable = 1 / all_histo_regions_enrich$FDR #invert for plotting purposes

for(i in 1:nrow(all_histo_regions_enrich)){
  if(all_histo_regions_enrich$Direction[i] == 'Down') all_histo_regions_enrich$sig_variable[i] = -1 * all_histo_regions_enrich$sig_variable[i]
}

write.table(all_histo_regions_enrich, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/supplementary_data/gct_cytoband_region_enrichment20210305.txt', col.names = T, row.names = F, quote = F, sep = '\t')

all_histo_regions_enrich_wide <- dcast(all_histo_regions_enrich, tissue ~ cytoband, value.var="sig_variable")
row.names(all_histo_regions_enrich_wide) = all_histo_regions_enrich_wide$tissue

all_histo_regions_enrich_wide = all_histo_regions_enrich_wide[,2:297]
all_histo_regions_enrich_wide_auto = all_histo_regions_enrich_wide[, !grepl('chrX|chrY', colnames(all_histo_regions_enrich_wide)),]

col_fun = colorRamp2(c(-100, -10, 0, 10, 100), c('purple', '#520160', 'black', '#435E1C', '#95D840'))
chrom_order = factor(unlist(lapply(strsplit(colnames(all_histo_regions_enrich_wide_auto), 'p|q'), '[[', 1)), levels = paste0('chr', 1:22))
levels(chrom_order) = 1:22

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/enriched_regions_gct_expression20210305.pdf', width = 15, height = 8)
Heatmap(all_histo_regions_enrich_wide_auto, cluster_columns = F, cluster_rows = F, col = col_fun, column_split = chrom_order, heatmap_legend_param = list(title = "FDR", at = c(-100, -10, 10, 100), labels = c(0.01, 0.1, 0.1, 0.01), legend_height = unit(4, "cm")), show_column_names = F, row_names_side = 'left', row_order = c('GCNIS', 'Dysgerminoma', 'Seminoma', 'Embryonal_carcinoma', 'Yolk_sac_tumour', 'Syncytiotrophoblasts', 'Neuroepithelium', 'Mature_glandular_teratoma', 'Other_epithelial_teratoma', 'Adipose_teratoma', 'Smooth_muscle_teratoma', 'Immature_cartilage_teratoma', 'Stroma'))
dev.off()
