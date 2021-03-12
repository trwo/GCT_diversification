#RNA data filtering & pre_processing
##Filtering
###Raw count table, preQC with genes mapped to nuclear genome only.

library(ggplot2)
library(ggpubr)

raw_data = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/filtering/germ_cell_lcm_rna_preQC.txt', header = T, sep = '\t', stringsAsFactors = F) #just normal testis and GCT raw count data, 500 microbiopsies

#remove transcriptomes that were:
#- identified as contaminated during microdissection
#- non-neoplastic tissues that are not healthy seminiferous tubules
#- uncertainty as to the histological categorisation of the microbiopsy

raw_data = raw_data[, !colnames(raw_data) %in% c("PD43298a_Z1_2_SL02_04_SEC02_04_F2", 'PR43299c_lo0009', 'PR43299c_lo0010', 'PR43299c_lo0011', 'PR43299c_lo0012', "PD43296a_SL02_SEC02_G1", "PD43296a_Z1_2_SL02_04_SEC02_04_E1A", "PD43298a_SL02_SEC02_B3", "PD43298a_SL02_SEC02_C3", "PD43298a_Z1_2_SL02_04_SEC02_04_G2", "PR43299c_lo0002", "PR43299c_lo0006", "PR42036a_SL0D_SEC02_A5", "PR42036a_SL0D_SEC02_G3", "PR42036a_SL0D_SEC02_H3", "PR43298a_SL04_SEC04_D5",             "PR43298a_SL04_SEC04_E5" ,"PR43298a_SL04_SEC04_G9", "PR45543a_SL02_SEC03_B4", "PR45543a_SL02_SEC03_D2", "PR45543a_SL02_SEC03_F2", "PR46269c_SL02_SEC03_E12", "PR46269c_SL02_SEC03_F12", "PR46270c_SL02_SEC03_B5", "PR46270c_SL02_SEC03_C5", "PR46270c_SL02_SEC03_D5", "PR46270c_SL02_SEC03_E3", "PR46270c_SL02_SEC03_F3",             "PR46270c_SL02_SEC03_G3", "PR46270c_SL02_SEC03_H3", "PR46967a_SL02_SEC03_D8", "PR46967a_SL02_SEC03_E4", "PR46967a_SL02_SEC03_G6", "PR46968d_SL02_SEC02_03_H5", "PR46969c_SL02_SEC03_C11", "PR46969c_SL02_SEC03_D11", "PR46969c_SL02_SEC03_E11", "PR46969c_SL02_SEC03_F11", "PR46969c_SL02_SEC03_G11", "PR46969c_SL02_SEC03_H11")] #460 samples left

#how many reads, across all 55502 mapped features do we have?
summary(colSums(raw_data[, 7:ncol(raw_data)])) #highly variable library sizes
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   682   86612  246980  502107  651605 3863511 

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/all_samples_read_count_QC.pdf', height = 6, width = 6, useDingbats = F)
ggplot() +
  geom_histogram(mapping = aes(colSums(raw_data[, 7:ncol(raw_data)])), bins = 100, fill = '#999999') + theme_pubr() + coord_cartesian(expand = F) + xlab('Number of reads across all mapped features') + ylab('Number of microbiopsies')
#spike in samples with under 1000 features expressed
dev.off()

#if we choose a threshold of 5 reads as the minimum to consider a gene truly expressed, how many are expressed per microbiopsy?
summary(colSums(raw_data[, 7:ncol(raw_data)] >= 5))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0    4150    8916    9077   13654   22865 

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/all_samples_feature_expression_QC.pdf', height = 6, width = 6, useDingbats = F)
ggplot() +
  geom_histogram(mapping = aes(colSums(raw_data[, 7:ncol(raw_data)] >= 5)), bins = 100, fill = '#999999') + theme_pubr() + coord_cartesian(expand = F) + geom_vline(xintercept = 1000, lty = 2, col = 'red') + xlab('Features expressed per microbiopsy to a depth of at least 5 reads') + ylim(0, 15) + ylab('Number of microbiopsies')
#spike in samples with under 1000 features expressed
dev.off()

###Pre-processing 

library(limma)
library(edgeR)
library(Seurat)

raw_data = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/filtering/germ_cell_lcm_rna_preQC.txt', header = T, sep = '\t', stringsAsFactors = F) #raw count data
metadata = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/supplementary_data/GCT_RNA_metadata.txt', header = T, sep = '\t', stringsAsFactors = F) #read in metadata for filtered, final cuts
row.names(metadata) = metadata$Sample_ID
gct_mat = raw_data[, 7:ncol(raw_data)]
row.names(gct_mat) = raw_data$Geneid

#remove cuts that failed QC
gct_mat = gct_mat[, colnames(gct_mat) %in% metadata$Sample_ID] #416 samples

#define factors that will be input into the model
histology = factor(metadata[colnames(gct_mat), ]$Updated.Description)
patient = factor(metadata[colnames(gct_mat), ]$Case_ID)
metadata = metadata[colnames(gct_mat),]

#pre-processing
y <- DGEList(counts = gct_mat, group = histology)
keep <- filterByExpr.default(y, group = histology)
keep <- filterByExpr(y, group = histology)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y) #33185 genes kept
y <- estimateCommonDisp(y) #to generate pseudocounts, adjusted for library size
pseudo_counts = y$pseudo.counts

#check that the tmm-normalised counts cluster by histology before performing DE analysis
gct_obj <- CreateSeuratObject(counts = y$pseudo.counts, meta.data = metadata, project = 'gct')
gct_obj <- NormalizeData(gct_obj, normalization.method = "LogNormalize", scale.factor = 10000)
gct_obj <- FindVariableFeatures(gct_obj, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(gct_obj )
gct_obj  <- ScaleData(gct_obj , features = all.genes)
gct_obj <- RunPCA(gct_obj, features = VariableFeatures(object = gct_obj), npcs = 100)
DimPlot(gct_obj, reduction = "pca", group.by = 'Updated.Description')
ElbowPlot(gct_obj, ndims = 100)
gct_obj <- RunUMAP(gct_obj, dims = 1:30)

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/gct_seurat_30_dims_histology_umap_20210304.pdf', width = 8, height = 6, useDingbats = F)
DimPlot(gct_obj, reduction = "umap", group.by = 'Updated.Description')
dev.off()

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/gct_seurat_30_dims_patient_umap_20210304.pdf', width = 8, height = 6, useDingbats = F)
DimPlot(gct_obj, reduction = "umap", group.by = 'Case_ID')
dev.off()

saveRDS(gct_obj, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/pre_processing/lcm_gct_seurat_object_20210403.rds')

#clusters by histology but also clearly some patient-specific effect which will need to be accounted for in the differential expression analysis

#save filtered raw counts and tmm-normalised counts, plus calculate tpm values
write.table(pseudo_counts, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/pre_processing/lcm_rna_tmm_norm_counts_20210304.txt', quote = F, col.names = T, row.names = T, sep = '\t')

write.table(y$counts, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/pre_processing/lcm_rna_raw_counts_20210304.txt', quote = F, col.names = T, row.names = T, sep = '\t')

row.names(raw_data) = raw_data$Geneid
edgeR_gene_data = raw_data[row.names(pseudo_counts),]
x = pseudo_counts
x = x / edgeR_gene_data$Length #normalise for transcript length
tpm <- t( t(x) * 1e6 / colSums(x)) #normalise to read depth
write.table(tpm, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/RNA/pre_processing/tmm_norm_gct_tpm_20210304.txt', col.names = T, row.names = T, sep = '\t', quote = F)

metadata$total_reads_filtered_genes = NA

for(i in 1:nrow(metadata)){
  metadata$total_reads_filtered_genes[i] = sum(y$counts[, row.names(metadata)[i]])
}

write.table(metadata, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/supplementary_data/GCT_RNA_metadata.txt', quote = F, col.names = T, row.names = F, sep = '\t')