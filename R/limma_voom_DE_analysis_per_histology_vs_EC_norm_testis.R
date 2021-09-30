setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/limma_voom/')

library(limma)
library(edgeR)

raw_data = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/original_nature_submission/analyses/RNA/filtering/germ_cell_lcm_rna_preQC.txt', header = T, sep = '\t', stringsAsFactors = F) #just normal testis and GCT raw count data, 500 microbiopsies
metadata = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/supplementary_data/mRNA_summary_metadata.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '') #read in metadata for filtered, final cuts
row.names(metadata) = metadata$Sample
gct_mat = raw_data[, 7:ncol(raw_data)]
row.names(gct_mat) = raw_data$Geneid

#remove cuts that failed QC
gct_mat = gct_mat[, colnames(gct_mat) %in% metadata$Sample] #416 samples

#define factors that will be input into the model
histology = factor(metadata[colnames(gct_mat), ]$Updated.Description)
patient = factor(metadata[colnames(gct_mat), ]$Case.ID)
metadata = metadata[colnames(gct_mat),]

#preprocessing
y <- DGEList(counts = gct_mat, group = histology)
#keep <- filterByExpr.default(y, group = histology)
keep <- filterByExpr(y, group = histology)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y) #33185 genes kept

plotMDS(y, col=as.numeric(patient), labels = histology) #no obvious outlier samples

#library sizes are variable for voom is preferred over limma-trend
#specify model to be fitted, no intercept to be fitted, each histology will have its own group mean
design = model.matrix(~ 0 + histology)
colnames(design) = levels(histology)

contr.matrix <- makeContrasts(
  fat_v_ec = Adipose_teratoma - Embryonal_carcinoma,
  dys_v_ec = Dysgerminoma - Embryonal_carcinoma,
  gcnis_v_ec = GCNIS - Embryonal_carcinoma,
  cartilage_v_ec = Cartilage_teratoma - Embryonal_carcinoma, 
  glandular_v_ec = Mature_glandular_teratoma - Embryonal_carcinoma, 
  ne_v_ec = Neuroepithelium - Embryonal_carcinoma, 
  other_ter_v_ec = Other_epithelial_teratoma - Embryonal_carcinoma, 
  muscle_v_ec = Smooth_muscle_teratoma - Embryonal_carcinoma, 
  sct_v_ec = Syncytiotrophoblasts - Embryonal_carcinoma, 
  yst_v_ec = Yolk_sac_tumour - Embryonal_carcinoma, 
  sem_v_ec = Seminoma - Embryonal_carcinoma,
  str_v_ec = Stroma - Embryonal_carcinoma,
  fat_v_st = Adipose_teratoma - Seminiferous_tubule,
  dys_v_st = Dysgerminoma - Seminiferous_tubule,
  ec_v_st = Embryonal_carcinoma - Seminiferous_tubule,
  gcnis_v_st = GCNIS - Seminiferous_tubule,
  cartilage_v_st = Cartilage_teratoma - Seminiferous_tubule, 
  glandular_v_st = Mature_glandular_teratoma - Seminiferous_tubule, 
  ne_v_st = Neuroepithelium - Seminiferous_tubule, 
  other_ter_v_st = Other_epithelial_teratoma - Seminiferous_tubule, 
  muscle_v_st = Smooth_muscle_teratoma - Seminiferous_tubule, 
  sct_v_st = Syncytiotrophoblasts - Seminiferous_tubule, 
  yst_v_st = Yolk_sac_tumour - Seminiferous_tubule, 
  sem_v_st = Seminoma - Seminiferous_tubule,
  str_v_st = Stroma - Seminiferous_tubule,
  levels = colnames(design))

v <- voom(y, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

saveRDS(v, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/limma_voom/limma_voom_de_LCM_GCT_DE_voom_20210925.rds')
saveRDS(efit, '/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/limma_voom/limma_voom_de_LCM_GCT_DE_model_fit_20210925.rds')

for(i in colnames(contr.matrix)){
  data = topTable(efit, coef = i, n = Inf, sort = "p")
  write.table(data, paste0(i, "_voom_de_20210925.txt"), col.names = T, row.names = T, sep = '\t', quote = F)
}