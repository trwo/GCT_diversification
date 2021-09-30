library(limma)
library(edgeR)

raw_data = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/original_nature_submission/analyses/RNA/filtering/germ_cell_lcm_rna_preQC.txt', header = T, sep = '\t', stringsAsFactors = F) #just normal testis and GCT raw count data, 500 microbiopsies
metadata = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/supplementary_data/mRNA_summary_metadata.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '') #read in metadata for filtered, final cuts
row.names(metadata) = metadata$Sample
gct_mat = raw_data[, 7:ncol(raw_data)]
row.names(gct_mat) = raw_data$Geneid

#remove cuts that failed QC
gct_mat = gct_mat[, colnames(gct_mat) %in% metadata[metadata$Updated.Description != 'GCNIS',]$Sample] #380 samples after removing GCNIS

#define factors that will be input into the model
histology = factor(metadata[colnames(gct_mat), ]$Updated.Description)
tumour = factor(metadata[colnames(gct_mat), ]$Normal.Tumour)
metadata = metadata[colnames(gct_mat),]

#preprocessing
y <- DGEList(counts = gct_mat, group = tumour)
keep <- filterByExpr(y, group = tumour)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y) #12394 genes kept
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
corfit$consensus #0.1987738

setwd('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/nature_resub/limma_voom')

saveRDS(v, 'LCM_invasive_GCT_v_testis_DE_voom_object_20210812.rds')
saveRDS(fit, 'LCM_invasive_GCT_v_testis_DE_model_fit_20210812.rds')

v = readRDS('LCM_invasive_GCT_v_testis_DE_voom_object_20210812.rds')
fit = readRDS('LCM_invasive_GCT_v_testis_DE_model_fit_20210812.rds')

cont_matrix = makeContrasts(Tumour - Normal, levels = colnames(coef(fit)))
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)
top.table <- topTable(fit2, number = nrow(fit$coefficients), sort.by = "P")

raw_data = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/original_nature_submission/analyses/RNA/filtering/germ_cell_lcm_rna_preQC.txt', header = T, sep = '\t', stringsAsFactors = F)
row.names(raw_data) = raw_data$Geneid
raw_data$Start = unlist(lapply(strsplit(raw_data$Start, ';'), '[[', 1))

top.table$chr = raw_data[row.names(top.table),]$Chr
top.table$start = raw_data[row.names(top.table),]$Start
top.table$start = as.numeric(top.table$start)

write.table(top.table, 'invasive_gct_v_normtestis_voom_de_20210812.txt', col.names = T, sep = '\t', row.names = T, quote = F)