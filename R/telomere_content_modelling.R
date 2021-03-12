#Telomere content modelling
##There are many biological and technical factors that influence telomere length estimates. We first review the raw content values and their association with metadata parameters before we construct a linear mixed effects model to adjust for factors that are not of interest.

library(lme4)
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(remef)

telomerehunter_df = read.table('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/analyses/telomerehunter/gct_telomerehunter_matrix_20201110.txt', sep = '\t', header = T, stringsAsFactors = F) #read in telomerehunter content values for cohort
telomerehunter_df = telomerehunter_df[telomerehunter_df$patientID != 'PD46967',] #excluded as DNA samples were not included in the study from this patient due to the lack of a matched normal

#make dysgerminoma/seminoma into one category so model as site is already another variable that will discriminate these tissues between testis and ovary
telomerehunter_df[telomerehunter_df$histology %in% c('Dysgerminoma', 'Seminoma'),]$histology = 'Dysgerminoma/Seminoma'

#plot the variables to assess whether they may impact telomere length
##platform effect
ggplot(telomerehunter_df[telomerehunter_df$normal.tumour == 'tumour',]) + geom_point(mapping = aes(x = histology, y = tel_content, col = platform)) #hiseq seems to have slightly higher telomere lengths
ggplot(telomerehunter_df[telomerehunter_df$normal.tumour == 'tumour',]) + geom_point(mapping = aes(x = histology, y = tel_content, col = platform, shape = site))
ggplot(telomerehunter_df[telomerehunter_df$histology == 'normal_ST',]) + geom_point(mapping = aes(x = patientID, y = tel_content, col = platform))
ggplot(telomerehunter_df[telomerehunter_df$histology == 'normal_ST',]) + geom_point(mapping = aes(x = age, y = tel_content, col = platform))
ggplot(telomerehunter_df[telomerehunter_df$normal.tumour == 'tumour' & telomerehunter_df$platform == 'NovaSeq',]) + geom_point(mapping = aes(x = histology, y = tel_content, col = site)) #even in Novaseq only, the telomeres of ovarian tumours were shorter

##age
ggplot(telomerehunter_df[telomerehunter_df$normal.tumour == 'tumour',]) + geom_point(mapping = aes(x = age, y = tel_content, col = histology)) #age effect not appreciable
ggplot(telomerehunter_df[telomerehunter_df$normal.tumour == 'normal',]) + geom_point(mapping = aes(x = age, y = tel_content, col = histology))

##coverage
ggplot(telomerehunter_df[telomerehunter_df$normal.tumour == 'tumour',]) + geom_point(mapping = aes(x = coverage, y = tel_content, col = histology)) #coverage effect not overtly appreciable
ggplot(telomerehunter_df[telomerehunter_df$normal.tumour == 'normal',]) + geom_point(mapping = aes(x = coverage, y = tel_content, col = histology))

##diagnosis
ggplot(telomerehunter_df[telomerehunter_df$normal.tumour == 'tumour',]) + geom_point(mapping = aes(x = histology, y = tel_content, col = diagnosis, shape = site)) #differences perhaps confounded by sequencing platform

#patient
ggplot(telomerehunter_df[telomerehunter_df$normal.tumour == 'tumour',]) + geom_point(mapping = aes(x = histology, y = tel_content, col = patientID, shape = platform)) #reasonable consistency in telomere lengths across histologies between patients

##site
ggplot(telomerehunter_df[telomerehunter_df$normal.tumour == 'tumour',]) + geom_point(mapping = aes(x = site, y = tel_content, col = patientID, shape = platform)) #all ovarian telomeres are considerably shorter

#create basic models to compare
##null model
null_model = lmer(formula = tel_content ~ 0 + (1 | patientID), data = telomerehunter_df, REML = F)

##full model
full_model = lmer(formula = tel_content ~ 0 + (1 | patientID) + coverage + age + platform + diagnosis + histology + site, data = telomerehunter_df, REML = F)

##compare null and full
anova(null_model, full_model) #full better than null, fixed effects do seem to impact telomere length, as expected

##optimise model using stepwise algorithm
step(full_model) #optimum model = tel_content ~ (1 | patientID) + platform + diagnosis + histology + site - 1, age and coverage not very informative
final_model = lmer(formula = tel_content ~ (1 | patientID) + platform + diagnosis + histology + site - 1, data = telomerehunter_df, REML = F)
anova(null_model, full_model, final_model) 
summary(final_model)

#partial regression - to remove the effect of patient and sequencing platform
full_model.no_patientID_platform_partial = remef(final_model, fix = c("platformHiSeq", "platformNovaSeq"), ran = 'all')
telomerehunter_df$tel_content_partial = full_model.no_patientID_platform_partial
telomerehunter_df$diagnosis_plot = telomerehunter_df$diagnosis
telomerehunter_df[telomerehunter_df$diagnosis %in% c('blood', 'normal_ST'),]$diagnosis_plot = 'Normal'

telomerehunter_df$histology = factor(telomerehunter_df$histology, levels = c('blood', 'normal_ST', 'GCNIS', 'Dysgerminoma/Seminoma', 'Embryonal_carcinoma', 'Yolk_sac_tumour', 'Neuroectoderm', 'Cartilage_teratoma', 'Glandular_teratoma', 'Malignant_stroma', 'Smooth_muscle_teratoma', 'Syncytiotrophoblasts'))
levels(telomerehunter_df$histology) = c('Blood', 'Seminiferous tubule', 'GCNIS', 'Dysgerminoma/Seminoma', 'Embryonal carcinoma', 'Yolk sac tumour', 'Neuroepithelium', 'Cartilage teratoma', 'Glandular teratoma', 'Malignant stroma', 'Smooth muscle teratoma', 'Syncytiotrophoblasts')

telomerehunter_df$diagnosis_plot_reannot = telomerehunter_df$diagnosis_plot
telomerehunter_df[telomerehunter_df$diagnosis_plot == 'Seminoma' & telomerehunter_df$histology == 'Dysgerminoma/Seminoma',]$diagnosis_plot_reannot = 'Seminoma, pure seminoma'
telomerehunter_df[telomerehunter_df$diagnosis_plot == 'Post-pubertal mixed germ cell tumour' & telomerehunter_df$histology == 'Dysgerminoma/Seminoma',]$diagnosis_plot_reannot = 'Dysgerminoma/Seminoma, NSGCT'
telomerehunter_df[telomerehunter_df$diagnosis_plot == 'Post-pubertal mixed germ cell tumour' & !(telomerehunter_df$histology %in% c('GCNIS', 'Dysgerminoma/Seminoma')),]$diagnosis_plot_reannot = 'Other tissues, NSGCT'
telomerehunter_df[telomerehunter_df$diagnosis_plot == 'Seminoma' & telomerehunter_df$histology == 'GCNIS',]$diagnosis_plot_reannot = 'GCNIS, pure seminoma'
telomerehunter_df[telomerehunter_df$diagnosis_plot == 'Post-pubertal mixed germ cell tumour' & telomerehunter_df$histology == 'GCNIS',]$diagnosis_plot_reannot = 'GCNIS, NSGCT'

telomerehunter_df$diagnosis_plot_reannot = factor(telomerehunter_df$diagnosis_plot_reannot, levels = c('Normal', 'GCNIS, pure seminoma', 'GCNIS, NSGCT', 'Seminoma, pure seminoma', 'Dysgerminoma/Seminoma, NSGCT', 'Pre-pubertal yolk sac tumour', 'Other tissues, NSGCT'))
levels(telomerehunter_df$diagnosis_plot_reannot) = c('Normal', 'GCNIS, pure seminoma', 'GCNIS, NSGCT', 'Seminoma, pure seminoma', 'Dysgerminoma/Seminoma, NSGCT', 'Prepubertal YST', 'Other tissues, NSGCT')
telomerehunter_df$site = as.factor(telomerehunter_df$site)
levels(telomerehunter_df$site) = c('Blood', 'Ovary', 'Testis')

pdf('/lustre/scratch119/casm/team294rr/to3/testes/tumour/manuscript/figures/gct_telomerecontent_by_tissue_by_tumour.pdf', width = 8, height = 6, useDingbats = F)
ggplot(telomerehunter_df) + 
  geom_quasirandom(mapping = aes(x = diagnosis_plot_reannot, y = tel_content_partial, col = diagnosis_plot_reannot), cex = 2) + 
  geom_pointrange(mapping = aes(x = diagnosis_plot_reannot, y = tel_content_partial),
                  stat = "summary", shape=19, 
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill="black", cex = 1) +
  facet_grid(. ~ site, scales = 'free_x', space = 'free_x') + theme_pubr(legend = 'right') + 
  labs(y = 'Normalised telomere content', x = '') + theme(axis.text.x = element_text(angle = 315, vjust = 1, hjust = 0), strip.background = element_blank(), strip.text.x = element_text(size = 12)) + scale_color_manual(values = c('#999999', as.vector(get_palette(palette = "Dark2", 6))), name = 'Category')
dev.off()

#get stats to place on plot
wilcox.test(telomerehunter_df[telomerehunter_df$diagnosis_plot_reannot == 'Dysgerminoma/Seminoma, NSGCT' & telomerehunter_df$site == 'Ovary', ]$tel_content_partial, telomerehunter_df[telomerehunter_df$diagnosis_plot_reannot == 'Dysgerminoma/Seminoma, NSGCT' & telomerehunter_df$site == 'Testis', ]$tel_content_partial, alternative = 'two.sided') #W = 0, p-value = 0.009524 **
wilcox.test(telomerehunter_df[telomerehunter_df$diagnosis_plot_reannot == 'Other tissues, NSGCT' & telomerehunter_df$site == 'Ovary', ]$tel_content_partial, telomerehunter_df[telomerehunter_df$diagnosis_plot_reannot == 'Other tissues, NSGCT' & telomerehunter_df$site == 'Testis', ]$tel_content_partial, alternative = 'two.sided') #W = 0, p-value = 3.052e-08 ***

#get number of genomes that support each category
aggregate(tel_content_partial ~ diagnosis_plot_reannot + site, telomerehunter_df, length)