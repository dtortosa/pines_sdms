##code for obtaining global metrics of model evaluation across species

#load evaluation metrics for all species
evaluations = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/evaluations/medians_evaluations_v2.csv", sep=",", header=TRUE)
str(evaluations)

#caculate the median and sd by columns considering all columns except the specie name column. This will include model metrics considering all species
medians_across_sp = apply(evaluations[,-1], 2, median)
sd_across_species = apply(evaluations[,-1], 2, sd)

#see differences auc, kappa and tss between models


#kappa
medians_across_sp["glm_kappa_median"]; medians_across_sp["gam_kappa_median"]; medians_across_sp["rf_kappa_median"]
sd_across_species["glm_kappa_median"]; sd_across_species["gam_kappa_median"]; sd_across_species["rf_kappa_median"]
wilcox.test(evaluations$glm_kappa_median, evaluations$gam_kappa_median)
wilcox.test(evaluations$glm_kappa_median, evaluations$rf_kappa_median)
wilcox.test(evaluations$gam_kappa_median, evaluations$rf_kappa_median)

#tss
medians_across_sp["glm_tss_median"]; medians_across_sp["gam_tss_median"]; medians_across_sp["rf_tss_median"]
sd_across_species["glm_tss_median"]; sd_across_species["gam_tss_median"]; sd_across_species["rf_tss_median"]
wilcox.test(evaluations$glm_tss_median, evaluations$gam_tss_median)
wilcox.test(evaluations$glm_tss_median, evaluations$rf_tss_median)
wilcox.test(evaluations$gam_tss_median, evaluations$rf_tss_median)

#auc
medians_across_sp["glm_auc_median"]; medians_across_sp["gam_auc_median"]; medians_across_sp["rf_auc_median"]
sd_across_species["glm_auc_median"]; sd_across_species["gam_auc_median"]; sd_across_species["rf_auc_median"]
wilcox.test(evaluations$glm_auc_median, evaluations$gam_auc_median)
wilcox.test(evaluations$glm_auc_median, evaluations$rf_auc_median)
wilcox.test(evaluations$gam_auc_median, evaluations$rf_auc_median)


#see summary of these metrics between different models
summary(evaluations$glm_kappa_median)
summary(evaluations$gam_kappa_median)
summary(evaluations$rf_kappa_median)

summary(evaluations$glm_tss_median)
summary(evaluations$gam_tss_median)
summary(evaluations$rf_tss_median)

summary(evaluations$glm_auc_median)
summary(evaluations$gam_auc_median)
summary(evaluations$rf_auc_median)