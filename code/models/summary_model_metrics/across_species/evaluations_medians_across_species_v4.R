##code for obtaining global metrics of model evaluation across species

#load evaluation metrics for all species
evaluations = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/evaluations/medians_evaluations_v4.csv", sep=",", header=TRUE)
str(evaluations)

#caculate the mean and sd by columns considering all columns except the specie name column. This will include model metrics considering all species
means_across_sp = apply(evaluations[,-1], 2, mean)
sd_across_species = apply(evaluations[,-1], 2, sd)

#see differences auc, kappa and tss between models


#kappa
means_across_sp["glm_kappa_mean"]; means_across_sp["gam_kappa_mean"]; means_across_sp["rf_kappa_mean"]
sd_across_species["glm_kappa_mean"]; sd_across_species["gam_kappa_mean"]; sd_across_species["rf_kappa_mean"]
wilcox.test(evaluations$glm_kappa_mean, evaluations$gam_kappa_mean)
wilcox.test(evaluations$glm_kappa_mean, evaluations$rf_kappa_mean)
wilcox.test(evaluations$gam_kappa_mean, evaluations$rf_kappa_mean)

#tss
means_across_sp["glm_tss_mean"]; means_across_sp["gam_tss_mean"]; means_across_sp["rf_tss_mean"]
sd_across_species["glm_tss_mean"]; sd_across_species["gam_tss_mean"]; sd_across_species["rf_tss_mean"]
wilcox.test(evaluations$glm_tss_mean, evaluations$gam_tss_mean)
wilcox.test(evaluations$glm_tss_mean, evaluations$rf_tss_mean)
wilcox.test(evaluations$gam_tss_mean, evaluations$rf_tss_mean)

#auc (we dont have auc for RF as its predictions were obteind direclty as binary, so no thresholds could be calculated)
means_across_sp["glm_auc_mean"]; means_across_sp["gam_auc_mean"]
sd_across_species["glm_auc_mean"]; sd_across_species["gam_auc_mean"]
wilcox.test(evaluations$glm_auc_mean, evaluations$gam_auc_mean)

#see summary of these metrics between different models
summary(evaluations$glm_kappa_mean)
summary(evaluations$gam_kappa_mean)
summary(evaluations$rf_kappa_mean)

summary(evaluations$glm_tss_mean)
summary(evaluations$gam_tss_mean)
summary(evaluations$rf_tss_mean)

summary(evaluations$glm_auc_mean)
summary(evaluations$gam_auc_mean)