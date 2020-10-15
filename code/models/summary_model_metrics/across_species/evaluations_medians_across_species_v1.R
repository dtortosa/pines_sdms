evaluations = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/evaluations/medians_evaluations.csv", sep=",", header=TRUE)
str(evaluations)


medians_across_sp = apply(evaluations[,-1], 2, median)
sd_across_species = apply(evaluations[,-1], 2, sd)


medians_across_sp["glm_kappa_median"]; medians_across_sp["gam_kappa_median"]
sd_across_species["glm_kappa_median"]; sd_across_species["gam_kappa_median"]
wilcox.test(evaluations$glm_kappa_median, evaluations$gam_kappa_median)


summary(evaluations$glm_kappa_median)
summary(evaluations$gam_kappa_median)

summary(evaluations$glm_tss_median)
summary(evaluations$gam_tss_median)

summary(evaluations$rf_oob_median)

