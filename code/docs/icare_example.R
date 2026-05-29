library(iCARE)
data("bc_data", package="iCARE")
set.seed(50)
res_snps_miss = computeAbsoluteRisk(model.snp.info = bc_72_snps,
                                    model.disease.incidence.rates = bc_inc,
                                    model.competing.incidence.rates = mort_inc,
                                    apply.age.start = 50, apply.age.interval.length = 30,
                                    return.refs.risk = TRUE)
summary(res_snps_miss$refs.risk)
res_snps_dat = computeAbsoluteRisk(model.snp.info = bc_72_snps,
                                   model.disease.incidence.rates = bc_inc,
                                   model.competing.incidence.rates = mort_inc,
                                   apply.age.start = 50, apply.age.interval.length = 30,
                                   apply.snp.profile = new_snp_prof,
                                   return.refs.risk = TRUE)
names(res_snps_dat)
plot(density(res_snps_dat$refs.risk),
     xlim = c(0.04,0.18), xlab = "Absolute Risk of Breast Cancer",
     main = "Referent SNP-only Risk Distribution: Ages 50-80 years")
abline(v = res_snps_dat$risk, col = "red")
legend("topright", legend = "New profiles", col = "red", lwd = 1)

############################################################

res_covs_snps = computeAbsoluteRisk(model.formula = bc_model_formula,
                                    + model.cov.info = bc_model_cov_info,
                                    + model.snp.info = bc_72_snps,
                                    + model.log.RR = bc_model_log_or,
                                    + model.ref.dataset = ref_cov_dat,
                                    + model.disease.incidence.rates = bc_inc,
                                    + model.competing.incidence.rates = mort_inc,
                                    + model.bin.fh.name = "famhist",
                                    + apply.age.start = 50,
                                    + apply.age.interval.length = 30,
                                    + apply.cov.profile = new_cov_prof,
                                    + apply.snp.profile = new_snp_prof,
                                    + return.refs.risk = TRUE)