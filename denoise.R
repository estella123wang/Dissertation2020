pro_tec = readRDS("~/hda_tds_Mechanomics_project/Proteins/Proteins_technical_covariates.rds")
covars = readRDS("~/dissertation/data/Covariates1.rds") 



library(plyr)
require(dplyr)
pro_tec$plate <- pro_tec$Plate.ID

pro_tec <- 
    pro_tec %>%
        mutate(plate = revalue(plate, c('20181212-001_SP190747_INFI'=1, 
                                        '20181212-002_SP190748_INFI'=2,
                                        '20181212-003_SP190749_INFI'=3,
                                        '20181212-004_SP190750_INFI'=4,
                                        '20181212-005_SP190751_INFI'=5,
                                        '20181212-006_SP190752_INFI'=6,
                                        '20181212-007_SP190753_INFI'=7,
                                        '20181212-008_SP190754_INFI'=8)))

covars<-merge(pro_tec,covars)
covars$subtype[is.na(covars$subtype)] <- 0



proteins =readRDS("~/dissertation/data/Proteins.rds")

proteins <- na.omit(proteins)

comparison <- compare(rownames(proteins),rownames(covars),allowAll=TRUE)
paired_id<-comparison$tM

proteins = proteins[paired_id,]


suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(RColorBrewer))



denoised = NULL 
Beta_pooled = NULL 
pvalue_pooled = NULL


f0 = "proteins[,k] ~ age.sample + gender  + centre + (1 | plate)"
f1 = paste(f0, "+ case")
as.formula(f1)




for (k in seq(1:ncol(proteins))) {
  # print(k)
  model = lmer(as.formula(f1), data = covars, REML = FALSE,
               control = lmerControl(check.conv.singular = .makeCC(action = "ignore",tol = 1e-04)))
  model0 = lmer(as.formula(f0), data = covars, REML = FALSE,control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  pvalue_pooled = c(pvalue_pooled, anova(model, model0)$"Pr(>Chisq)"[2])
  beta = fixef(model)[c("(Intercept)", "case1")]
  Beta_pooled = c(Beta_pooled, fixef(model)["case1"])
  X = cbind(rep(1, length(covars$case)), as.numeric(covars$case))
  denoised = cbind(denoised, (X %*% beta + resid(model)))
}




colnames(denoised) = colnames(proteins) 
rownames(denoised) = rownames(proteins) 
saveRDS(denoised, "~/dissertation/data/test1_denoised.rds")


