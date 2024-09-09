# run cross-validation for the spatial transcriptomics data

library(dplyr)

# rm(list = ls())

source("rpPCA_profile.R")
source("pred_pca.R")
source("pca_funcs.R")

load("tumor_data_norm.RData")

r = 20
cv_rslt = cv_pca(x = NULL, y = y, coords = coords, r = r,
                 gms = c(.05, seq(.1,.9,.2), 1, 1.5, 2, 5, 10, 50),
                 lambdas = c(.05, .1, .2, .5, 1, 2),
                 lambda.ratios = 1,
                 d = 1, m = 2, x.dim.reduction = FALSE,
                 x.dim.reduction2 = FALSE, x.r = 15, K = 10,
                 compact = FALSE)

saveRDS(cv_rslt, file = "tumor_cv_ft.rda")
save.image("tumor_cv_ft_full.RData")