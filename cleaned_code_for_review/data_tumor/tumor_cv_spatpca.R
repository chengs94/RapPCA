library(SpatialPCA)
library(mgcv)
library(pbmcapply)

load("Tumor_data.RData")
load("tumor_data_norm.RData")

# make sure we use the same configurations (e.g. CV fold assignment) for all PCA methods
rappca_ft = readRDS("tumor_cv_ft.rda")
cv_idx = rappca_ft$cv_idx
K = length(cv_idx)
r = ncol(rappca_ft$utildes)
rm(rappca_ft)
m = 2

# predictions made by SpatialPCA's own function
pred_spatpca = function(k){
  tst_idx = cv_idx[[k]]
  uhats_spatpca = uhats_spatpca2 = matrix(NA, nrow = length(tst_idx), ncol = r)
  
  object <- new(
    Class = "SpatialPCA", counts = rawcount[,-tst_idx], project = "SpatialPCA",
    location = coords[-tst_idx,],  normalized_expr = t(y[-tst_idx,])
  )
  ST = SpatialPCA_buildKernel(object, kerneltype="gaussian", bandwidthtype="SJ")
  ST = SpatialPCA_EstimateLoading(ST,fast=FALSE,SpatialPCnum=r)
  ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)
  
  utildes = t(ST@SpatialPCs)
  x1train = coords[-tst_idx,1]
  x2train = coords[-tst_idx,2]
  for (j in 1:ncol(uhats_spatpca)){
    mod <- mgcv::gam(utildes[,j] ~ s(x1train, x2train, bs="tp", k=nrow(utildes), m=m, fx=FALSE))
    uhats_spatpca[,j] = predict(mod, data.frame(x1train = coords[tst_idx,1], 
                                                       x2train = coords[tst_idx,2]))
    message(j)
  }
  spatpca_highres = SpatialPCA_highresolution(ST, platform="ST",newlocation=coords[tst_idx,])
  uhats_spatpca2 = t(spatpca_highres@highPCs)
  
  return(list(uhats_sp_tprs = uhats_spatpca, uhats_sp_default = uhats_spatpca2))
}

if (!file.exists("tumor_cv_spatpca.rda")){
  cv_rslt_sp = pbmclapply(1:K, pred_spatpca, mc.cores = 5)
  
  saveRDS(cv_rslt_sp, file = "tumor_cv_spatpca.rda")
  save.image("tumor_cv_spatpca_full.RData")
} else {
  cv_rslt_sp = readRDS("tumor_cv_spatpca.rda")
}

mspe_sp = msre_sp = tmse_sp = mspe_sp2 = tmse_sp2 = rep(NA, K)
mse_pc_sp = matrix(NA, nrow = K, ncol = r)
mse_pc_sp2 = matrix(NA, nrow = K, ncol = r)
utildes_sp = matrix(NA, nrow = nrow(y), ncol = r)
for (k in 1:K){
  tst_idx = cv_idx[[k]]
  uhats_spatpca = uhats_spatpca2 = matrix(NA, nrow = length(tst_idx), ncol = r)
  
  object <- new(
    Class = "SpatialPCA", counts = rawcount[,-tst_idx], project = "SpatialPCA",
    location = coords[-tst_idx,],  normalized_expr = t(y[-tst_idx,])
  )
  ST = SpatialPCA_buildKernel(object, kerneltype="gaussian", bandwidthtype="SJ")
  ST = SpatialPCA_EstimateLoading(ST,fast=FALSE,SpatialPCnum=r)
  ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)
  
  XB = colMeans(y[-tst_idx,])
  y_centered = scale(y[tst_idx,], center = XB, scale = FALSE)
  utildes.tst = y_centered %*% ST@W
  uhats_sp_tprs = cv_rslt_sp[[k]]$uhats_sp_tprs
  uhats_sp_default = cv_rslt_sp[[k]]$uhats_sp_default
  
  mse_pc_sp[k,] = colMeans((utildes.tst - uhats_sp_tprs)^2)
  mse_pc_sp2[k,] = colMeans((utildes.tst - uhats_sp_default)^2)
  tmse_sp[k] = sum((y_centered - uhats_sp_tprs%*%t(ST@W))^2) / length(tst_idx) 
  tmse_sp2[k] = sum((y_centered - uhats_sp_default%*%t(ST@W))^2) / length(tst_idx) 
  mspe_sp[k] = sum(((uhats_sp_tprs - utildes.tst) %*% t(ST@W))^2) / length(tst_idx)
  mspe_sp2[k] = sum(((uhats_sp_default - utildes.tst) %*% t(ST@W))^2) / length(tst_idx)
  msre_sp[k] = sum((y_centered - utildes.tst %*% t(ST@W))^2) / length(tst_idx)
  
  utildes_sp[tst_idx,] = utildes.tst
}

save(mspe_sp, mspe_sp2, msre_sp, tmse_sp, tmse_sp2,
     mse_pc_sp, mse_pc_sp2, utildes_sp,
     file = "tumor_cv_spatpca_metrics.RData")
