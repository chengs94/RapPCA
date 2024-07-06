setwd("/home/users/chengsi/Desktop/multi_pollutants")

rm(list = ls())

source("rpPCA_profile.R")
source("pred_pca.R")
source("pca_funcs.R")

# temp = one_sim(dim_x = 10, n_pc = 2, n_noise = 3, dim_y = 15, seed = 1, d = 1, r = 3)
# case1 = mclapply(X = 1:100, FUN = one_sim, mc.cores = detectCores() - 1,
#                  dim_x = 10, n_pc = 3, n_noise = 3, dim_y = 15, r = 3, #d = 2,
#                  cov_range = .5, cov_sill = .5, cov_nugget = .5, y_sig = .1, 
#                  gms = c(.05, .1, .2, .5, 1, 2, 5),
#                  lambdas = c(.05, .1, .2, .5, 1, 2, 5),
#                  lambda.ratios = c(.05, .1, .2, .5, 1, 2, 5),
#                  krige = TRUE, decay_v = FALSE,
#                  linear = TRUE, mean_function = sin
# )
# saveRDS(case1, file = "sim_case1_updated.rda")
# 
# case2 = mclapply(X = 1:100, FUN = one_sim, mc.cores = detectCores() - 1,
#                  dim_x = 10, n_pc = 3, n_noise = 3, dim_y = 15, r = 3, #d = 2,
#                  cov_range = .5, cov_sill = .5, cov_nugget = .5, y_sig = .1, 
#                  gms = c(.05, .1, .2, .5, 1, 2, 5),
#                  lambdas = c(.05, .1, .2, .5, 1, 2, 5),
#                  lambda.ratios = c(.05, .1, .2, .5, 1, 2, 5),
#                  krige = TRUE, decay_v = TRUE,
#                  linear = TRUE, mean_function = sin
# )
# saveRDS(case2, file = "sim_case2_updated.rda")

fun_int = function(x){
  odd.idx = seq(1, ncol(x), length.out = floor(ncol(x)/2))
  even.idx = seq(2, ncol(x), length.out = floor(ncol(x)/2))
  temp = x[,odd.idx]*x[,even.idx]
  # c = rep(1, ncol(temp))
  # c[seq(2,ncol(temp),2)] = -1
  b1 = matrix(runif(ncol(x)*3, -1, 1), nrow = ncol(x))
  b2 = matrix(runif(floor(ncol(x)/2)*3, -1, 1), nrow = floor(ncol(x)/2))
  return(2*temp%*%b2 + x^2%*%b1)
}
case3 = mclapply(X = 1:100, FUN = one_sim, mc.cores = detectCores() - 1,
                 dim_x = 10, n_pc = 3, n_noise = 3, dim_y = 15, r = 3, # d = 2,
                 cov_range = .5, cov_sill = .5, cov_nugget = .5, y_sig = .1, 
                 gms = c(.05, .1, .2, .5, 1, 2, 5),
                 lambdas = c(.05, .1, .2, .5, 1, 2, 5),
                 lambda.ratios = c(.05, .1, .2, .5, 1, 2, 5),
                 krige = TRUE, decay_v = FALSE,
                 linear = FALSE, mean_function = fun_int# function(x) x^4 / 25
)
saveRDS(case3, file = "sim_case4_interaction2.rda")
# save(case1, case2, case3, file = "sim_case123_updated.RData")
save.image("case4_interaction2_image.RData")

# lowdim_pred = mclapply(X = 1:100, FUN = one_sim, mc.cores = detectCores() - 1,
#                         dim_x = 10, n_pc = 4, n_noise = 1, dim_y = 15, r = 3, #d = 2,
#                         cov_range = .5, cov_sill = .1, cov_nugget = .1, y_sig = .15, 
#                         gms = c(.05, .1, .2, .5, 1, 2, 5),
#                         lambdas = c(.05, .1, .2, .5, 1, 2, 5),
#                         lambda.ratios = c(.05, .1, .2, .5, 1, 2, 5),
#                         # lambda2s = c(.05, .1, .2, .5, 1, 2),
#                         # nvars = c(2, 5, 8, 10, 15),
#                         linear = TRUE
#                         # gms = c(seq(.05,.2,.05), seq(.5,2,.5), seq(4,10,2)),
#                         # lambdas = c(seq(.05,.2,.05), seq(.5,2,.5), 5, 8),
#                         # lambda.ratios = c(seq(.05,.2,.05), seq(.5,2,.5), 5, 8)
# )
# save.image("fullsim_lowdim_cache.RData")
# 
# lowdim_unpred = mclapply(X = 1:100, FUN = one_sim, mc.cores = detectCores() - 1,
#                           dim_x = 10, n_pc = 1, n_noise = 4, dim_y = 15, r = 3, #d = 2,
#                           cov_range = .5, cov_sill = .05, cov_nugget = .3, y_sig = .5, 
#                           gms = c(.05, .1, .2, .5, 1, 2, 5),
#                           lambdas = c(.05, .1, .2, .5, 1, 2, 5),
#                           lambda.ratios = c(.05, .1, .2, .5, 1, 2, 5),
#                           # lambda2s = c(.05, .1, .2, .5, 1, 2),
#                           # nvars = c(2, 5, 8, 10, 15),
#                           linear = TRUE
#                           # gms = c(seq(.05,.2,.05), seq(.5,2,.5), seq(4,10,2)),
#                           # lambdas = c(seq(.05,.2,.05), seq(.5,2,.5), 5, 8),
#                           # lambda.ratios = c(seq(.05,.2,.05), seq(.5,2,.5), 5, 8)
# )
# save.image("fullsim_lowdim.RData")
# 
# highdim_pred = mclapply(X = 1:100, FUN = one_sim, mc.cores = detectCores() - 1,
#                        dim_x = 100, n_pc = 4, n_noise = 1, dim_y = 15, r = 3, #d = 2,
#                        cov_range = .5, cov_sill = .1, cov_nugget = .1, y_sig = .15, 
#                        gms = c(.05, .1, .2, .5, 1, 2, 5),
#                        lambdas = c(.05, .1, .2, .5, 1, 2, 5),
#                        lambda.ratios = c(.05, .1, .2, .5, 1, 2, 5),
#                        # lambda2s = c(.05, .1, .2, .5, 1, 2),
#                        # nvars = c(2, 5, 8, 10, 15),
#                        linear = FALSE, mean_function = function(x){x^2}
#                        # gms = c(seq(.05,.2,.05), seq(.5,2,.5), seq(4,10,2)),
#                        # lambdas = c(seq(.05,.2,.05), seq(.5,2,.5), 5, 8),
#                        # lambda.ratios = c(seq(.05,.2,.05), seq(.5,2,.5), 5, 8)
#                        )
# save.image("fullsim_highdim_cache.RData")
# 
# highdim_unpred = mclapply(X = 1:100, FUN = one_sim, mc.cores = detectCores() - 1,
#                        dim_x = 100, n_pc = 1, n_noise = 4, dim_y = 15, r = 3, #d = 2,
#                        cov_range = .5, cov_sill = .05, cov_nugget = .3, y_sig = .5, 
#                        gms = c(.05, .1, .2, .5, 1, 2, 5),
#                        lambdas = c(.05, .1, .2, .5, 1, 2, 5),
#                        lambda.ratios = c(.05, .1, .2, .5, 1, 2, 5),
#                        # lambda2s = c(.05, .1, .2, .5, 1, 2),
#                        # nvars = c(2, 5, 8, 10, 15),
#                        linear = FALSE, mean_function = function(x){x^2}
#                        # gms = c(seq(.05,.2,.05), seq(.5,2,.5), seq(4,10,2)),
#                        # lambdas = c(seq(.05,.2,.05), seq(.5,2,.5), 5, 8),
#                        # lambda.ratios = c(seq(.05,.2,.05), seq(.5,2,.5), 5, 8)
#                        )
# save(lowdim_pred, lowdim_unpred, highdim_pred, highdim_unpred,
#      file = "fullsim.RData")
# save.image("fullsim_image.RData")
