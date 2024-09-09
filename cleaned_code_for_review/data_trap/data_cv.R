# run cross-validation for the Seattle TRAP data

library(dplyr)

source("rpPCA_profile.R")
source("pred_pca.R")
source("pca_funcs.R")

load("pca_data.RData")

coords = cbind(dat_x$lambert_x, dat_x$lambert_y)
x = dat_x[,-(1:6)]
y = dat_y %>% 
  mutate_at(vars(c("bc_uv")), function(x) x - dat_y$bc_ir) %>%
  mutate(ufp_pnc.noscreen = ufp_pnc.noscreen - ufp_pnc.screen) %>%
  mutate_at(vars(-c("location", "longitude", "latitude", 
                    "no2", "ufp_pmdisc.size")),
            function(x) x / dat_y$no2) %>% 
  select(-c(pm2.5, co, co2, ufp_total.conc,
            location, longitude, latitude))

cv_rslt = cv_pca(x, y, coords, r = 3, 
                 gms = c(.1, .5, seq(1,9,2), 12, 15, 20, 30, 50, 70, 90),
                 lambdas = c(.05, .1, .2, .5, 1, 2, 5),
                 lambda.ratios = c(.05, .1, .2, .5, 1, 2, 5),
                 d = 1, m = 2, x.dim.reduction = FALSE, 
                 x.dim.reduction2 = FALSE, x.r = 15, K = 10)

save(cv_rslt, file = "cv_rslt_u_w_int_v2.rda")
save.image("cv_rslt_u_w_int_v2_full.RData")