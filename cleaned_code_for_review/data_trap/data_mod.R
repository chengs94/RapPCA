# make predictions on the Seattle TRAP data
# at a finer grid of locations

library(dplyr)

source("pred_pca.R")
# source("pca_funcs.R")

# load("data/cv_rslt_u_full.RData")
load("data/cv_rslt_u_w_int_v2_full.RData")
source("rpPCA_profile.R")
cov_grid = readRDS("cov_grid.rda")
xtst = cov_grid[,-(1:6)]
coords.tst = cbind(cov_grid$lambert_x, cov_grid$lambert_y)

set.seed(233)

r = 3 # number of PCs
x.r = 15
d = 1
m = 2
p = ncol(y)
pars = cv_rslt$par_best

x.scaled = scale(x)
xtst = xtst[, colnames(x)] %>% 
  scale(center = attr(x.scaled, 'scaled:center'), 
        scale = attr(x.scaled, 'scaled:scale'))
n.tst = nrow(xtst)
coords = scale(coords)
coords.tst = scale(coords.tst, 
                   center = attr(coords, 'scaled:center'),
                   scale = attr(coords, 'scaled:scale'))
x.pca = prcomp(x, center = TRUE, scale. = TRUE)
.x = x.pca$x[,1:x.r]
Btrn = smooth.construct(s(coord1, coord2, bs="tp", k=4, m=m),
                        data = list(coord1 = coords[,1], coord2 = coords[,2]),
                        knots = NULL)
xx = Btrn$X[,-which(apply(Btrn$X, 2, sd) < .01)]
x.scaled2 = scale(cbind(.x, xx))
b.init = runif(ncol(x.scaled2), -1, 1)

x1train = coords[,1]
x2train = coords[,2]

y.scaled0 = y.scaled = y.scaled.pca = y.scaled.pred = scale(y)
loadings = uhats = utildes = y.res = NULL
loadings.pca = uhats.pca = utildes.pca = y.res.pca = NULL
loadings.pred = uhats.pred = utildes.pred = y.res.pred = NULL
.rslt = .rslt.pca = list(v = rep(0, ncol(y)))
.v = rep(0, ncol(y))
.r = 0
# tmse = mspe = msre = 0

while (r > 0){
  .r = .r + 1
  r = r-1
  y.scaled = y.scaled - y.scaled %*% .rslt$v %*% t(.rslt$v)
  y.scaled.pca = y.scaled.pca - y.scaled.pca %*% .rslt.pca$v %*% t(.rslt.pca$v)
  y.scaled.pred = y.scaled.pred - y.scaled.pred %*% .v %*% t(.v)
  y.svd = svd(t(y.scaled)) 
  dmat = diag(y.svd$d)
  U = y.svd$u
  V = y.svd$v
  y.svd1 = svd(t(y.scaled.pca)) 
  dmat1 = diag(y.svd1$d)
  U1 = y.svd1$u
  V1 = y.svd1$v
  
  .rslt = rep_predPCA(gm = pars[.r,1], lambda = pars[.r,2], lambda.ratio = pars[.r,3],
                          xtrn = x.scaled, xtst = xtst, 
                          ytrn = y.scaled, ytst = matrix(0, nrow = n.tst, ncol = p), 
                          coords.trn = coords, coords.tst = coords.tst,
                          d = d, m = m, dmat = dmat, U = U, V = V)
  .rslt.pca = rep_predPCA(gm = -1, lambda = 0, lambda.ratio = 0,
                          xtrn = x.scaled, xtst = xtst, 
                          ytrn = y.scaled.pca, ytst = matrix(0, nrow = n.tst, ncol = p), 
                          coords.trn = coords, coords.tst = coords.tst,
                          d = d, m = m, dmat = dmat1, U = U1, V = V1)
  .rslt.pred = space.sparse.pca(Z = x.scaled2, X = y.scaled.pred,
                           beta.init = b.init,
                           v.init = runif(ncol(y.scaled.pred), -1, 1),
                           lambda2 = 1, varnum1 = 1,
                           niter = 500, err = 1e-4)
  
  y.res = cbind(y.res, y.scaled)
  loadings = cbind(loadings, .rslt$v)
  uhats = cbind(uhats, .rslt$uhat.tst)
  utildes = cbind(utildes, .rslt$utilde.tst)
  
  y.res.pca = cbind(y.res.pca, y.scaled.pca)
  loadings.pca = cbind(loadings.pca, .rslt.pca$v)
  uhats.pca = cbind(uhats.pca, .rslt.pca$uhat.tst)
  utildes.pca = cbind(utildes.pca, .rslt.pca$utilde.tst)
  
  # .v = .rslt$v/norm(.rslt$v, '2') #!!!! wrong
  .v = .rslt.pred$v/norm(.rslt.pred$v, '2')
  .uhat.trn = y.scaled.pred%*%.v
  rf <- randomForest(scale(x), .uhat.trn,
                     xtest = xtst, nodesize=5, importance = FALSE)
  mod <- mgcv::gam((.uhat.trn - rf$predicted)~s(x1train, x2train,
                                                bs="tp", k=160, m=m, fx=FALSE))
  
  .uhat.tst <- rf$test$predicted + predict(mod,
                                           data.frame(x1train = coords.tst[,1],
                                                      x2train = coords.tst[,2]))
  .utilde = y.scaled.pred %*%.v
  
  y.res.pred = cbind(y.res.pred, y.scaled.pred)
  loadings.pred = cbind(loadings.pred, .v)
  utildes.pred = cbind(utildes.pred, .utilde)
  uhats.pred = cbind(uhats.pred, .uhat.tst)
}

save(list = c(grep("y.res", ls(), value=TRUE),
              grep("loadings", ls(), value = TRUE),
              grep("uhats", ls(), value = TRUE)), 
     file = "mod_grid_u_v2_corrected.RData")
