library(SpatialPCA)
library(mgcv)


set.seed(233)

load("tumor_cv_ft_full.RData")
coords.tst = coords
pars = cv_rslt$par_best
d = 1
m = 2
r = 20
p = ncol(y)
n = nrow(y)
xtst = x.scaled = NULL

# ST = CreateSpatialPCAObject(counts=rawcount, location=location, project = "SpatialPCA",
#                             gene.type="spatial",sparkversion="spark", gene.number=3000,
#                             customGenelist=NULL,min.loctions = 20, min.features=20)
# y = t(ST@normalized_expr)
# coords = coords.tst = scale(location)
# save(y, coords, file = "tumor_data_norm.RData")

x1train = coords[,1]
x2train = coords[,2]
Btrn = smooth.construct(s(coord1, coord2, bs="tp", k=10, m=m),
                        data = list(coord1 = coords[,1], coord2 = coords[,2]),
                        knots = NULL)
xx = Btrn$X[,-which(apply(Btrn$X, 2, sd) < .01)]
x.scaled2 = cbind(1, scale(xx))
b.init = runif(ncol(x.scaled2), -1, 1)

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
  
  .rslt = rep_predPCA(gm = pars[.r,1], lambda = pars[.r,2], lambda.ratio = 1,
                      xtrn = x.scaled, xtst = xtst, 
                      ytrn = y.scaled, ytst = y.scaled, 
                      coords.trn = coords, coords.tst = coords.tst,
                      d = d, m = m, dmat = dmat, U = U, V = V)
  .rslt.pca = rep_predPCA(gm = -1, lambda = 0, lambda.ratio = 0,
                          xtrn = x.scaled, xtst = xtst, 
                          ytrn = y.scaled.pca, ytst = y.scaled, 
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
  
  .v = .rslt.pred$v/norm(.rslt.pred$v, '2')
  .uhat.trn = y.scaled.pred%*%.v
  
  mod <- mgcv::gam((.uhat.trn)~s(x1train, x2train, bs="tp", k=160, m=m, fx=FALSE))
  .uhat.tst <- predict(mod, data.frame(x1train = coords.tst[,1],
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
     file = "STmod_rslt_ft.RData")
save.image("STmod_rslt_ft_full.RData")
