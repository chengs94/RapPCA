# verifying the optimality of our solution
library(pbmcapply)

source("rpPCA_profile.R")
load("rslts/sim_case123_updated.RData")
dat = case1[[1]]$data
rm(case1, case2, case3)

y = scale(dat$y)
x = dat$x
x.scaled = scale(x)
coords = scale(dat$coords)

# set.seed(233)
# tst.ind = sample(1:nrow(x), floor(nrow(x) * .2), replace = FALSE)

y.svd = svd(t(y)) 
dm = y.svd$d[1]
dmat = diag(y.svd$d)
U = y.svd$u
V = y.svd$v

# pars = expand.grid(lambda = exp(seq(-3,5,.5)),
#                    lambda.ratio = exp(seq(-2,2,.5)),
#                    gm = c(0, exp(seq(-4,7,.5))))
pars = expand.grid(lambda = c(seq(.05,.2,.05), seq(.5,2,.5), seq(4,16,2)),
                   lambda.ratio = c(seq(.05,.2,.05), seq(.5,2,.5), seq(4,16,2)),
                   gm = c(0, seq(.05,.2,.05), seq(.5,2,.5), seq(4,16,2)))

rslts = pbmcmapply(FUN = rep_predPCA, mc.cores = detectCores() - 2,
                   lambda = pars$lambda, lambda.ratio = pars$lambda.ratio, gm = pars$gm,
                   MoreArgs = list(xtrn = x.scaled, xtst = x.scaled,
                                   ytrn = y, ytst = y,
                                   coords.trn = coords, coords.tst = coords,
                                   dmat = dmat, U = U, V = V, return.condn = FALSE,
                                   delta = 0.05, d = 1, m = 2))

save(rslts, dat, pars, file = "rslts/case1_suppl.RData")
save.image("rslts/case1_suppl_image.RData")
