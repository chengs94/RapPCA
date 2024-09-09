# auxiliary functions for RapPCA

# generate data for simulation
sim_y = function(dim_x, n_pc, n_noise, dim_y, 
                 cov_range = 1, cov_sill = .5, cov_nugget = .5, 
                 y_sig = .5, n = 200, seed = 1, 
                 mean_function = NULL, linear = TRUE,
                 krige = TRUE, decay_v = FALSE){
  
  set.seed(seed)
  x = matrix(runif(n*dim_x, -1, 1), nrow = n)
  b = matrix(runif(dim_x*n_pc, -1, 1), nrow = dim_x)
  if (linear){
    f_xb = x%*%b
  } else{
    # f_xb = mean_function(x%*%b)
    f_xb = mean_function(x) # changed for the interaction terms only
  }
  
  coord_x = runif(n, 0, 1)
  coord_y = runif(n, 0, 1)
  if (krige){
    d_sq = outer(coord_x, coord_x, '-')^2 + outer(coord_y, coord_y, '-')^2
    pc_cov = cov_sill*(exp(-d_sq/cov_range)) + cov_nugget*diag(n)
    pcs = cbind(f_xb,
                t(mvtnorm::rmvnorm(n_noise, sigma = pc_cov)))
    bb = NULL
  } else {
    Btrn = smooth.construct(s(coord1, coord2, bs="tp", k=4, m=2),
                            data = list(coord1 = coord_x, coord2 = coord_y),
                            knots = NULL)
    bb = matrix(runif(4*(n_noise), 0, 1), nrow = 4)
    pcs = cbind(f_xb, 
                cov_sill* Btrn$X %*% bb +
                  matrix(rnorm(n*(n_noise), 0, cov_nugget), nrow = n))
    pc_cov = NULL
  }
  
  pcs[,(n_pc+1):(n_pc+n_noise)] = pcs[,(n_pc+1):(n_pc+n_noise)] *
    mean(apply(as.matrix(pcs[,1:n_pc]), 2, sd))/mean(apply(as.matrix(pcs[,(n_pc+1):(n_pc+n_noise)]), 2, sd))
  
  v_raw = matrix(runif((n_pc+n_noise)*dim_y, -1, 1), 
                 nrow = n_pc+n_noise)
  if (decay_v) v_raw = v_raw * seq((n_pc+n_noise),1,-1)^2
  v = scale(v_raw, center = FALSE,
            scale = apply(v_raw, 2, function(x) norm(x, '2')))
  y = pcs %*% v + matrix(rnorm(dim_y*n, 0, y_sig), nrow = n)
  
  return(list(x = x, b = b, f_xb = f_xb, coords = cbind(coord_x, coord_y),
              pc_cov = pc_cov, pcs = pcs, v = v, y = y))
}

# run RapPCA with a given dataset and a specified r (# PCs)
prpca_one = function(x, y, coords, r, lambdas, lambda.ratios, gms, 
                     d = 1, m = 2, x.dim.reduction = FALSE, x.r = 5,
                     track.each.pc = TRUE, tst.ind = NULL, compact = FALSE){
  
  pars = expand.grid(gm = gms, lambda = lambdas, lambda.ratio = lambda.ratios)
  pars = rbind(#c(-1,0,0), 
               c(0,0,0), pars)
  
  coords = scale(coords)
  if (is.null(tst.ind))
    tst.ind = 1:floor(nrow(x)*.2)
  
  if (!is.null(x)){
    x.scaled = scale(x)
    if (x.dim.reduction){
      x.pca = prcomp(x.scaled)
      x.scaled = x.pca[, 1:x.r]
    }
    .xtrn = x.scaled[-tst.ind,]
    .xtst = x.scaled[tst.ind,] 
  } else {
    .xtrn = .xtst = NULL
  }
  
  y.scaled0 = y.scaled = y.scaled1 =  scale(y)
  rslt.pc1 = list()
  rslts = list()
  .rslt = .rslt.pca = list(v = rep(0, ncol(y)))
  loadings = uhats = utildes = NULL
  loadings.pca = uhats.pca = utildes.pca = NULL
  .r = 0
  # tmse = mspe = msre = 0
  
  while (r > 0){
    print(paste(r, 'PCs left for PCA'))
    .r = .r + 1
    r = r-1
    y.scaled = y.scaled - y.scaled %*% .rslt$v %*% t(.rslt$v)
    y.scaled1 = y.scaled1 - y.scaled1 %*% .rslt.pca$v %*% t(.rslt.pca$v)
    y.svd = svd(t(y.scaled[-tst.ind,])) 
    # dm = y.svd$d[1]
    dmat = diag(y.svd$d)
    U = y.svd$u
    V = y.svd$v
    y.svd1 = svd(t(y.scaled1[-tst.ind,])) 
    # dm = y.svd$d[1]
    dmat1 = diag(y.svd1$d)
    U1 = y.svd1$u
    V1 = y.svd1$v
    
    .rslts = mapply(rep_predPCA, gm = pars$gm, lambda = pars$lambda, 
                    lambda.ratio = pars$lambda.ratio,
                    MoreArgs = list(# xtrn = x.scaled[-tst.ind,], xtst = x.scaled[tst.ind,], 
                                    xtrn = .xtrn, xtst = .xtst, 
                                    ytrn = y.scaled[-tst.ind,], ytst = y.scaled[tst.ind,], 
                                    coords.trn = coords[-tst.ind,], coords.tst = coords[tst.ind,],
                                    d = d, m = m, dmat = dmat, U = U, V = V))
    .idx = which.min(unlist(.rslts[7,]))
    .rslt = .rslts[,.idx]
    ####.rslt.pca = .rslts[,2]
    #.rslt.pca = .rslts[,1]
    .rslt.pca = rep_predPCA(gm = -1, lambda = 0, lambda.ratio = 0,
                            # xtrn = x.scaled[-tst.ind,], xtst = x.scaled[tst.ind,], 
                            xtrn = .xtrn, xtst = .xtst, 
                            ytrn = y.scaled1[-tst.ind,], ytst = y.scaled1[tst.ind,], 
                            coords.trn = coords[-tst.ind,], coords.tst = coords[tst.ind,],
                            d = d, m = m, dmat = dmat1, U = U1, V = V1)
    
    if (track.each.pc){
      # rslt.pc1 = .rslts
      rslt.pc1[[.r]] = rbind(cbind(pars, unlist(.rslts[7,]), unlist(.rslts[8,]), 
                             unlist(.rslts[9,]), unlist(.rslts[10,]), 
                             unlist(.rslts[11,]), unlist(.rslts[12,]), unlist(.rslts[13,])),
                             c(-1,0,0, .rslt.pca$tmse, .rslt.pca$mspe, 
                               .rslt.pca$msre, .rslt.pca$msre.trn, 
                               .rslt.pca$r2, .rslt.pca$u.diff, .rslt.pca$condn))
      colnames(rslt.pc1[[.r]])[4:10] = c("tmse", "mspe", "msre", "msre_trn", 
                                        "r2", "mse_pc", "cond_n")
      # track.first.pc = FALSE
    } 
    
    rslts = append(rslts, list(.rslt))
    loadings = cbind(loadings, .rslt$v)
    uhats = cbind(uhats, .rslt$uhat.tst)
    utildes = cbind(utildes, .rslt$utilde.tst)
    
    loadings.pca = cbind(loadings.pca, .rslt.pca$v)
    uhats.pca = cbind(uhats.pca, .rslt.pca$uhat.tst)
    utildes.pca = cbind(utildes.pca, .rslt.pca$utilde.tst)
  }
  
  tmse = sum((y.scaled0[tst.ind,] - uhats%*%t(loadings))^2) / length(tst.ind) 
  mspe = sum(((uhats - utildes) %*% t(loadings))^2) / length(tst.ind)
  msre = sum((y.scaled0[tst.ind,] - utildes %*% t(loadings))^2) / length(tst.ind)
  u.diff = colMeans((utildes - uhats)^2)
  
  tmse.pca = sum((y.scaled0[tst.ind,] - uhats.pca%*%t(loadings.pca))^2) / length(tst.ind) 
  mspe.pca = sum(((uhats.pca - utildes.pca) %*% t(loadings.pca))^2) / length(tst.ind)
  msre.pca = sum((y.scaled0[tst.ind,] - utildes.pca %*% t(loadings.pca))^2) / length(tst.ind)
  u.diff.pca = colMeans((utildes.pca - uhats.pca)^2)
  
  if (compact){
    loadings = loadings.pca = NULL
  }
  
  return(list(tmse = tmse, mspe = mspe, msre = msre, u.diff = u.diff,
              tmse.pca = tmse.pca, mspe.pca = mspe.pca, 
              msre.pca = msre.pca, u.diff.pca = u.diff.pca,
              rslt.pc1 = rslt.pc1, 
              uhats = uhats, uhats.pca = uhats.pca,
              loadings = loadings, loadings.pca = loadings.pca,
              utildes = utildes, utildes.pca = utildes.pca))
}

# run Predictive PCA with a given dataset and a specified r (# PCs)
predpca_one = function(x, y, coords, r, #lambda2s, nvars,
                       m = 2, x.dim.reduction = TRUE, x.r = 10,
                       tst.ind = NULL, b.init = NULL, compact = FALSE){
  
  cov.given = !is.null(x)
  if (is.null(tst.ind)) tst.ind = 1:floor(nrow(x)*.2)
  if (cov.given & x.dim.reduction){
    Btrn = smooth.construct(s(coord1, coord2, bs="tp", k=100, m=m),
                            data = list(coord1 = coords[,1], coord2 = coords[,2]),
                            knots = NULL)
    xx = Btrn$X[,-which(apply(Btrn$X, 2, sd) < .01)]
    x.scaled = scale(cbind(x, xx))
    x.pca = prcomp(x.scaled)
    x.scaled = cbind(1,x.pca$x[, 1:x.r])
  } else if (cov.given){
    if (ncol(x) > 20){
      x.pca = prcomp(x, center = TRUE, scale. = TRUE)
      .x = x.pca$x[,1:x.r]
    } else{
      .x = scale(x)
    }
    Btrn = smooth.construct(s(coord1, coord2, bs="tp", k=4, m=m),
                            data = list(coord1 = coords[,1], coord2 = coords[,2]),
                            knots = NULL)
    xx = Btrn$X[,-which(apply(Btrn$X, 2, sd) < .01)]
    # x.scaled = scale(cbind(.x, xx))
    x.scaled = cbind(1, .x, scale(xx))
  } else {
    Btrn = smooth.construct(s(coord1, coord2, bs="tp", k=10, m=m),
                            data = list(coord1 = coords[,1], coord2 = coords[,2]),
                            knots = NULL)
    xx = Btrn$X[,-which(apply(Btrn$X, 2, sd) < .01)]
    x.scaled = cbind(1, scale(xx))
  }
  
  if (is.null(b.init)) b.init = runif(ncol(x.scaled), -1, 1)
  # px = ifelse(x.dim.reduction, x.r, ncol(x))
  coords = scale(coords)
  
  y.scaled0 = y.scaled = scale(y)
  x1train = coords[-tst.ind,1]
  x2train = coords[-tst.ind,2]
  # rslts = list()
  .v = rep(0, ncol(y))
  # .u = rep(0, nrow(y))
  loadings = uhats = utildes = msre.trn = NULL
  tmse1 = mspe1 = msre1 = r2 = u.diff1 = NULL # metrics for each PC, individually
  # tmse = mspe = msre = 0
  
  while (r > 0){
    print(paste(r, 'PCs left for PredPCA'))
    r = r-1
    y.scaled = y.scaled - y.scaled %*% .v %*% t(.v)
    # y.scaled = y.scaled - .u %*% t(.v)
    rslts = .uhats = list()
    
    .rslt = space.sparse.pca(Z = x.scaled[-tst.ind,], X = y.scaled[-tst.ind,],
                             beta.init = b.init,
                             v.init = runif(ncol(y.scaled), -1, 1),
                             lambda2 = 1, varnum1 = 1,
                             niter = 500, err = 1e-4)
    
    .v = .rslt$v/norm(.rslt$v, '2')
    # .v = .rslt$v
    # .u = .rslt$u
    # .beta = .rslt$beta
    # .u = x.scaled %*% .beta
    # .u = .u / norm(.u, '2')
    .uhat.trn = y.scaled[-tst.ind,]%*%.v
    # .uhat.trn = .rslt$u
    
    if (cov.given){
      rf <- randomForest(scale(x[-tst.ind,]), .uhat.trn,
                         xtest = scale(x[tst.ind,]), nodesize=5, importance = FALSE)
      # rf <- randomForest(.x[-tst.ind,], .uhat.trn,
      #                    xtest = .x[tst.ind,], nodesize=5, importance = FALSE)
      mod <- mgcv::gam((.uhat.trn - rf$predicted)~s(x1train, x2train,
                                                    bs="tp", k=160, m=m, fx=FALSE))
      .uhat.tst <- rf$test$predicted + predict(mod,
                                               data.frame(x1train = coords[tst.ind,1],
                                                          x2train = coords[tst.ind,2]))
    } else {
      mod <- mgcv::gam(.uhat.trn~s(x1train, x2train, bs="tp", k=300, m=m, fx=FALSE))
      .uhat.tst <- predict(mod, data.frame(x1train = coords[tst.ind,1],
                                           x2train = coords[tst.ind,2]))
    }
    
    .utilde = y.scaled[tst.ind,]%*%.v
    # .utilde = x.scaled[tst.ind,] %*% .beta
    # .utilde = .utilde / norm(.utilde, '2')
    
    loadings = cbind(loadings, .v)
    utildes = cbind(utildes, .utilde)
    uhats = cbind(uhats, .uhat.tst)
    
    msre.trn = c(msre.trn, 
                 sum((y.scaled[-tst.ind,] - .uhat.trn %*% t(.v))^2) / length(.uhat.trn))
    tmse1 = c(tmse1,
              sum((y.scaled[tst.ind,] - .uhat.tst%*%t(.v))^2) / length(tst.ind))
    mspe1 = c(mspe1,
              sum(((c(.uhat.tst) - c(.utilde)) %*% t(.v))^2) / length(tst.ind))
    msre1 = c(msre1, 
              sum((y.scaled[tst.ind,] - .utilde %*% t(.v))^2) / length(tst.ind))
    r2 = c(r2, 1-mean((c(.uhat.tst) - c(.utilde))^2)/var(c(.utilde)))
    u.diff1 = c(u.diff1, mean((c(.uhat.tst) - c(.utilde))^2))
  }
  
  tmse = sum((y.scaled0[tst.ind,] - uhats%*%t(loadings))^2) / length(tst.ind) 
  mspe = sum(((uhats - utildes) %*% t(loadings))^2) / length(tst.ind)
  msre = sum((y.scaled0[tst.ind,] - utildes %*% t(loadings))^2) / length(tst.ind)
  u.diff = colMeans((utildes - uhats)^2)
  
  if (compact)
    loadings.pred = NULL
  
  return(list(tmse = tmse, mspe = mspe, msre = msre, u.diff = u.diff,
              metrics.pc = cbind(tmse1, mspe1, msre1, msre.trn, r2, u.diff1),
              loadings.pred = loadings, utildes.pred = utildes, uhats.pred = uhats))
}

solve_one = function(x, y, coords, r, gms, lambdas, lambda.ratios,
                     # lambda2s, nvars,
                     d = 1, m = 2, x.dim.reduction = FALSE, 
                     x.dim.reduction2 = FALSE, x.r = 10, 
                     tst.ind = NULL, b.init = NULL, compact = FALSE){
  
  print("Step 1 started")
  rslts_prpca = prpca_one(x, y, coords, r, 
                          lambdas = lambdas, lambda.ratios = lambda.ratios, gms = gms,
                          d = d, m = m, x.dim.reduction = x.dim.reduction,
                          x.r = x.r, track.each.pc = TRUE, tst.ind = tst.ind,
                          compact = compact)
  print("Step 1 finished, step 2 started")
  rslts_predpca = predpca_one(x, y, coords, r, 
                              # lambda2s = lambda2s, nvars = nvars,
                              m = m, x.dim.reduction = x.dim.reduction2, 
                              x.r = x.r, tst.ind = tst.ind, b.init = b.init,
                              compact = compact)
  print("one fold done")
  # df = cbind(pars, unlist(rslts_prpca$rslt.pc1[7,]), unlist(rslts_prpca$rslt.pc1[8,]),
  #            unlist(rslts_prpca$rslt.pc1[9,]))
  # colnames(df)[4:6] = c("tmse", "mspe", "msre")

  return(list(rslt_prpc1 = rslts_prpca$rslt.pc1, 
              rslt_pr = c(rslts_prpca$tmse, rslts_prpca$mspe, rslts_prpca$msre, 
                          rslts_prpca$u.diff),
              rslt_pca = c(rslts_prpca$tmse.pca, rslts_prpca$mspe.pca, rslts_prpca$msre.pca, 
                           rslts_prpca$u.diff.pca),
              rslt_pred = c(rslts_predpca$tmse, rslts_predpca$mspe, rslts_predpca$msre, 
                            rslts_predpca$u.diff),
              metric1_pred = rslts_predpca$metrics.pc,
              loadings = rslts_prpca$loadings, utildes = rslts_prpca$utildes, 
              uhats = rslts_prpca$uhats,
              loadings.pca = rslts_prpca$loadings.pca, utildes.pca = rslts_prpca$utildes.pca, 
              uhats.pca = rslts_prpca$uhats.pca,
              loadings.pred = rslts_predpca$loadings.pred, utildes.pred = rslts_predpca$utildes.pred,
              uhats.pred = rslts_predpca$uhats.pred))
}

one_sim = function(dim_x, n_pc, n_noise, dim_y, seed,
                   cov_range = 1, cov_sill = .2, cov_nugget = .2, 
                   y_sig = .5, n = 200, d = 1, # seed = 1, 
                   mean_function = NULL, linear = TRUE, r = 5,
                   krige = TRUE, decay_v = FALSE,
                   gms = c(.1, .5, 1), lambdas = c(.05, .1, .5), lambda.ratios = c(.5, 1, 2),
                   #lambda2s = c(.1, .5, 1), nvars = c(3, 5)
                   tst.ind = NULL, b.init = NULL
                   ){
  
  dat = sim_y(dim_x, n_pc, n_noise, dim_y, 
              cov_range, cov_sill, cov_nugget, 
              y_sig, n, seed, mean_function, linear,
              krige, decay_v)
  rslt = try(solve_one(x = dat$x, y = dat$y, coords = dat$coords, r = r, d = d,
                   gms = gms, lambdas = lambdas, lambda.ratios = lambda.ratios,
                   tst.ind = tst.ind, b.init = b.init
                   #lambda2s = lambda2s, nvars = nvars
                   ))
  if (class(rslt)[1] == 'try-error')
    return(list(data = dat, rslt_pr = NA, 
                rslt_pred = NA, rslt_pca = NA,
                metric1_pr = NA, metric1_pred = NA))
  
  return(list(data = dat, 
              rslt_pr = rslt$rslt_pr, 
              rslt_pred = rslt$rslt_pred,
              rslt_pca = rslt$rslt_pca,
              metric1_pr = rslt$rslt_prpc1,
              metric1_pred = rslt$metric1_pred))
}

# run cross-validation
cv_pca = function(x, y, coords, r, gms, lambdas, lambda.ratios,
                  d = 1, m = 2, x.dim.reduction = FALSE, 
                  x.dim.reduction2 = FALSE, x.r = 10,
                  idx = NULL, seed = 1, K = 10, b.init = NULL, 
                  compact = FALSE){
  set.seed(seed)
  
  n = nrow(y)
  if (is.null(idx)){
    idx = sample(1:n, n, replace = FALSE)
  }
  
  .idx = list()
  n.each = round(n/K)
  for (j in 1:(K-1)){
    .idx[[j]] = idx[(n.each*(j-1)+1):(n.each*j)]
  }
  .idx[[K]] = idx[(n.each*(K-1)+1):n]
  
  cv_metrics = mclapply(X = .idx, FUN = solve_one,
                        x = x, y = y, coords = coords, r = r, 
                        gms = gms, lambdas = lambdas, lambda.ratios = lambda.ratios,
                        d = d, m = m, x.dim.reduction = x.dim.reduction,
                        x.dim.reduction2 = x.dim.reduction2,
                        x.r = x.r, b.init = b.init, compact = compact,
                        mc.cores = min(detectCores() - 1, K+1))
  
  pars = cv_metrics[[1]]$rslt_prpc1[[1]][,1:3]
  print(pars)
  npars = nrow(pars)
  pr_tmse = matrix(0, nrow = npars, ncol = r)
  utildes = utildes.pca = utildes.pred = matrix(NA, nrow = n, ncol = r)
  uhats = uhats.pca = uhats.pred = matrix(NA, nrow = n, ncol = r)
  for (k in 1:K){
    utildes[.idx[[k]],] = cv_metrics[[k]]$utildes
    utildes.pca[.idx[[k]],] = cv_metrics[[k]]$utildes.pca
    utildes.pred[.idx[[k]],] = cv_metrics[[k]]$utildes.pred
    uhats[.idx[[k]],] = cv_metrics[[k]]$uhats
    uhats.pca[.idx[[k]],] = cv_metrics[[k]]$uhats.pca
    uhats.pred[.idx[[k]],] = cv_metrics[[k]]$uhats.pred
    for (j in 1:r){
      pr_tmse[,j] = pr_tmse[,j] + cv_metrics[[k]]$rslt_prpc1[[j]][,'tmse']
    } 
  }
  idx_best = apply(pr_tmse, 2, which.min)
  print(idx_best)
  temp = pars[idx_best[1],]
  if (r > 1){
    for (.r in 2:r){
      temp = rbind(temp, pars[idx_best[.r],])
    }
  }
  
  print(temp)
  par_best = cbind(temp, 1:r)
  # par_best = cbind(rbind(pars[idx_best[1],],
  #                        pars[idx_best[2],],
  #                        pars[idx_best[3],]),
  #                  1:r)
  colnames(par_best)[4] = "pc"
  
  metric_pr = cv_metrics[[1]]$rslt_pr
  metric_pca = cv_metrics[[1]]$rslt_pca
  metric_pred = cv_metrics[[1]]$rslt_pred
  for (k in 2:K){
    metric_pr = metric_pr + cv_metrics[[k]]$rslt_pr
    metric_pca = metric_pca + cv_metrics[[k]]$rslt_pca
    metric_pred = metric_pred + cv_metrics[[k]]$rslt_pred
  }
  metrics = rbind(metric_pr/K, metric_pca/K, metric_pred/K)
  colnames(metrics) = c("tmse", "mspe", "msre", paste("pc", 1:r, sep=""))
  rownames(metrics) = c("proposed", "pca", "pred_pca")
  
  return(list(cv_raw = cv_metrics, par_best = par_best,
              metrics = metrics, 
              utildes = utildes, utildes.pca = utildes.pca, utildes.pred = utildes.pred,
              uhats = uhats, uhats.pca = uhats.pca, uhats.pred = uhats.pred,
              cv_idx = .idx))
}
