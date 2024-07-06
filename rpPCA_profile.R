# optimization by setting u = Yv and profiling out beta

library(mgcv)
library(Matrix)
library(randomForest)
library(parallel)
# library(RSpectra)
library(pracma)
## library(reticulate)

## np = import("numpy")

polyk = function(x, y, d){
  k = (1 + x %*% y)^d
  return(k)
}
polyK = function(X, Y, d){
  return(outer(X = 1:nrow(X), Y = 1:nrow(Y), 
               FUN = Vectorize(function(i, j) polyk(X[i,], Y[j,], d)))
  )
}

rep_predPCA = function(xtrn, xtst, ytrn, ytst, coords.trn, coords.tst,
                       lambda, lambda.ratio = 1, gm, d = 1, m = 2,
                       dmat = parent.frame()$dmat, U = parent.frame()$U,
                       V = parent.frame()$V, return.condn = FALSE,
                       delta = 0.05){
  
  # print(paste('calculation started, gm =', round(gm,2), 'lambda =', round(lambda,2)))
  # this function extracts one PC only
  n.trn = nrow(ytrn)
  n.tst = nrow(ytst)
  cov.given = !is.null(xtrn)
  
  Btrn = smooth.construct(s(coord1, coord2, bs="tp", k=n.trn, m=m),
                          data = list(coord1 = coords.trn[,1], coord2 = coords.trn[,2]),
                          knots = NULL)
  S = Btrn$S[[1]]
  if (cov.given){
    Xtrn = polyK(xtrn, xtrn, d)
    Ztrn = cbind(1, Xtrn, Btrn$X/lambda.ratio)
    pen = as.matrix(bdiag(1, Xtrn, S/lambda.ratio))
    # Ztrn = cbind(Xtrn, Btrn$X/lambda.ratio)
    # pen = as.matrix(bdiag(Xtrn, S/lambda.ratio))
  } else {
    Ztrn = cbind(1, Btrn$X)
    pen = bdiag(1, S)
  }
  # Xtst = polyK(xtst, xtrn, d) # X capitalized: kernel; x not capitalized: raw
  p = ncol(Ztrn)
  D = dmat
  
  if (gm > 0){
    # A = (gm-1) * D^2 - gm^2*D%*%t(V)%*%Ztrn%*%solve(gm*t(Ztrn)%*%Ztrn+diag(lambda, p))%*%t(Ztrn)%*%V%*%D
    A = (gm-1) * D^2 - gm^2*D%*%t(V)%*%Ztrn%*%solve(gm*t(Ztrn)%*%Ztrn+lambda*pen+diag(delta, p))%*%t(Ztrn)%*%V%*%D
    # eig = eigs_sym(-A, 1)
    # q = eig$vectors
    ## eig = np$linalg$eig(-A)
    ## m_idx = which.max(eig[[1]])
    ## q = eig[[2]][,m_idx]
    eig = eigen(-A, symmetric = TRUE)
    m_idx = which.max(eig$values)
    q = eig$vectors[,m_idx]
    ab = solve(gm*t(Ztrn)%*%Ztrn + lambda*pen + diag(delta, p), gm*t(Ztrn)%*%V%*%D%*%q)
    # ab = solve(gm*t(Ztrn)%*%Ztrn + lambda*diag(p), gm*t(Ztrn)%*%V%*%D%*%q)
    v = U%*%q
    condn = ifelse(return.condn, cond(A, 2), NA)
  } else if (gm == 0) {
    v = U[,1] 
    q = t(U) %*% v
    ab = rep(0, p)
    condn = 0
  } else {
    .pca = prcomp(ytrn, center = FALSE, scale. = FALSE)
    v = .pca$rotation[,1]
    q = t(U) %*% v
    ab = rep(0, p)
    condn = 0
  }
  
  u.orig = ytrn%*%v
  # obj.orig = sum((ytrn - u.orig%*%t(v))^2) + gm*sum((u.orig - Ztrn%*%ab)^2) + lambda*sum(ab^2)
  obj.orig = sum((ytrn - u.orig%*%t(v))^2) + gm*sum((u.orig - Ztrn%*%ab)^2) + delta*sum(ab^2) + lambda*t(ab)%*%pen%*%ab
  
  # evaluate on test data
  uhat.trn = u.orig # ytrn %*% v
  
  x1train = coords.trn[,1]
  x2train = coords.trn[,2]
  if (cov.given){
    rf <- randomForest(xtrn, uhat.trn, xtest = xtst, nodesize=5, importance = FALSE)
    mod <- mgcv::gam((uhat.trn - rf$predicted)~s(x1train, x2train, 
                                                 bs="tp", k=n.trn, m=m, fx=FALSE))
    uhat.tst <- rf$test$predicted + predict(mod, data.frame(x1train = coords.tst[,1], x2train = coords.tst[,2]))
  } else {
    mod <- mgcv::gam(uhat.trn~s(x1train, x2train, bs="tp", k=n.trn, m=m, fx=FALSE))
    uhat.tst = predict(mod, data.frame(x1train = coords.tst[,1], x2train = coords.tst[,2]))
  }

  uhat.tst = matrix(uhat.tst, ncol = 1)
  
  tmse = sum((ytst - uhat.tst%*%t(v))^2) / n.tst # overall mse
  utilde.tst = ytst %*% v
  
  mspe = sum(((uhat.tst - utilde.tst) %*% t(v))^2) / n.tst # mean squared prediction error
  msre = sum((ytst - utilde.tst %*% t(v))^2) / n.tst # mean squared representation error
  
  msre.trn = sum((ytrn - u.orig %*% t(v))^2) / n.trn
  r2 = 1 - mean((uhat.tst - utilde.tst)^2) / var(utilde.tst)
  u.diff = mean((uhat.tst - utilde.tst)^2)
  
  return(list(q = q, v = v, ab = ab, 
              uhat.tst = uhat.tst, utilde.tst = utilde.tst, 
              obj.orig = obj.orig,
              tmse = tmse, mspe = mspe, msre = msre,
              msre.trn = msre.trn, r2 = r2, u.diff = u.diff,
              condn = condn))
}
