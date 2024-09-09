library(tidyverse)
library(SpatialPCA)
library(nnet) # multinom function
library(DescTools) # PseudoR2 function
library(mclust) # ARI
library(cluster)
library(pbmcapply)

source("rpPCA_profile.R")
load("Tumor_data.RData")
load("tumor_data_norm.RData")
multinom = nnet::multinom

# annotated data as ground truth for domain detection
# H1 sample was used in Shang(2022)
anno = read.csv("H1_labeled_coordinates.tsv" ,sep="\t") %>%
  mutate(rown = paste(round(x), "x", round(y), "_1", sep = "")) %>%
  filter(rown %in% rownames(location))
loc_anno = location %>% as.data.frame %>% rownames_to_column("rown") %>%
  left_join(anno, by = "rown")

rap_clust = function(y, coords, gm, lambda, r, true_lbl, 
                     n_clust = 7, knn = NULL, d = 1, m = 2, na_ind = NULL){
  p = ncol(y)
  n = nrow(y)
  x1train = coords[,1]
  x2train = coords[,2]
  
  y.scaled0 = y.scaled = scale(y)
  loadings = uhats = utildes = y.res = NULL
  .rslt = list(v = rep(0, ncol(y)))
  .r = 0
  
  while (r > 0){
    .r = .r + 1
    r = r-1
    y.scaled = y.scaled - y.scaled %*% .rslt$v %*% t(.rslt$v)
    y.svd = try(svd(t(y.scaled)))
    if (class(y.svd)[1] == 'try-error'){
      return(list(ARI = NA, pseudoR2 = NA, silhouette_score = NA,
                  gm = gm, lambda = lambda,
                  loadings = loadings, utildes = utildes))
    }
    dmat = diag(y.svd$d)
    U = y.svd$u
    V = y.svd$v
    
    .rslt = rep_predPCA(gm = gm, lambda = lambda, lambda.ratio = 1,
                        xtrn = NULL, xtst = NULL, 
                        ytrn = y.scaled, ytst = y.scaled, 
                        coords.trn = coords, coords.tst = coords,
                        d = d, m = m, dmat = dmat, U = U, V = V)
    
    y.res = cbind(y.res, y.scaled)
    loadings = cbind(loadings, .rslt$v)
    uhats = cbind(uhats, .rslt$uhat.tst)
    utildes = cbind(utildes, .rslt$utilde.tst)
  }
  
  if (is.null(knn)) knn = round(sqrt(n))
  clust_lbl = walktrap_clustering(n_clust, t(utildes), knn)
  if (!is.null(na_ind)){
    true_lbl1 = true_lbl[-na_ind]
    clust_lbl1 = clust_lbl[-na_ind]
    utildes1 = utildes[-na_ind,]
  } else {
    true_lbl1 = true_lbl
    clust_lbl1 = clust_lbl
    utildes1 = utildes
  }
  
  ari = adjustedRandIndex(true_lbl1, clust_lbl1)
  fit = multinom(true_lbl ~ utildes, maxit = 1000, MaxNWts = 2000, model = TRUE)
  pseudoR2 = PseudoR2(fit, c("McFaddenAdj"))
  
  temp = silhouette(x = as.integer(clust_lbl1), dmatrix = dist(utildes1))
  silhouette_score = summary(temp)$avg.width
  
  return(list(ARI = ari, pseudoR2 = pseudoR2, silhouette_score = silhouette_score,
              gm = gm, lambda = lambda,
              loadings = loadings, utildes = utildes))
}

gms = seq(.5, 5, .5)
lambdas = 10^(seq(-1, .8, .2))

pars = expand.grid(gms, lambdas)
clust_tuning = pbmcmapply(rap_clust, gm = pars[,1], lambda = pars[,2], #r = pars[,3],
                        MoreArgs = list(y = y, coords = coords, true_lbl = loc_anno$label,
                                        r = 20, n_clust = 7, knn = NULL, 
                                        na_ind = which(loc_anno$label == "undetermined")),
                        mc.cores = 5)
save(pars, clust_tuning, file = "cl_tuning_demo_finer.RData")