# Predictive PCA (Jandarov 2017)
# Credit to the author: Roman Jandarov
#####################################################################################################

thresh <- function(z, type=c("soft", "hard", "SCAD"), delta, a=3.7){

  ### z: argument
  ### type: thresholding rule type
  ### delta: thresholding level
  ### a: default choice for SCAD penalty

  if(type=="soft"){
    return(sign(z)*(abs(z)>=delta)*(abs(z)-delta))
  }
  
  if(type=="hard"){
    return(z*(abs(z)>delta))
  }
  
  if(type=="SCAD"){
    return(sign(z)*(abs(z)>=delta)*(abs(z)-delta)*(abs(z)<=2*delta)+
             ((a-1)*z-sign(z)*a*delta)/(a-2)*(2*delta<abs(z))*(abs(z)<=a*delta)+
             z*(abs(z)>a*delta))
  }
  
}

#####################################################################################################

frob.norm <- function(X, Z, beta, v, lambda2){
	# beta = 1 : p
	# v = 1 : p
	# lambda2 = penalty term

	u = c(Z%*%beta)
	u = u/norm.vec(u)

	tmp = (X - u%*%t(v))
	sum(diag(tmp%*%t(tmp))) + 2*lambda2*sum(abs(beta))
}
#frob.norm(X, Z, beta, v, lambda2 = 0)

#####################################################################################################

#optim.beta <- function(X, Z, v, lambda2, beta.init){
#	tmp.fun <- function(beta){
#		frob.norm(X, Z, beta, v, lambda2)
#	}
#	optim(beta.init, tmp.fun)$par
#}
#beta1 = optim.beta(X, Z, v = rep(1,17), lambda2 = 0, beta.init = #beta2)

optim.beta <- function(X, Z, v, lambda2, beta.init){
	Y = X%*%v/c((t(v)%*%v))
	lm(Y~-1+Z)$coef
}
#beta2 = optim.beta(X, Z, v = rep(1,17), lambda2 = 0, beta.init = #rep(2,dim(Z)[2]))

#frob.norm(X, Z, beta1, v, lambda2 = 0)
#frob.norm(X, Z, beta2, v, lambda2 = 0)

#cor(beta1, beta2)

#####################################################################################################

norm.vec <- function(vec){
	sqrt(sum(vec^2))
}
#norm.vec(rep(1,3))

#####################################################################################################

space.sparse.pca <-function(X, Z, beta.init, v.init, lambda2, varnum1, 
                            type = "soft", niter=100, err=10^(-3), trace=FALSE){

	u.d <- v.d <- 1
	iter <- 0
	while(u.d>err | v.d>err){

		u = c(Z%*%beta.init)
		u = u/norm.vec(u)

		iter <- iter+1

		v1 <- t(X)%*%u

		# if (varnum1 < length(v1)){
		# 	lambda1 <- sort(abs(v1))[length(v1)-varnum1]
		# } else {
		# 	lambda1 <- min(abs(v1))/2
		# }
		# 
		# v1 <- thresh(v1, type, lambda1, 3.7)

		beta1 = optim.beta(X, Z, v1, lambda2, beta.init)

		u.d <- sqrt(sum((beta1-beta.init)^2))
		v.d <- sqrt(sum((v1-v.init)^2))

		if(iter > niter){
			print("Fail to converge! Increase the niter!")
			break
		}
		
		beta.init <- beta1
		v.init <- v1
	}

	if(trace){
		print(paste("iter:", iter, "u.d:", u.d, "v.d:", v.d))
	}

	return(list(beta=cbind(beta1), v = cbind(v1), u = cbind(u)))
}
