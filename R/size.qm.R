#
#    Copyright (C) 2021  David Preinerstorfer
#    david.preinerstorfer@ulb.ac.be
#
#    This file is a part of hrt.
#
#    hrt is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/

# size.qm is only used locally in the function: critical.value, size
# hence no default inputs are provided

size.qm <-   function(
C,	 		      #critical value
alpha,        #targeted level of significance
R, 				    #restriction matrix (q times k, rank q)
X, 				    #design matrix (n times k, rank k)
hcmethod,   	#-1:4; -1 = classical F-test without df adjustment; 0 = HC0, 1 = HC1, etc.
restr.cov,  	#Covariance Matrix estimator computed with restricted residuals (TRUE or FALSE)
Mp, 			    #Starting values
M1,				    #Number of initial candidates in 1st step optimization 
M2,				    #Number of initial candidates in 2nd step optimization
N0, 		      #Replications in starting value generation (needed only if q > 1)
N1, 		      #Replications in 1st step optimization (needed only if q > 1)
N2, 		      #Replications in 2nd step optimization (needed only if q > 1)
tol,          #tolerance parameter used in checking invertibility of VCESTIMATOR in test statistic
control.1,  	#controls for constrOptim function in 1st step optimization
control.2,   	#controls for constrOptim function in 2nd step optimization
cores,		    #number of cores used in computations
lower,		    #lower bound for variances
eps.close,    #closeness to boundary variance in starting value search
lim,          #input for davies function
acc,          #input for davies function
Bfac,		      #R(X'X)^{-1}X'
Bfac2, 		    #(R(X'X)^{-1}R')^{-1}
qrM0lin,	  	#qr decomp of M0lin (only if restr.cov = TRUE)
RF.loc,		    #auxiliary function for computing restricted estimators (only if restr.cov = TRUE)
qrX,  		    #qr decomp of X
size.tol,     #tolerance level in final size
only.second.stage = FALSE
){

###########################################################################
#run input checks 
###########################################################################

check.C(C)

###########################################################################
# Elementary quantities
###########################################################################

n <- dim(X)[1]                        #sample size
k <- dim(X)[2]                        #number of regressors
q <- dim(R)[1]                        #number of restrictions
rloc <- rep(0, length = q)

###########################################################################
# Objective function (depends also on ``sample'' y if q > 1) - minimize 
# negative rejection probability
###########################################################################

if(q > 1){

 objective.fun <- function(z){
  variances <- lower + (1-n*lower)*c(z, 1-sum(z))
  sqvariances <- rep(sqrt(variances), length = length(y))
  ytransf <- matrix(sqvariances*y, nrow = n, byrow = FALSE)
  Rbmat <- R %*% qr.coef(qrX, ytransf)

  #evaluate test statistic at transformed values

  vals <- F.wrap(ytransf, R, rloc, X, n, k, q, qrX, hcmethod, cores, 
 	restr.cov, Bfac, Bfac2, qrM0lin, RF.loc, tol)
  return(- mean(vals >= C))
  
 }
 
} else {

  if(restr.cov){
  Projmat <- qr.Q(qrM0lin)
  l <- (k-q)
  qrXQ <- qrM0lin
  } else {
  Projmat <- qr.Q(qrX)
  l <- k
  qrXQ <- qrX
  }

  Projmat <- diag(n) - tcrossprod(Projmat)

  v <- c(Bfac)

  if(hcmethod != -1){
  Wvecprod <- sqrt(wvec(l, n, qrXQ, hcmethod))*v    
  MQF <- tcrossprod(v) - C* tcrossprod(Projmat%*% diag(Wvecprod))
  TCP <- tcrossprod(Projmat%*% diag(Wvecprod))
  } else {
  fact <- backsolve(qr.R(qrX), diag(k)) 
  fact <- sum((R%*%fact)^2)
  MQF <- tcrossprod(v) - (C*fact/(n-l))*Projmat
  } 

 objective.fun <- function(z){
  sqvariances <- sqrt(lower + (1-n*lower)*c(z, 1-sum(z)))
  fullmiddle <- MQF*tcrossprod(sqvariances)
  eivec <- eigen(fullmiddle, only.values = TRUE, symmetric = TRUE)$values
  return(-davies(0, lambda = eivec, lim = lim, acc = acc)$Qq)
 }
  
 objective.fun2 <- function(C2){
  sqvariances <- sqrt(lower + (1-n*lower)*c(para.max, 1-sum(para.max)))
  
  if(hcmethod != -1){
   MQF2 <- tcrossprod(v) - C2*TCP
  } else {
   MQF2 <- tcrossprod(v) - (C2*fact/(n-l))*Projmat
  } 

  fullmiddle <- MQF2*tcrossprod(sqvariances)
  eivec <- eigen(fullmiddle, only.values = TRUE, symmetric = TRUE)$values
  return(davies(0, lambda = eivec, lim = lim, acc = acc)$Qq - alpha)
  }
}


###########################################################################
#create storage space for results
###########################################################################

fs.results <- matrix(NA, nrow = M1, ncol = 3*(n+1)+1)

###########################################################################
#generate starting values
###########################################################################

if(q == 1 & lower == 0){
 s0 <- -1
} else {
 s0 <- 0
}

for(i in s0:n){
 if(i < 0){
  tpara <- rep(eps.close, length = n)
  maxdiag <- (diag(MQF) == max(diag(MQF)))
  nbrmax <- sum(maxdiag)
  tpara[maxdiag] <- (1-(n-nbrmax)*eps.close)/nbrmax
 }
 if(i == 0){
  tpara <- rep(n^(-1), length = n)
 } 
 if(i > 0){
  tpara <- rep(lower + eps.close, length = n)
  tpara[i] <- 1-(n-1)*(lower + eps.close)
 }
 if(q > 1){
 y <- rnorm(N0 * n)
 }
 val <- objective.fun(tpara[1:(n-1)])
 M <- rbind(fs.results[, 1:(n+1)], c(tpara, val))
 M <- M[order(M[,n+1]),]
 fs.results[, 1:(n+1)] <- M[1:M1,]
}

for(i in 1:Mp){
 tpara <- rexp(n, rate = 1)
 tpara <- tpara/(sum(tpara))
 if(i >= Mp/4){
 tpara <- lower + (1-n*lower)*tpara^2/sum(tpara^2)
 } else {
 tpara <- lower + (1-n*lower)*tpara
 }
 if(q > 1){
 y <- rnorm(N0 * n)
 }
 val <- objective.fun(tpara[1:(n-1)])
 M <- rbind(fs.results[, 1:(n+1)], c(tpara, val))
 M <- M[order(M[,n+1]),]
 fs.results[, 1:(n+1)] <- M[1:M1,]
}

###########################################################################
#initiate first stage optimizations
###########################################################################

#restriction parameters for constrained optimization

ui <- rbind(diag(n-1), rep(-1, length = n-1))
ci <- c(rep(lower, length = n-1), -(1-lower))

for(i in 1:M1){
 if(q > 1){
 y <- rnorm(N1 * n)
 }
 startval <- fs.results[i,1:(n-1)]
 max.rjp <- constrOptim(startval, objective.fun, NULL, ui, ci, control = control.1)
 fs.results[i,(n+2):(2*n +1)] <-  c(max.rjp$par, 1-sum(max.rjp$par))
 fs.results[i,2*n+2] <- max.rjp$value
}

###########################################################################
#initiate second stage optimizations
###########################################################################

index.top <- order(fs.results[,2*n+2])[1:M2]

for(i in index.top){
 if(q > 1){
 y <- rnorm(N2 * n)
 }
 startval <- fs.results[i,(n+2):(2*n)]
 max.rjp <- constrOptim(startval, objective.fun, NULL, ui, ci, control = control.2)
 fs.results[i,(2*n+3):(3*n +2)] <-  c(max.rjp$par, 1-sum(max.rjp$par))
 fs.results[i,3*n+3] <- max.rjp$value
 fs.results[i,3*n+4] <- max.rjp$convergence
}

if(only.second.stage){
M <- fs.results[index.top, (2*n+3):(3*n +2), drop = FALSE]
return(list("second.stage.parameters" = M))
}

quant.max <- NA

if(max(- fs.results[index.top, 3*n+3]) > alpha){

#obtain worst case parameter configurations

M <- fs.results[index.top, (2*n+3):(3*n +1), drop = FALSE]

#compute quantile based on worst-case parameter configuration

quant.max <- 0

if(q > 1){

 for(i in 1:M2){
  para.max <- M[i,]
  variances <- lower + (1-n*lower)*c(para.max, 1-sum(para.max))
  sqvariances <- rep(sqrt(variances), length = length(y))
  ytransf <- matrix(sqvariances*y, nrow = n, byrow = FALSE)
  Rbmat <- R %*% qr.coef(qrX, ytransf)
  vals <- F.wrap(ytransf, R, rloc, X, n, k, q, qrX, hcmethod, cores, 
  restr.cov, Bfac, Bfac2, qrM0lin, RF.loc, tol)
  quant.max <- max(quant.max, quantile(vals, probs = 1-alpha))
 }
 
} else {

 for(i in 1:M2){
  para.max <- M[i,]
  mqcand <- uniroot(objective.fun2, interval = c(0, 3*qchisq(1-alpha, df = 1)), extendInt = "downX", tol = size.tol)
  quant.max <- max(quant.max, mqcand$root)
 }
 
}

}

#return outputs

return(list(
"quant.max" = quant.max,
"size" = max(- fs.results[index.top, 3*n+3])
		)
	)
}