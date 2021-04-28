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

size <-   function(
C,	 		        #critical value
R, 				      #restriction matrix (q times k, rank q)
X, 				      #design matrix (n times k, rank k)
hcmethod,       #-1:4; -1 = classical F-test without df adjustment; 0 = HC0(R), 1 = HC1(R), etc.
restr.cov,      #Covariance Matrix estimator computed with restricted residuals (TRUE or FALSE)
Mp, 			      #Starting values drawn from uniform distribution on the simplex
M1,				      #Number of initial candidates in 1st step optimization 
M2,				      #Number of initial candidates in 2nd step optimization
N0 = NULL, 			#Replications in Monte-Carlo for starting value generation (needed only if q > 1)
N1 = NULL, 			#Replications in Monte-Carlo for 1st step optimization (needed only if q > 1)
N2 = NULL, 			#Replications in Monte-Carlo for 2nd step optimization (needed only if q > 1)
tol = 1e-08,    #tolerance parameter used in checking invertibility of VCESTIMATOR in test statistic
control.1 = 	  #controls for constrOptim function in 1st step optimization
	list("reltol" = 1e-02, "maxit" = dim(X)[1]*20),
control.2 = 	  #controls for constrOptim function in 2nd step optimization
	list("reltol" = 1e-03, "maxit" = dim(X)[1]*30),
cores = 1,		  #number of cores used in computations
lower = 0,		  #lower bound for variances
eps.close = .0001, #closeness to boundary variance in starting value search
lim = 30000,    #input for davies function
acc = 0.001,   #input for davies function
levelCl = 0,    #level in [0, 1) that should be used in the starting value search in case of large C
                #if levelCl = 0, then no additional starting values are generated taking care of potentially
                #large critical values
LBcheck = FALSE, #compare C to theoretical lower bounds, and return 1 if C is smaller than the bound
as.tol = 1e-08  #tolerance level in checking rank conditions for verifying Assumptions 1, 2, non-constancy,
                #and also for computing lower bounds for size-controlling cvs 
){

###########################################################################
#run input checks
###########################################################################

check.C(C)
check.X.R(X, R)
check.hcmethod(hcmethod)
check.restr.cov(restr.cov)

  if(dim(R)[1] > 1){
  check.N.M(N0, N1, N2, Mp, M1, M2)
  } else {
  check.M(Mp, M1, M2)
  }

check.tol(tol)
check.cores(cores)
check.lower(lower, dim(X)[1])
check.eps.close(eps.close)
check.eps.close.lower(eps.close, lower, dim(X)[1])
check.levelCl(levelCl)
check.LBcheck(LBcheck)
check.as.tol(as.tol)

###########################################################################
# Elementary quantities
###########################################################################

n <- dim(X)[1]                        #sample size
k <- dim(X)[2]                        #number of regressors
q <- dim(R)[1]                        #number of restrictions
qrX <- qr(X)                          #qr decomposition of X
rloc <- rep(0, length = q)

if(restr.cov){
 qrM0lin <- qr(M0lin(X, R))            #qr decomp of basis of M0lin = M0-mu0
 factor.tmp2.loc <- backsolve(qr.R(qrX), diag(k))
 factor.tmp2.loc <- tcrossprod(factor.tmp2.loc)
 RF.loc <- factor.tmp2.loc%*%t(R)%*%Bfactor.matrix2(qrX, R)
 } else {
 qrM0lin <- NULL
 RF.loc <- NULL
}

Bfac <- Bfactor.matrix(qrX, n, R)     #R(X'X)^{-1}X'
Bfac2 <- Bfactor.matrix2(qrX, R)      #(R(X'X)^{-1}R')^{-1}

###########################################################################
# Check Assumption 1 (if hcmethod >= 0 & restr.cov = FALSE) or
# Check Assumption 2 (if hcmethod >= 0 & restr.cov = TRUE)

# If the respective assumption does not hold STOP
###########################################################################

if(hcmethod >= 0){

 if(!restr.cov){
 if(!As.check(qrX, Bfac, q, as.tol)){
  stop("Assumption 1 seems to be violated. Perhaps change 'as.tol' value.")
 }
 } 
 
 if(restr.cov){
 if(!As.check(qrM0lin, Bfac, q, as.tol)){
  stop("Assumption 2 seems to be violated. Perhaps change 'as.tol' value.")
 }
 }
 
}

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
}

###########################################################################
# If hcmethod >=0 and restr.cov = TRUE, check that the test statistic
# is not constant on the complement of B tilde
###########################################################################

if( hcmethod >= 0 & restr.cov == TRUE){

if( q > 1 ){

  Yconst <- matrix(rnorm(n*100), nrow = n, ncol = 100)
  
  V <- F.wrap(Yconst, R, rloc, X, n, k, q, qrX, hcmethod, cores, 
 	restr.cov, Bfac, Bfac2, qrM0lin, RF.loc, tol)
 	
 	if(max(V)-min(V) < as.tol){
  stop("Test statistic seems to be constant on the complement of Btilde. 
  Perhaps decrease 'as.tol' value.")
  } 

} else {

 C0 <- sum(v^2)/sum((v*Wvecprod)^2)
 
 if(abs(C - C0) < as.tol){
 if(max(abs(MQF)) < as.tol){
  stop("Test statistic seems to be constant on the complement of Btilde. 
  Perhaps decrease 'as.tol' value.") 
  }
 }

}

}

###########################################################################
# Make use of theoretical lower bounds (that imply size = 1) if LBcheck = 
# TRUE
###########################################################################

if(LBcheck){

 LB <- LB.cv(X, R, hcmethod, restr.cov, as.tol, checks = FALSE)

 if(C < LB){
 stop("Theoretical results imply that the size equals 1, as C is smaller than
 lower bounds for size-controlling critical values. Set 'LBcheck = FALSE' to
 compute the size numerically nevertheless.")
 }

}

###########################################################################
#create storage space for results
###########################################################################

fs.results <- matrix(NA, nrow = M1, ncol = 3*(n+1)+1)

###########################################################################
#icheck if C is larger than 5 times the critical value (init.homo in the following)
#for which the rejection probability under homoskedasticity equals levelCl
###########################################################################

Clstart <- NULL

if(levelCl > 0){

init.homo <- explore.qm(levelCl, R, X, hcmethod, restr.cov, N2, 
tol, cores, lower, lim, acc, Bfac, Bfac2, qrM0lin, RF.loc, qrX, as.tol,
AScheck = FALSE)$quant.max

 if(C > 5*init.homo){
 
 #run size function (recursively) to obtain second stage paramters
 #which realize the size for the smaller critical value init.homo
 #these are used further below as an additional set of starting values
 
 #size.tol parameter is not used in this call
 
 Clstart <- c(size.qm(init.homo, 0, R, X, hcmethod, restr.cov, Mp, M1, 1, N0, N1, 
 N2, tol, control.1, control.2, cores, lower, eps.close, lim, acc, Bfac, Bfac2, 
 qrM0lin, RF.loc, qrX, size.tol=.0001, only.second.stage = TRUE)$second.stage.parameters)
 
 } 
}


###########################################################################
#generate starting values
###########################################################################

if(q == 1 & lower == 0){
 s0 <- -1
} else {
 s0 <- 0
}

if(!is.null(Clstart)){
 N <- n+1
} else {
 N <- n
}

for(i in s0:N){
 #variances maximizing the expectation of the quadratic form
 if(i < 0){
  tpara <- rep(eps.close, length = n)
  maxdiag <- (diag(MQF) == max(diag(MQF)))
  nbrmax <- sum(maxdiag)
  tpara[maxdiag] <- (1-(n-nbrmax)*eps.close)/nbrmax
 }
 #variances corresponding to homoskedasticity
 if(i == 0){
  tpara <- rep(n^(-1), length = n)
 } 
 #variances with a dominant coordinate
 if( i > 0 & (i <= n ) ){
  tpara <- rep(lower + eps.close, length = n)
  tpara[i] <- 1-(n-1)*(lower + eps.close)
 }
 #variances from a large C search 
 if(i == n+1){
 tpara <- Clstart
 }
 
 if(q > 1){
 y <- rnorm(N0 * n)
 }
 val <- objective.fun(tpara[1:(n-1)])
 M <- rbind(fs.results[, 1:(n+1)], c(tpara, val))
 M <- M[order(M[,n+1]),]
 fs.results[, 1:(n+1)] <- M[1:M1,]
}

#variances drawn randomly from a uniform distribution on the standard simplex
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

#return outputs

return(list(
"starting.parameters" = fs.results[,1:n],
"starting.rejection.probs" = - fs.results[,n+1],
"first.stage.parameters" = fs.results[, (n+2):(2*n +1)],
"first.stage.rejection.probs" = - fs.results[,2*n+2],
"second.stage.parameters" = fs.results[index.top, (2*n+3):(3*n +2)], 
"second.stage.rejection.probs" = - fs.results[index.top,3*n+3],
"convergence" = fs.results[index.top,3*n+4],
"size" = max(- fs.results[index.top, 3*n+3])
		)
	)
}