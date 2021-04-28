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

# explore.qm is only used locally in the functions: critical.value, size
# therefore no default inputs are given if not necessary

# it computes the 1-alpha quantile of the distribution of the test
# statistic under homoskedasticity

explore.qm <-   function(
alpha,          #targeted level of significance
R, 				      #restriction matrix (q times k, rank q)
X, 				      #design matrix (n times k, rank k)
hcmethod,     	#-1:4; -1 = classical F-test without df adjustment; 0 = HC0, 1 = HC1, etc.
restr.cov,  	  #Covariance Matrix estimator computed with restricted residuals (TRUE or FALSE)
N2, 		      	#Replications in 2nd step optimization
tol,            #tolerance parameter used in checking invertibility of VCESTIMATOR in test statistic
cores,    		  #number of cores used in computations
lower,		      #lower bound for variances
lim,            #input for davies function
acc,            #input for davies function
Bfac,		        #R(X'X)^{-1}X'
Bfac2, 		      #(R(X'X)^{-1}R')^{-1}
qrM0lin,	    	#qr decomp of M0lin (only if restr.cov = TRUE)
RF.loc,		      #auxiliary function for computing restricted estimators (only if restr.cov = TRUE)
qrX, 			      #qr decomp of X
as.tol = 1e-08,  #tolerance level in checking rank conditions for verifying Assumptions 1, 2, non-constancy,
                #and also for computing lower bounds for size-controlling cvs 
AScheck = FALSE #should Assumptions be checked (Assumption 1 or 2, constancy, existence of size-controlling cv
){


###########################################################################
# Elementary quantities
###########################################################################

n <- dim(X)[1]                        #sample size
k <- dim(X)[2]                        #number of regressors
q <- dim(R)[1]                      #number of restrictions
rloc <- rep(0, length = q)

###########################################################################
# Objective function (depends also on ``sample'' y if q > 1) - minimize 
# negative rejection probability
###########################################################################

if(q == 1){

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
  TCP <- tcrossprod(Projmat%*% diag(Wvecprod))  
  } else {
  fact <- backsolve(qr.R(qrX), diag(k)) 
  fact <- sum((R%*%fact)^2)
  } 
  
  objective.fun2 <- function(C2){
  sqvariances <- sqrt(lower + (1-n*lower)*c(para.max, 1-sum(para.max)))
  
  if(hcmethod != -1){
   MQF2 <- tcrossprod(v) - C2* TCP
  } else {
   MQF2 <- tcrossprod(v) - (C2*fact/(n-l))*Projmat
  } 

  fullmiddle <- MQF2*tcrossprod(sqvariances)
  eivec <- eigen(fullmiddle, only.values = TRUE, symmetric = TRUE)$values
  return(davies(0, lambda = eivec, lim = lim, acc = acc)$Qq - alpha)
  }
}

###########################################################################
#CHECK ASSUMPTIONS 1 or 2 and non-constancy condition (if necessary).
###########################################################################

if(AScheck){

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
  #check sufficient condition for existence of size controlling critical value
  ###########################################################################

  sc <- SC.check(X, R, hcmethod, restr.cov, as.tol, checks = FALSE)

  if(!sc){
  warning("The output of SC.check suggests that the sufficient conditions for
  the existence of a size-controlling critical value (obtained in the paper 
  this package relies on) do not seem to hold! This can be ignored in case
  'lower' is nonzero.")
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
   MQFT <- tcrossprod(v) - C0* TCP
   if(max(abs(MQFT)) < as.tol){
    stop("Test statistic seems to be constant on the complement of Btilde. 
    Perhaps decrease 'as.tol' value.") 
    }
   }
  }
}

#homoskedastic configuration

para.max <- rep(1/n, length = n-1)

#compute quantile based on homoskedastic configuration

if(q > 1){

  y <- rnorm(N2 * n)
  variances <- lower + (1-n*lower)*c(para.max, 1-sum(para.max))
  sqvariances <- rep(sqrt(variances), length = length(y))
  ytransf <- matrix(sqvariances*y, nrow = n, byrow = FALSE)
  Rbmat <- R %*% qr.coef(qrX, ytransf)
  vals <- F.wrap(ytransf, R, rloc, X, n, k, q, qrX, hcmethod, cores, 
  restr.cov, Bfac, Bfac2, qrM0lin, RF.loc, tol)
  quant.max <- max(quant.max, quantile(vals, probs = 1-alpha))
 
} else {

  mqcand <- uniroot(objective.fun2, interval = c(0, 3*qchisq(1-alpha, df = 1)), extendInt = "downX", tol = alpha/100)
  quant.max <- mqcand$root

}

#return quantile

return(list("quant.max" = quant.max))
}