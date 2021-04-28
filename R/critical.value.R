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

critical.value <- function(
alpha,	 	    #significance level
R, 				    #restriction matrix (q times k, rank q)
X, 				    #design matrix (n times k, rank k)
hcmethod,   	#-1:4; -1 = classical F-test without df adjustment; 0 = HC0, 1 = HC1, etc.
restr.cov,  	#Covariance Matrix estimator computed with restricted residuals (TRUE or FALSE)
Mp, 			    #Starting values drawn from uniform distribution on the simplex
M1,				    #Number of initial candidates in 1st step optimization 
M2,				    #Number of initial candidates in 2nd step optimization
N0 = NULL, 		#Replications in Monte-Carlo for starting value generation (needed only if q > 1)
N1 = NULL, 		#Replications in Monte-Carlo for 1st step optimization (needed only if q > 1)
N2 = NULL, 		#Replications in Monte-Carlo for 2nd step optimization (needed only if q > 1)
tol = 1e-08,  #tolerance parameter used in checking invertibility of VCESTIMATOR in test statistic
control.1 = 	#controls for constrOptim function in 1st step optimization
	list("reltol" = 1e-02, "maxit" = dim(X)[1]*20),
control.2 = 	#controls for constrOptim function in 2nd step optimization
	list("reltol" = 1e-03, "maxit" = dim(X)[1]*30),
cores = 1,		#number of cores used in computations
lower = 0,		#lower bound for variances
eps.close = .0001, #closeness to boundary variance in starting value search in size computations
lim = 30000,  #input for davies function
acc = 0.001, #input for davies function
size.tol = .001, #tolerance level in final size
maxit = 25,    #maximum number of iterations
as.tol = 1e-08  #tolerance level in checking rank conditions for verifying Assumptions 1, 2, non-constancy,
                #and also for computing lower bounds for size-controlling cvs 
){

###########################################################################
#run input checks
###########################################################################

check.alpha(alpha)
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
check.hcmethod(hcmethod)
check.restr.cov(restr.cov)
check.as.tol(as.tol)
check.size.tol(size.tol)
check.maxit(maxit)

###########################################################################
# Elementary quantities
###########################################################################

n <- dim(X)[1]                        #sample size
k <- dim(X)[2]                        #number of regressors
q <- dim(R)[1]                        #number of restrictions
qrX <- qr(X)                          #qr decomposition of X

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
#compute lower bound for size controlling critical value
###########################################################################

if(lower == 0){
lbcv <- LB.cv(X, R, hcmethod, restr.cov, as.tol, checks = FALSE)
} else {
lbcv <- 0
}

###########################################################################
#auxiliary function that returns the size for a given critical value C
###########################################################################

fuC <- function(C){

size.qm(C, alpha, R, X, hcmethod, restr.cov, Mp, M1, M2, N0, N1, N2, 
tol, control.1, control.2, cores, lower, eps.close, lim, acc, Bfac, Bfac2,
qrM0lin, RF.loc, qrX, size.tol, only.second.stage = FALSE)

}

###########################################################################
#size of initial cv (maximum of 1-alpha quantile of distribution of test
#stat under homoskedasticity, and the lower bounds)

#in the course of explore.qm also the Assumptions are checked via AScheck:
#ASSUMPTIONS 1 or 2 and non-constancy condition (if necessary).
###########################################################################

init.homo <- explore.qm(alpha, R, X, hcmethod, restr.cov, N2, 
tol, cores, lower, lim, acc, Bfac, Bfac2, qrM0lin, RF.loc, qrX, as.tol, 
AScheck = TRUE)

cvrun <- max(init.homo$quant.max, lbcv)

sizerun <- fuC(cvrun)
it <- 0

###########################################################################
# iteration
###########################################################################

while((sizerun$size > alpha + size.tol) & (it < maxit)){

 show(paste("Iteration Number:", it))

 it <- it + 1
 cvrun <- sizerun$quant.max
 sizerun <- fuC(cvrun)

}

return(list(
"critical.value" = cvrun,
"approximate.size" = sizerun$size,
"iter" = it
 )
)

}