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

test.stat <- function(
y,	 		      #Observed vector(s) (n times l matrix)
R, 				    #restriction matrix (q times k, rank q)
r,            #right-hand side in restriction (q-vector)
X, 				    #design matrix (n times k, rank k)
hcmethod,   	#-1:4; -1 = classical F-test without df adjustment; 0 = HC0, 1 = HC1, etc.
restr.cov,  	#Covariance Matrix estimator computed with restricted residuals (TRUE or FALSE)
tol = 1e-08,  #tolerance parameter used in checking invertibility of VCESTIMATOR in test statistic
cores = 1		  #number of cores used in computations
){

#transform y to a matrix in case it is a vector

if(is.vector(y) == TRUE){
y <- as.matrix(y, nrow = length(y))
} 

###########################################################################
#run input checks
###########################################################################

check.y(y, X)
check.X.R(X, R)
check.r(r, R)
check.hcmethod(hcmethod)
check.restr.cov(restr.cov)
check.tol(tol)
check.cores(cores)

###########################################################################
# Elementary quantities
###########################################################################

n <- dim(X)[1]                        #sample size
k <- dim(X)[2]                        #number of regressors
q <- length(r)                        #number of restrictions
qrX <- qr(X)                          #qr decomposition of X

###########################################################################
# prepare the input for F.wrap
###########################################################################

if(restr.cov){
 qrM0lin <- qr(M0lin(X, R))           #qr decomp of basis of M0lin = M0-mu0
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
# apply F.wrap to obtain the output
###########################################################################

test.val <- F.wrap(y, R, r, X, n, k, q, qrX, hcmethod, cores, 
	restr.cov, Bfac, Bfac2, qrM0lin, RF.loc, tol)
	
return(list("test.val" = test.val))
}