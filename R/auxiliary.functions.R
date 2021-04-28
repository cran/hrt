#
#    Copyright (C) 2021 David Preinerstorfer
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

###########################################################################
###########################################################################
#
# Wrapper function for test statistic of the form (in case y is vector):
#
# ``(R\hat{\beta}(y))' VCESTIMATOR^{-1} (R\hat{\beta}(y))''
#
# here \hat{\beta} denotes the OLS estimator
# and VCESTIMATOR is a covariance matrix estimator that might be based
# on restricted or unrestricted residuals (restriction: R\beta = r);
# (the choice of restricted or unrestricted residuals also has an effect
# on the weights in the construction of the HC0 - HC4 estimators)
#
###########################################################################
###########################################################################

  F.wrap <- function(
  y,          #matrix (n rows) of observations 
  R,          #restriction matrix (q times k, rank q)
  r,          #right-hand side in restriction (q-vector)
  X,          #design matrix (n times k, rank k)
  n,          #sample size
  k,          #number of columns of X
  q,          #number of rows of R
  qrX,        #qr decomposition of X
  hcmethod,   #-1:4; -1 = classical F-test without df adjustment; 0 = HC0, 1 = HC1, etc.
  cores,      #number of cores used in computation
  restr.cov,  #Covariance Matrix estimator computed with restricted residuals (TRUE or FALSE)  
  Bfac,       #R(X'X)^{-1}X' (cf. function Bfactor.matrix below)
  Bfac2,      #(R(X'X)^(-1)R')^(-1) (cf. function Bfactor.matrix2 below)
  qrM0lin,    #qr decomposition of M0lin (cf. function M0lin below)
  RF,         #RF = (X'X)^(-1)R'(R(X'X)^(-1)R')^(-1) used in res.OLSRF (cf. below)
  tol         #tolerance parameter used in checking invertibility of VCESTIMATOR in test statistic
              #nonpositive numbers are reset to tol = machine epsilon (.Machine$double.eps)
              #if VCESTIMATOR is not invertible, -1 is returned 
              #if VCESTIMATOR is invertible, the test statistic is always nonnegative
  ){

  #computation of residuals and weights used in the construction of the covariance matrix estimators

  if(restr.cov == TRUE){
   umat <- res.OLSRF(y, qrX, R, r, RF)$Res.res
   Wmat <- wvec(k-q, n, qrM0lin, hcmethod)
  }

  if(restr.cov ==  FALSE){
   umat <- qr.resid(qrX, y)
   Wmat <- wvec(k, n, qrX, hcmethod)
  }

  Rbmat <- R %*% qr.coef(qrX, y) - matrix(r, byrow = FALSE, nrow = q, ncol = dim(y)[2])

  if(hcmethod == -1){
   if(tol <= 0){
   tol <- .Machine$double.eps
   }
   df <- n-(k - restr.cov*q)
   var.est <- apply(umat^2,2,sum)/df
   nonzero <- (var.est > tol)
   test.val <- nonzero*1
   CBFac <- chol(Bfac2)
   Wmat <- CBFac%*%Rbmat
   Wmat <- apply(Wmat^2, 2, sum)
   test.val[test.val == 1] <- (Wmat[test.val == 1]/var.est[test.val == 1])
   test.val[test.val == 0] <- 0
  } else {
  test.val <- .Call('_hrt_ctest', PACKAGE = 'hrt', umat, Rbmat, Wmat, Bfac, cores, tol)
  }

  return(test.val)
  }

###########################################################################
###########################################################################
#
# In case q = 1 the rejection event can be characterized by a quadratic
# form being non-negative (under Assumptions). The following function 
# returns the matrix defining this quadratic form.
#
###########################################################################
###########################################################################


  A.mat <-   function(
  C,	 		  #critical value
  R, 				#restriction matrix (1 times k, rank q)
  X, 				#design matrix (n times k, rank k)
  hcmethod, #-1:4; -1 = classical F-test without df adjustment; 0 = HC0, 1 = HC1, etc.
  restr.cov #Covariance Matrix estimator computed with restricted residuals (TRUE or FALSE)
  ){

  ###########################################################################
  #run input checks 
  ###########################################################################

  check.C(C)
  check.X.R(X, R)

  if(dim(R)[1] != 1){
   stop("R needs to be 1xk dimensional for A.mat")
  }

  check.hcmethod(hcmethod)
  check.restr.cov(restr.cov)

  ###########################################################################
  # Elementary quantities
  ###########################################################################

  n <- dim(X)[1]                        #sample size
  k <- dim(X)[2]                        #number of regressors
  q <- dim(R)[1]                        #number of restrictions
  rloc <- rep(0, length = q)            #set r = 0 (result independent of r)
  qrX <- qr(X)

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
  # quadratic form
  ###########################################################################

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

  return(list("A" = MQF))
  }

###########################################################################
###########################################################################
#
# Check if a given design matrix X satisfies Assumption 1 or 2
# whether Assumption 1 or 2 is checked depends on the input qr decomposition
#
# if the qr decomposition of X is provided, Assumption 1 is checked
# if the qr decomposition of a basis of M0lin is provided, Assumption 2 is
# checked
# 
# Bfac = R(X'X)^{-1} X' (cf. the function Bfactor.matrix below)
#
###########################################################################
###########################################################################

  As.check <- function(
  qrXM,   #qr decomposition X (n times k, rank k), or of a basis of 
          #M0lin (n times (k-q), rank (k-q), cf. function M0lin below)
          #in the function we abbreviate ``k-q'' by l
  Bfac,   #R(X'X)^{-1}X' (cf. function Bfactor.matrix below)
  q,      #number of rows of restriction matrix R
  as.tol  #tolerance parameter used in checking Assumptions 1 or 2, respectively
          #nonpositive numbers are reset to as.tol = machine epsilon (.Machine$double.eps)
  ){

  XM <- qr.X(qrXM)
  n <- dim(XM)[1]
  l <- dim(XM)[2]

  #if l == 0, the Assumption to be checked automatically satisfied

  if(l == 0){
  return(TRUE)
  } else {

  #if l != 0, do the following:

    e0 <- matrix(0, nrow = n, ncol = 1)  
    ind <- rep(NA, length = n) 
    
    #check whether elements of the canonical basis are in the span of XM

    for(i in 1:n){
      ei <- e0
      ei[i,1] <- 1
      ind[i] <- ( rrank(cbind(XM, ei), as.tol) > l )
    }
    
  #if all entries of ind are TRUE, the Assumption is automatically
  #satisfied, because of the rank assumptions on R and X
  #if there are entries of ind that are FALSE, the Assumption
  #has to be checked via a rank computation

    if(sum(ind) == n){
     return ( TRUE )
    } else {
     return( rrank(Bfac[, ind, drop=FALSE], as.tol) == q )
    }
  }
  }


###########################################################################
###########################################################################
#
# Compute the 
#
#   rank of a matrix based on the Eigen FullPivLU decomposition
#
# in case input is a matrix with one column or with one row, the 
# function checks if the maximal absolute entry of the input exceeds tol 
#
############################################################################
############################################################################

  rrank <- function(
  A,    #input matrix
  tol   #tolerance parameter passed to Eigen FullPivLU decomposition;
        #if tolerance parameter is nonpositive, tol is reset to machine epsilon
  ){
  if(min(dim(A)) == 1){
    if(tol <= 0){
    tol <- .Machine$double.eps
    }
   (max(abs(A)) > tol)*1
  } else {
   .Call('_hrt_rrank', PACKAGE = 'hrt', A, tol)
  }
  }

###########################################################################
###########################################################################
#
# Compute the 
#
#     kernel of a matrix based on the Eigen FullPivLU decomposition
#
###########################################################################
###########################################################################

  rkernel <- function(
  A,  #input matrix
  tol #tolerance parameter passed to Eigen FullPivLU decomposition;
      #if tolerance parameter is nonpositive, tol is reset to machine epsilon
  ){
  .Call('_hrt_rkernel', PACKAGE = 'hrt', A, tol)
  }

###########################################################################
###########################################################################
#
# Basis of $M_0^{lin} = \{y \in \mathbb{R}^n: y = Xb, Rb = 0\}$
#
#Computation is based on function rkernel with tol = -1 
#
# In case q = k, i.e., M_0^{lin} = (0), output vector of dimension n times 0
#
###########################################################################
###########################################################################

  M0lin <- function(
  X,  #input matrix (n times k, rank k)
  R   #restriction matrix (q times k, rank q)
  ){
  if(dim(R)[1] == dim(X)[2]){
    return(X[,-(1:dim(X)[2])]) 
  } else {
    X%*%rkernel(R, -1)
  }
  }

###########################################################################
###########################################################################
#
# Compute the matrix
#
#     R(X'X)^(-1)X' 
#
###########################################################################
###########################################################################

  Bfactor.matrix <- function(
  qrX,        #qr decomposition of X (n times k, rank k)
  n,          #number of rows of X
  R = FALSE   #restriction matrix R (q times k, rank q)
              #if R = FALSE, then (X'X)^(-1)X' is returned
  ){
  A <- qr.coef(qrX, diag(n))
  if(is.matrix(R) == TRUE){
  return(R%*%A)
  } else {
  return(A)
  }
  }

###########################################################################
###########################################################################
#
# Compute the matrix 
#
#     (R(X'X)^(-1)R')^(-1) 
#
###########################################################################
###########################################################################

  Bfactor.matrix2 <- function(
  qrX,      #qr decomposition of X (n times k, rank k)
  R         #restriction matrix R (q times k, rank q)
  ){
  factor.tmp <- qr.R(qrX)
  factor.tmp <- R%*% backsolve(qr.R(qrX), diag(dim(factor.tmp)[1]))
  Bfactor2 <- solve(tcrossprod(factor.tmp))
  return(Bfactor2)
  }

###########################################################################
###########################################################################
#
# Generate the 
#
#     weights d_i (unrestricted case) or \tilde{d}_i (restricted case)
#
# used in the construction of the covariance estimators H0 - H4 
#
# Whether d_i or \tilde{d}_i is generated is via qrXM, which determines the
# hat matrix used
#
###########################################################################
###########################################################################

  wvec <- function(
  l,        #see description of qrXM
  n,        #see description of qrXM
  qrXM,     #qr decomposition of matrix XM (n times l, rank l)
            #weights are based on the diag of hat matrix corresp. to XM
  hcmethod  #0-4 (HC0, HC1, HC2, HC3, HC4)
  ){
  if(hcmethod == 0 | l == 0){
   return(rep(1, length = n))
  }

  if(hcmethod == 1){
   return(rep(n/(n-l), length = n))
  }

  Q <- qr.Q(qrXM)
  H <- diag(tcrossprod(Q,Q))
  H[H == 1] <- 0
  h <- 1/(1-H)

  if(hcmethod == 2){
   return(h)
  }

  if(hcmethod == 3){
   return(h^2)
  }

  if(hcmethod == 4){
   delta <-  sapply(n*H/l, function(x) {min(4, x)})
   return(h^delta)
  }

  }

###########################################################################
###########################################################################
#
# Function that computes 
#
#     restricted OLS estimators and restricted residuals 
#
###########################################################################
###########################################################################

  res.OLSRF <- function(
  y,     #matrix (n rows) of observations 
  qrX,   #qr decomposition of X (n times k, rank k)
  R,     #restriction matrix (q times k, rank q)
  r,     #right-hand side in the restriction (q-vector)
  RF     #RF = (X'X)^(-1)R'(R(X'X)^(-1)R')^(-1)
  ){

  OLS.coef <- qr.coef(qrX, y)
  OLS.coefs.restr <- OLS.coef - RF%*%R%*%OLS.coef
  OLS.coefs.restr <- OLS.coefs.restr + matrix(RF%*%r, nrow = dim(R)[2], ncol = dim(y)[2])
  OLS.res.restr <- y - qr.X(qrX)%*%OLS.coefs.restr

  return(list("Coef.restr" = OLS.coefs.restr, 
      "Res.res" = OLS.res.restr))
  }

###########################################################################
###########################################################################
#
# Checks if a given design matrix X satisfies the sufficient conditions for
# the existence of a size-controlling critical value
#
# Bfac = R(X'X)^{-1} X' (cf. the function Bfactor.matrix below)
#
# Returns TRUE if a size-controlling CV exists, FALSE else.
#
###########################################################################
###########################################################################

  SC.check <- function(
  X,            #X
  R,            #R
  hcmethod,   	#-1:4; -1 = classical F-test without df adjustment; 0 = HC0, 1 = HC1, etc.
  restr.cov,  	#Covariance Matrix estimator computed with restricted residuals (TRUE or FALSE)
  as.tol,        #tolerance parameter used in checking rank conditions
  checks = TRUE #input checks are run if TRUE
  ){

  #input checks

  if(checks == TRUE){
  check.X.R(X, R)
  check.hcmethod(hcmethod)
  check.restr.cov(restr.cov)
  check.as.tol(as.tol)
  }

  #check sufficient condition

  if( (hcmethod == -1) & (restr.cov == TRUE) ){
   return(TRUE)
  } else {

   n <- dim(X)[1]
   k <- dim(X)[2]
   q <- dim(R)[1]
   e0 <- matrix(0, nrow = n, ncol = 1)  
   XM <- M0lin(X, R)
   l <- dim(XM)[2]
   qrX <- qr(X)
   
   Bfac <- Bfactor.matrix(qrX, n, R)
    
   #check condition for every index such that ei notin span(XM)

   v <- c()
   
   if( hcmethod > -1 ){
   
    if( restr.cov ) {
     qrXM <- qr(XM)  
    } else {
     qrXM <- qrX
    }
   
    for(i in 1:n){
     ei <- e0
     ei[i,1] <- 1

     if( rrank(cbind(XM, ei), as.tol) > l ){
      A <- Bfac %*% diag(c(qr.resid(qrXM, ei)))
      v <- c(v, rrank(A, as.tol))
     }
    }
   return(!( length(v[v < q])>0 ))
   } 
   

   if( (hcmethod == -1) & (restr.cov == FALSE) ){
   
    for(i in 1:n){
     ei <- e0
     ei[i,1] <- 1

     if( rrank(cbind(XM, ei), as.tol) > l ){
      A <- cbind(X, ei)
      v <- c(v, rrank(A, as.tol))
     }
    }
    return(!( length(v[v < k])>0 ))
    }
   
   }
  }

###########################################################################
###########################################################################
#
# Lower bound for size-controlling critical value 
#
###########################################################################
###########################################################################

  LB.cv <- function(
  X,            #X
  R,            #R
  hcmethod,   	#-1:4; -1 = classical F-test without df adjustment; 0 = HC0, 1 = HC1, etc.
  restr.cov,  	#Covariance Matrix estimator computed with restricted residuals (TRUE or FALSE)
  as.tol,        #tolerance parameter used in checking rank conditions
  checks = TRUE #input checks are run if TRUE
  ){

  #input checks
  if(checks == TRUE){
  check.X.R(X, R)
  check.hcmethod(hcmethod)
  check.restr.cov(restr.cov)
  check.as.tol(as.tol)
  }

  #basic quantities

  n <- dim(X)[1]
  k <- dim(X)[2]
  q <- dim(R)[1]
  e0 <- matrix(0, nrow = n, ncol = 1)  
  XM <- M0lin(X, R)
  l <- dim(XM)[2]

  #compute lower bound (evaluate test statistic
  #at all points e_i in I_1, and take the maximum
  #of those values)

  y <- c()

  for(i in 1:n){
   ei <- e0
   ei[i,1] <- 1
   
   if( rrank(cbind(XM, ei), as.tol) > l ){
    y <- cbind(y, ei)
   }
   
  }

  if(!is.matrix(y)){
   stop("The set I1(M0lin) seems to be the whole space ... check your inputs!")
   } else {
   m <- test.stat(y, R, rep(0, length = q), 
   X, hcmethod, restr.cov)$test.val
   return(max(m))
  }

  }

###########################################################################
###########################################################################
#
# Input checks 
#
###########################################################################
###########################################################################

  check.alpha <- function(alpha){
  if( !is.numeric(alpha) | alpha <= 0 | alpha >= 1 ) {
   stop("Invalid 'alpha' value - 'alpha' must be in the interval (0,1)")
  }
  }

  check.N.M <- function(N0, N1, N2, Mp, M1, M2){

  if( N0%%1 != 0 | N0 < 0 ){
   stop("Invalid 'N0' value - 'N0' must be a positive integer")
  }

  if( N1%%1 != 0 | N1 < 0 ){
   stop("Invalid 'N1' value - 'N1' must be a positive integer")
  }

  if( N2%%1 != 0 | N2 < 0 ){
   stop("Invalid 'N2' value - 'N2' must be a positive integer")
  }

  if( Mp%%1 != 0 | Mp < 0 ){
   stop("Invalid 'Mp' value - 'Mp' must be a positive integer")
   if(Mp < 100){
   warning("The chosen value for 'Mp' seems very low, you may want to 
   increase 'Mp' for more accurate results")
   }
  }

  if( M1%%1 != 0 | M1 < 0 ){
   stop("Invalid 'M1' value - 'M1' must be a positive integer")
   if(M1 < 10){
   warning("The chosen value for 'M1' seems very low, you may want to 
   increase 'M1' for more accurate results")
   }
  }

  if( M2%%1 != 0 | M2 < 0 ){
   stop("Invalid 'M2' value - 'M2' must be a positive integer")
  }

  if( N0 > N1 ) {
   warning("'N1' should be greater than 'N0'")
  }

  if( N1 > N2 ) {
   warning("'N2' should be greater than 'N1'")
  }

  if( M1 > Mp ) {
   stop("Invalid 'M1' value - 'M1' can not be greater than 'Mp'")
  }

  if( M2 > M1 ) {
   stop("Invalid 'M2' value - 'M2' can not be greater than 'M1'")
  }
                  
  }


  check.M <- function(Mp, M1, M2){

  if( Mp%%1 != 0 | Mp < 0 ){
   stop("Invalid 'Mp' value - 'Mp' must be a positive integer")
   if(Mp < 100){
   warning("The chosen value for 'Mp' seems very low, you may want to 
   increase 'Mp' for more accurate results")
   }
  }

  if( M1%%1 != 0 | M1 < 0 ){
   stop("Invalid 'M1' value - 'M1' must be a positive integer")
   if(M1 < 10){
   warning("The chosen value for 'M1' seems very low, you may want to 
   increase 'M1' for more accurate results")
   }
  }
  if( M2%%1 != 0 | M2 < 0 ){
   stop("Invalid 'M2' value - 'M2' must be a positive integer")
  }

  if( M1 > Mp ) {
   stop("Invalid 'M1' value - 'M1' can not be greater than 'Mp'")
  }

  if( M2 > M1 ) {
   stop("Invalid 'M2' value - 'M2' can not be greater than 'M1'")
  }
                  
  }

  check.y <- function(y, X){
  if( !is.numeric(y) | !is.matrix(y) | dim(y)[2] == 0 | dim(y)[1] != dim(X)[1]){
   stop("Invalid 'y' value - 'y' must either be a real vector of length 
   dim(X)[1], or a real matrix with dim(X)[1] rows with more than 0 columns")
  }
  }

  check.X.R <- function(X, R){

  if( !is.matrix(X) ){
   stop("Invalid 'X' value - must be a matrix")
  }

  if( !is.matrix(R) ) {
   stop("Invalid 'R' value - must be a matrix")
  }

  if( dim(X)[2] >= dim(X)[1] ){
   stop("Number of columns of 'X' is not smaller than its number of rows")
  }

  if( dim(X)[2] == 0 ){
   stop("Number of rows of 'X' must be greater than 0")
  }

  if( dim(X)[2] != dim(R)[2] ) {
   stop("Matrices 'X' and 'R' have different numbers of columns")
  }

  if( rrank(X, -1) < dim(X)[2] ) {
   stop("The matrix 'X' is numerically of rank < k")
  }

  if( rrank(R, -1) < dim(R)[1] ) {
   stop("The matrix 'R' is numerically of rank < q")
  }

  }

  check.r <- function(r, R){
  if( !is.vector(r) | !is.numeric(r) | (length(r) != dim(R)[1])){
   stop("Invalid 'r' value - 'r' must be a real vector
   the length of which coincides with the number of rows of R") 
  }
  }

  check.hcmethod <- function(hcmethod){
  if(length(hcmethod) != 1 | sum(c(-1:4) == hcmethod) != 1) {
   stop("Invalid 'hcmethod' value - 'hcmethod' must be -1, 0, 1, 2, 3, or 4")
  }
  }

  check.restr.cov <- function(restr.cov){
   if(!is.logical(restr.cov) | length(restr.cov) != 1){
   stop("Invalid 'restr.cov' value - 'restr.cov' must be TRUE or FALSE")
  }
  }

  check.tol <- function(tol){
  if(!is.numeric(tol) | length(tol) != 1){
   stop("Invalid 'tol' value - 'tol' must be a real number")
  }
  if(tol > .1){
  warning("Tolerance parameter 'tol' should typically be chosen small, e.g., 1e-07; 
  your choice seems unusually large, and might yield wrong results.")
  }
  if(tol <= 0){
  warning("Your non-positive 'tol' parameter was converted to tol = machine epsilon")
  }
  }

  check.as.tol <- function(as.tol){
  if(!is.numeric(as.tol) | length(as.tol) != 1){
   stop("Invalid 'as.tol' value - 'as.tol' must be a real number")
  }
  if(as.tol > .1){
  warning("Tolerance parameter 'as.tol' should typically be chosen small, e.g., 1e-07; 
  your choice seems unusually large, and might yield wrong results.")
  }
  if(as.tol <= 0){
  warning("Your non-positive 'as.tol' parameter was converted to as.tol = machine epsilon")
  }
  }

  check.in.tol <- function(in.tol){
  if(!is.numeric(in.tol) | length(in.tol) != 1){
   stop("Invalid 'in.tol' value - 'in.tol' must be a real number")
  }
  if(in.tol > .1){
  warning("Tolerance parameter 'in.tol' should typically be chosen small, e.g., 1e-07; 
  your choice seems unusually large, and might yield wrong results.")
  }
  if(in.tol <= 0){
  warning("Your non-positive 'in.tol' parameter was converted to in.tol = machine epsilon")
  }
  }
  
  check.size.tol <- function(size.tol){
  if(!is.numeric(size.tol) | length(size.tol) != 1){
   stop("Invalid 'size.tol' value - 'size.tol' must be a real number")
  }
  if(size.tol > .1){
  warning("Tolerance parameter 'size.tol' should typically be chosen small; 
  your choice seems large, and might yield wrong results.")
  }
  if(size.tol <= 0){
  warning("Your non-positive 'size.tol' parameter was converted to size.tol = .0001")
  }
  }
  
  check.eps.close <- function(eps.close){
  if(!is.numeric(eps.close) | length(eps.close) != 1){
   stop("Invalid 'eps.close' value - 'in.tol' must be a real number")
  }
  if(eps.close > .2){
  warning("Parameter 'eps.close' should typically be chosen small, e.g., 0.001; 
  your choice seems unusually large, and might yield wrong results.")
  }
  if(eps.close <= 0){
  warning("Your non-positive 'eps.close' parameter was converted to eps.close = 0.001")
  }
  }


  check.checks <- function(checks){
  if(!is.logical(checks) | length(checks) != 1){
   stop("Invalid 'checks' value - 'checks' must be TRUE or FALSE")
  }
  }

  check.cores <- function(cores){
  if( cores%%1 != 0 | cores <= 0) {
   stop("Invalid 'cores' value - 'cores' must be a positive integer")
  }
  }

  check.C <- function(C){
  if( !is.numeric(C) | length(C)!=1 | C <= 0) {
   stop("Invalid 'C' value - 'C' must be a positive real number")
  }
  }

  check.lower <- function(lower, n){
  if(lower < 0 | lower >= n^(-1) ){
   stop("Invalid 'lower' value - must be in (0, n^(-1))")
   }
  }

  check.levelCl <- function(levelCl){
  if( !is.numeric(levelCl) | length(levelCl)!=1 | levelCl < 0 | levelCl >= 1) {
   stop("Invalid 'levelCl' value - 'levelCl' must be a number in [0, 1)")
  }
  }
  
  check.LBcheck <- function(LBcheck){
  if( !is.logical(LBcheck) ){
   stop("Invalid 'LBcheck' value - 'LBcheck' must be TRUE or FALSE")
  }
  }
  
  check.maxit <- function(maxit){
  if( maxit%%1 != 0 | maxit <= 0) {
   stop("Invalid 'maxit' value - 'maxit' must be a positive integer")
  }
  }

  check.eps.close.lower <- function(eps.close, lower, n){
  if(lower != 0){
   if(lower + eps.close >= 1/(n-1)){
    stop("Invalid 'lower' and 'eps.close' combination. lower + eps.close
    must be smaller than 1/(n-1)")
   }
  }
  }