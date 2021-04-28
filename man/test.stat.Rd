\name{test.stat}
\alias{test.stat}
\encoding{latin1}
\title{Computation of the test statistics}
\description{This function computes the test statistics \eqn{T_{uc}}, \eqn{T_{Het}} 
(with HC0-HC4 weights), \eqn{\tilde{T}_{uc}}, or \eqn{\tilde{T}_{Het}} 
(with HC0R-HC4R weights) as defined in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021).}

\usage{
test.stat(y, R, r, X, hcmethod, restr.cov, tol = 1e-08, cores = 1)
}

\arguments{
\item{y}{Either an observation vector, or a matrix the columns of which are
observation vectors. The
number of rows of an observation vector must coincide with the number of rows
of the design matrix \code{X}.}

\item{R}{The restriction matrix. \code{test.stat} computes a test 
statistic for the hypothesis \eqn{R \beta = r}. \code{R} needs to be of
full row rank, and needs to have the same number of columns as \code{X}.}

\item{r}{The restriction vector. \code{test.stat} computes a test 
statistic for the hypothesis \eqn{R \beta = r}. \code{r} needs to be a
vector with the same number of coordinates as the number of rows of \code{R}.}
                
\item{X}{The design matrix \code{X} needs to be of full column rank. 
The number of columns of \code{X} must be smaller than the number of rows of 
\code{X}.}

\item{hcmethod}{Integer in [-1, 4]. Determines the method applied in the construction of the covariance estimator
used in the test statistic. The value -1 corresponds to unadjusted (i.e., classical) F statistic without df adjustment; the value 0 
corresponds to the HC0 estimator; ...; the value 4 corresponds to the HC4 estimator. Note that in case \code{restr.cov}
is TRUE the null-restricted versions of the covariance estimators are computed. Cf. \enc{Pötscher}{Poetscher} 
and Preinerstorfer (2021) and the references there for details.}

\item{restr.cov}{TRUE or FALSE. Covariance matrix estimator based on null-restricted (TRUE) or unrestricted (FALSE) residuals.}

\item{tol}{(Small) positive real number. Tolerance parameter used in checking invertibility of the covariance matrix in the test statistic. 
Default is 1e-08.}

\item{cores}{The number of CPU cores used in the (parallelized)
computations. Default is 1. Parallelized computation is enabled only
if the compiler used to build \pkg{hrt} supports OpenMP.}
}

\details{For details concerning the test statistics please
see the relevant sections in \enc{Pötscher}{Poetscher} and Preinerstorfer (2021)
.
}

\value{The function returns a list consisting of:
\item{test.val}{Either a vector the entries of which correspond to the values 
of the
test statistic evaluated at each column of the input matrix \code{y}, or, 
if \code{y} is a vector, the test statistic evaluated at \code{y}.}
}

\examples{
n <- 20
y <- rnorm(n)
X <- cbind(rep(1, length = n), rnorm(n))
R <- matrix(1, nrow = 1, ncol = 2)
r <- 0
hcmethod <- 4
restr.cov <- FALSE
test.stat(y, R, r, X, hcmethod, restr.cov)
}

\references{
\enc{Pötscher}{Poetscher}, B. M. and Preinerstorfer, D. (2021). Valid Heteroskedasticity Robust Testing. <arXiv:2104.12597>
}
