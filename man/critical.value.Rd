\name{critical.value}

\alias{critical.value}

\title{Critical Values for Heteroskedasticity Robust Testing}
\encoding{latin1}
\description{
This function provides an implementation of 
Algorithm 3 in \enc{Pötscher}{Poetscher} and Preinerstorfer (2021), based on 
Algorithm 1 (if \eqn{q = 1}) or Algorithm 2 
(if \eqn{q > 1}) in the same reference as the auxiliary algorithm \eqn{\mathsf{A}}.
Which of the two algorithms is used is automatically 
determined as a function of \eqn{q}, the number of rows of \code{R}.
 
The user is referred to
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021) for definitions, 
a detailed description of the problems solved
by the algorithms, and for a detailed description of the algorithms themselves.
 
Most of the input parameters to \code{critical.value} are actually used
in the auxiliary Algorithm 1 or 2, respectively. 
 Algorithm 1 is based on the function 
\code{\link{davies}} from the package \pkg{CompQuadForm}. The parameters
\code{lim} and \code{acc} for \code{\link{davies}} can be supplemented by the user.
Algorithms 1 and 2 are implemented using the function \code{\link{constrOptim}} 
from \pkg{stats} in Stages 1 and 2; this function
is used with default parameters, but control parameters can be supplied by the user.

After determining a critical value for a given testing 
problem via the function \code{critical.value}, it is recommended that: (i) the user
applies the function \code{size} to compute the size of the test corresponding to the critical value obtained; 
and (ii) to check whether the size obtained does coincide with (or is close to) the targeted level of 
significance (that is \code{alpha}). If (ii) is not the case, this is an indication
of numerical issues, which potentially can be avoided by changing the input parameters
responsible for the accuracy of the computations.
}

\usage{
critical.value(alpha, R, X, hcmethod, restr.cov, Mp, M1, M2, 
N0 = NULL, N1 = NULL, N2 = NULL, tol = 1e-08, 
control.1 = list("reltol" = 1e-02, "maxit" = dim(X)[1]*20),
control.2 = list("reltol" = 1e-03, "maxit" = dim(X)[1]*30),
cores = 1, lower = 0, eps.close = .0001, lim = 30000, acc = 0.001, 
size.tol = .001, maxit = 25, as.tol = 1e-08)
}

\arguments{
\item{alpha}{Significance level. A real number in the interval \eqn{(0, 1)}.}

\item{R}{The restriction matrix. \code{critical.value} computes the (smallest) size-controlling critical value
for a test of the hypothesis \eqn{R \beta = r}. \code{R} needs to be of
full row rank, and needs to have the same number of columns as \code{X}.}
                 
\item{X}{The design matrix \code{X} needs to be of full column rank. 
The number of columns of \code{X} must be smaller than the number of rows of 
\code{X}.}

\item{hcmethod}{Integer in [-1, 4]. Determines the method applied in the construction of the covariance estimator
used in the test statistic. The value -1 corresponds to unadjusted (i.e., classical) F statistic without df adjustment; the value 0 
corresponds to the HC0 estimator; ...; the value 4 corresponds to the HC4 estimator. Note that in case \code{restr.cov}
is TRUE the null-restricted versions of the covariance estimators are computed. Cf. \enc{Pötscher}{Poetscher} 
and Preinerstorfer (2021) and the references there for details.}

\item{restr.cov}{TRUE or FALSE. Covariance matrix estimator based on null-restricted (TRUE) or unrestricted (FALSE) residuals.}

\item{Mp}{This input is used in Algorithm 1 or 2, respectively. \code{Mp} is a positive integer 
(should be chosen large, e.g., 50000; but the feasibility depends on the dimension of \code{X}, etc). 
\code{Mp} determines \eqn{M_0} in Algorithm 1 or 2 (i.e., \eqn{\mathsf{A}}), respectively, that is, 
the number of initial values chosen in Stage 0 of that algorithm. The way initial values  
(i.e., the sets of variance covariance matrices \eqn{\Sigma_j} in Stage 0 of the algorithm; 
 the diagonal entries of each \eqn{\Sigma_j} sum up to 1) are chosen is as follows:

\enumerate{
\item If \eqn{q = 1} and \eqn{lower = 0}, one of the initial values \eqn{\Sigma_j} is a matrix which maximizes the expectation of the quadratic form \eqn{y \mapsto y'\Sigma^{1/2}
A_C \Sigma^{1/2}y} under an n-variate standard normal distribution. Here, \eqn{A_C} is a matrix that is defined \enc{Pötscher}{Poetscher} and Preinerstorfer (2021).
If diagonal entries of this maximizer are 0, then they are replaced by the value of \code{eps.close} (and the other values are adjusted
so that the diagonal sums up to 1).

\item One starting value \eqn{\Sigma_j} is a diagonal matrix with constant diagonal entries.

\item If \code{lower} is zero, then (i) \eqn{\lceil Mp/4 \rceil - 1} covariance matrices \eqn{\Sigma_j} are drawn by sampling their diagonals \eqn{\tau_1^2, ..., \tau_n^2} from a uniform distribution 
on the unit simplex in \eqn{R^n}; and (ii)  the remaining \eqn{M_p - (\lceil Mp/4 \rceil - 1)} covariance matrices \eqn{\Sigma_j} are each drawn by first sampling a vector \eqn{(t_1, ..., t_n)'} from a uniform distribution 
on the unit simplex in \eqn{R^n}, and by then obtaining the diagonal \eqn{\tau_1^2, ..., \tau_n^2} of \eqn{\Sigma_j} via \eqn{(t_1^2, ..., t_n^2)/\sum_{i = 1}^n t_i^2}. If \code{lower} is nonzero,
then the initial values are drawn analogously, but from a uniform distribution on the subset of the unit simplex in \eqn{R^n} corresponding to the restriction imposed
by the lower bound \code{lower}.

\item \eqn{n} starting values equal to covariance matrices with a single dominant diagonal entry and all other diagonal entries constant. The size of the dominant
diagonal entry is regulated via the input parameters \code{eps.close} and \code{lower}. In case \code{lower} is nonzero,
the size of the dominant diagonal entry equals \eqn{1-(n-1)*(lower + eps.close)}. In case \code{lower} is zero, the size of the dominant diagonal
entry equals \eqn{1-eps.close}.

}
}

\item{M1}{This input is used in Algorithm 1 or 2, respectively. A positive integer 
(should be chosen large, e.g., 500; but the feasibility depends on the dimension of \code{X}, etc). 
Corresponds to \eqn{M_1} in the description of 
Algorithm 1 and 2 in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021). M1 must not exceed Mp.}

\item{M2}{This input is used in Algorithm 1 or 2, respectively.
A positive integer. Corresponds to \eqn{M_2} in the description of 
Algorithm 1 and 2 in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021). M2 must not exceed M1.}

\item{N0}{This input is needed in Algorithm 2.
Only used in case \eqn{q > 1} (i.e., when Algorithm 2 is used). A positive integer. 
Corresponds to \eqn{N_0} in the description of Algorithm 2 in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021).}

\item{N1}{This input is needed in Algorithm 2.
Only used in case \eqn{q > 1} (i.e., when Algorithm 2 is used). A positive integer. 
Corresponds to \eqn{N_1} in the description of Algorithm 2 in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021). N1 should be greater than 
N0.}     

\item{N2}{This input is needed in Algorithm 2.
Only used in case \eqn{q > 1} (i.e., when Algorithm 2 is used). A positive integer. 
Corresponds to \eqn{N_2} in the description of Algorithm 2 in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021). N2 should be greater than
N1.} 

\item{tol}{This input is used in Algorithm 1 or 2, respectively.
(Small) positive real number. Tolerance parameter used in checking invertibility of the 
covariance matrix in the test statistic. Default is 1e-08.}

\item{control.1}{This input is used in Algorithm 1 or 2, respectively.
Control parameters passed to the \code{\link{constrOptim}} function in 
Stage 1 of Algorithm 1 or 2, respectively. Default is \code{control.1 = list("reltol" = 1e-02, 
"maxit" = dim(X)[1]*20)}.}

\item{control.2}{This input is used in Algorithm 1 or 2, respectively.
Control parameters passed to the \code{\link{constrOptim}} function in 
Stage 2  of Algorithm 1 or 2, respectively. Default is \code{control.2 = list("reltol" = 1e-03, 
"maxit" = dim(X)[1]*30)}.}
 
\item{cores}{The number of CPU cores used in the (parallelized)
computations. Default is 1. Parallelized computation is enabled only
if the compiler used to build \pkg{hrt} supports OpenMP.}

\item{lower}{Number in \eqn{[0, n^{-1})} (note that the diagonal of \eqn{\Sigma} is normalized
to sum up to \eqn{1}; if lower > 0, then lower corresponds to what is denoted \eqn{\tau_*} in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021)).
\code{lower} specifies a lower bound on each diagonal entry of the (normalized) covariance matrix
in the covariance model for which the user wants to obtain a critical value that achieves size control.
If this lower bound is nonzero (which is the non-standard choice), then the size is only computed over all covariance matrices,
which are restricted such that their minimal diagonal entry is not smaller than \code{lower}. 
The relevant optimization problems in Algorithm 1 and 2 are then carried out only over this
restricted set of covariance matrices. The size will then in general depend on \code{lower}.
See the relevant discussions concerning restricted heteroskedastic covariance models in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021). Default is \eqn{0}, which is the 
recommended choice, unless there are strong reasons implying a specific lower bound on 
the variance in a given application.
}

\item{eps.close}{(Small) positive real number. This determines the size of the dominant entry in the choice of the
initial values as discussed in the description of the input \code{Mp} above.
Default is 1e-4.}

\item{lim}{This input is needed in Algorithm 1. Only used in case \eqn{q = 1} (i.e., when Algorithm 1 is used).
Input parameter for the function \code{davies}. Default is 30000.}

\item{acc}{This input is needed in Algorithm 1. Only used in case \eqn{q = 1} (i.e., when Algorithm 1 is used).
Input parameter for the function \code{davies}. Default is 1e-3.}

\item{size.tol}{(Small) positive real number. \eqn{\epsilon} in Algorithm 3. Default is 1e-3.}

\item{maxit}{Maximum number of iterations in the while loop of Algorithm 3. Default is 25.}

\item{as.tol}{(Small) positive real number. Tolerance parameter used in checking rank 
conditions for verifying Assumptions 1, 2, and for checking a non-constancy condition
on the test statistic in case \code{hcmethod} is not \eqn{-1} and \code{restr.cov} is
TRUE. \code{as.tol} is also used in the rank computations required for
computing lower bounds for size-controlling critical values. 
Furthermore, \code{as.tol} is used in checking the sufficient conditions
for existence of a size-controlling critical value provided in \enc{Pötscher}{Poetscher} 
and Preinerstorfer(2021). Default is 1e-08.}

}

\details{For details see the relevant sections in \enc{Pötscher}{Poetscher} 
and Preinerstorfer 
(2021), in particular the description of Algorithms 1 and 2 in the Appendix.
}

\value{The output of \code{critical.value} is the following:

\item{critical.value}{The critical value obtained by Algorithm 3.}

\item{approximate.size}{The approximate size of the test based on the 
returned critical value.}

\item{iter}{The number of iterations performed. If \code{iter} is smaller than \code{maxit}, 
then the algorithm determined because the required level of accuracy was achieved.}
}

\examples{

#critical value for the classical (uncorrected) F-test in a location model
#with unrestricted heteroskedasticity

#it is known that (in this very special case) the conventional critical value 
#C <- qt(.975, df = 9)^2
#is size-controlling (thus the resulting size should be 5% (approximately))

R <- matrix(1, nrow = 1)
X <- matrix(rep(1, length = 10), nrow = 10, ncol = 1)
hcmethod <- -1
restr.cov <- FALSE
Mp <- 1000
M1 <- 5
M2 <- 1

#here, the parameters are chosen such that the run-time is low
#to guarantee a high accuracy level in the computation, 
#Mp, M1 and M2 should be chosen much higher

critical.value(alpha = .05, R, X, hcmethod, restr.cov, Mp, M1, M2)
}

\references{
\enc{Pötscher}{Poetscher}, B. M. and Preinerstorfer, D. (2021). Valid Heteroskedasticity Robust Testing. <arXiv:2104.12597>
}

\seealso{
\code{\link{davies}}, \code{\link{constrOptim}}.
}                  