\name{size}

\alias{size}

\title{Computing the Size of Heteroskedasticity Robust Tests}
\encoding{latin1}
\description{
This function provides an implementation of Algorithm 1 (if \eqn{q = 1})  or 2 (if \eqn{q > 1}), respectively, in \enc{Pötscher}{Poetscher} and Preinerstorfer 
(2021). Which of the two algorithms is applied is automatically determined as a function of \eqn{q}. 

The user is referred to the just-mentioned article for definitions, a detailed description of the problem solved
the algorithms, and for a detailed description of the algorithms themselves.

Algorithm 1 is based on the function \code{\link{davies}} from the package \pkg{CompQuadForm}. The parameters
\code{lim} and \code{acc} for \code{\link{davies}} can be supplemented by the user.
Algorithms 1 and 2 are implemented using the function \code{\link{constrOptim}} from \pkg{stats} in Stages 1 and 2; this function
is used with default parameters, but control parameters can be supplied by the user.
}

\usage{
size(C, R, X, hcmethod, restr.cov, Mp, M1, M2, 
N0 = NULL, N1 = NULL, N2 = NULL, tol = 1e-08, 
control.1 = list("reltol" = 1e-02, "maxit" = dim(X)[1]*20),
control.2 = list("reltol" = 1e-03, "maxit" = dim(X)[1]*30),
cores = 1, lower = 0, eps.close = .0001, lim = 30000, acc = 0.001, 
levelCl = 0, LBcheck = FALSE, as.tol = 1e-08)
}

\arguments{
\item{C}{Critical value. A positive real number (for negative critical values
the size of the test equals \eqn{1}).}
                          
\item{R}{The restriction matrix. \code{size} computes the size of a test 
for the hypothesis \eqn{R \beta = r}. \code{R} needs to be of
full row rank, and needs to have the same number of columns as \code{X}.}
                 
\item{X}{The design matrix \code{X} needs to be of full column rank. 
The number of columns of \code{X} must be smaller than the number of rows of 
\code{X}.}

\item{hcmethod}{Integer in [-1, 4]. Determines the method applied in the construction of the covariance estimator
used in the test statistic. The value -1 corresponds to the unadjusted (i.e., classical) F statistic without df adjustment; the value 0 
corresponds to the HC0 estimator; ...; the value 4 corresponds to the HC4 estimator. Note that in case \code{restr.cov}
is TRUE the null-restricted versions of the covariance estimators are computed. Cf. \enc{Pötscher}{Poetscher} 
and Preinerstorfer (2021) and the references there for details.}

\item{restr.cov}{TRUE or FALSE. Covariance matrix estimator based on null-restricted (TRUE) or unrestricted (FALSE) residuals.}

\item{Mp}{A positive integer (should be chosen large, e.g., 50000; but the feasibility depends on the dimension of \code{X}, etc). 
\code{Mp} determines \eqn{M_0} in Algorithm 1 or 2, respectively, that is, 
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
the size of the dominant diagonal entry equals \eqn{1-(n-1)*(lower+eps.close)}. In case \code{lower} is zero, the size of the dominant diagonal
entry equals \eqn{1-eps.close}.

\item If \code{levelCl} is nonzero (see the description of \code{levelCl} below for details concerning this input), then
one further initial value may be obtained by: (i) checking whether \code{C} exceeds 5 times the critical value \eqn{C_H}, say,
for which the rejection probability under homoskedasticity equals \eqn{1-levelCl}; and (ii) if this is the case, running the
function \code{size} (with the same input parameters, but with \code{levelCl} set to \eqn{0} and \code{M2} set to \eqn{1}) on the critical value \eqn{C_H},
and then using the output \code{second.stage.parameter} as a further initial value.
}
}

\item{M1}{A positive integer (should be chosen large, e.g., 500; but the feasibility depends on the dimension of \code{X}, etc). Corresponds to \eqn{M_1} in the description of 
Algorithm 1 and 2 in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021). M1 must not exceed Mp.}

\item{M2}{A positive integer. Corresponds to \eqn{M_2} in the description of 
Algorithm 1 and 2 in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021). M2 must not exceed M1.}

\item{N0}{Only used in case \eqn{q > 1} (i.e., when Algorithm 2 is used). A positive integer. 
Corresponds to \eqn{N_0} in the description of Algorithm 2 in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021).}

\item{N1}{Only used in case \eqn{q > 1} (i.e., when Algorithm 2 is used). A positive integer. 
Corresponds to \eqn{N_1} in the description of Algorithm 2 in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021). N1 should be greater than 
N0.}     

\item{N2}{Only used in case \eqn{q > 1} (i.e., when Algorithm 2 is used). A positive integer. 
Corresponds to \eqn{N_2} in the description of Algorithm 2 in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021). N2 should be greater than
N1.} 

\item{tol}{(Small) positive real number. Tolerance parameter used in checking invertibility of the 
covariance matrix in the test statistic. Default is 1e-08.}

\item{control.1}{Control parameters passed to the \code{\link{constrOptim}}  function in 
Stage 1 of Algorithm 1 or 2, respectively. Default is \code{control.1 = list("reltol" = 1e-02, 
"maxit" = dim(X)[1]*20)}.}

\item{control.2}{Control parameters passed to the \code{\link{constrOptim}}  function in 
Stage 2 of Algorithm 1 or 2, respectively. Default is \code{control.2 = list("reltol" = 1e-03, 
"maxit" = dim(X)[1]*30)}.}
 
\item{cores}{The number of CPU cores used in the (parallelized)
computations. Default is 1. Parallelized computation is enabled only
if the compiler used to build \pkg{hrt} supports OpenMP.}

\item{lower}{Number in \eqn{[0, n^{-1})} (note that the diagonal of \eqn{\Sigma} is normalized
to sum up to \eqn{1}; if lower > 0, then lower corresponds to what is denoted \eqn{\tau_*} in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021)).
\code{lower} specifies a lower bound on each diagonal entry of the (normalized) covariance matrix
in the covariance model for which the user wants to compute the size.
If this lower bound is nonzero, then the size is only computed over all covariance matrices,
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

\item{levelCl}{Number in \eqn{[0, 1)}. This enters via the choice of the initial values as discussed
in the input \code{Mp} above. \code{levelCl} should be used in case \code{C}
is unusually large. In this case, the additional set of starting values
provided may help to increase the accuracy of the size computation. Default is 0.}
 
\item{LBcheck}{Either FALSE (default) or TRUE. If TRUE, then \code{C} is compared
to the theoretical lower bounds on size-controlling critical values in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021). If the supplemented \code{C}
is smaller than the respective lower bound, theoretical results imply that the
size equals 1 and the function \code{size} is halted.}

\item{as.tol}{(Small) positive real number. Tolerance parameter used in checking rank 
conditions for verifying Assumptions 1, 2, and for checking a non-constancy condition
on the test statistic in case \code{hcmethod} is not \eqn{-1} and \code{restr.cov} is
TRUE. Furthermore, \code{as.tol} is used in the rank computations required for
computing lower bounds for size-controlling critical values (in case \code{LBcheck} is
TRUE or \code{levelCl} is nonzero). Default is 1e-08.}
}

\details{For details see the relevant sections in \enc{Pötscher}{Poetscher} 
and Preinerstorfer 
(2021), in particular the description of Algorithms 1 and 2 in the Appendix.
}

\value{The output of \code{size} is the following:

\item{starting.parameters}{The rows of this matrix are the
initial values (diagonals of covariance matrices) that were used in Stage 1 of the algorithm, and which were 
chosen from the pool of initial values in Stage 0.}

\item{starting.rejection.probs}{The null-rejection
probabilities corresponding to the initial values used in Stage 1.}

\item{first.stage.parameters}{The rows of this matrix are the 
parameters (diagonals of covariance matrices) that were obtained in Stage 1 of the algorithm.}

\item{first.stage.rejection.probs}{The
null-rejection
probabilities corresponding to the \code{first.stage.parameters}.} 

\item{second.stage.parameters}{The rows of this matrix are the 
parameters (diagonals of covariance matrices) that were obtained in Stage 2 of the algorithm.}

\item{second.stage.rejection.probs}{The null-rejection probabilities
corresponding to the \code{second.stage.parameters}.}

\item{convergence}{Convergence codes returned from \code{\link{constrOptim}}
in Stage 2 of the algorithm for each initial value.}

\item{size}{The size computed by the algorithm, i.e., the maximum of the
\code{second.stage.rejection.probs}.}
}

\examples{

#size of the classical (uncorrected) F-test in a location model
#with conventional t-critical value (5% level)

#it is known that (in this very special case) the conventional critical value 
#is size-controlling (i.e., the resulting size should be 5% (approximately))

C <- qt(.975, df = 9)^2
R <- matrix(1, nrow = 1)
X <- matrix(rep(1, length = 10), nrow = 10, ncol = 1)
hcmethod <- -1
restr.cov <- FALSE
Mp <- 100
M1 <- 5
M2 <- 1

#here, the parameters are chosen such that the run-time is low
#to guarantee a high accuracy level in the computation, 
#Mp, M1 and M2 should be chosen much higher

size(C, R, X, hcmethod, restr.cov, Mp, M1, M2)
}

\references{
\enc{Pötscher}{Poetscher}, B. M. and Preinerstorfer, D. (2021). Valid Heteroskedasticity Robust Testing. <arXiv:2104.12597>
}

\seealso{
\code{\link{davies}}, \code{\link{constrOptim}}.
}                  