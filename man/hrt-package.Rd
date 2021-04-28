\name{hrt-package}

\alias{hrt-package}

\alias{hrt}

\docType{package}

\title{
\packageTitle{hrt}
}
\encoding{latin1}

\description{
The package \pkg{hrt} provides three functions in the context of testing affine restrictions on 
the regression coefficient vector in linear models with heteroskedastic (but independent) errors.
The methods implemented in \pkg{hrt} are based on the article 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021).
The package can be used to compute various heteroskedasticity robust 
test statistics; to numerically determine size-controlling critical values 
when the error vector is heteroskedastic and Gaussian (or, more generally, elliptically symmetric);
and to compute the size of a test that is 
obtained from a heteroskedasticity robust test statistic and a user-supplied critical value.
}

\details{\pkg{hrt} provides three functions: 
\enumerate{
\item The function \code{test.stat} can be used to evaluate the
test statistics \eqn{T_{uc}}, \eqn{T_{Het}} (with HC0-HC4 weights), 
\eqn{\tilde{T}_{uc}}, or \eqn{\tilde{T}_{Het}} (with HC0R-HC4R weights), as defined in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2021).
\item The function \code{critical.value} provides an implementation of 
Algorithm 3 in \enc{Pötscher}{Poetscher} and Preinerstorfer (2021), based 
on the auxiliary algorithm \eqn{\mathsf{A}} equal to Algorithm 1 (if \eqn{q = 1}) or Algorithm 2 
(if \eqn{q > 1}) in the same reference. This function can be 
used to determine size-controlling critical values for the test statistics 
\eqn{T_{uc}}, \eqn{T_{Het}} (with HC0-HC4 weights), \eqn{\tilde{T}_{uc}}, or 
\eqn{\tilde{T}_{Het}} (with HC0R-HC4R weights), whenever such critical values 
exist (which is checked numerically when the algorithm is applied).
\item The function \code{size} provides an implementation of Algorithm 1 or 2, 
respectively, in \enc{Pötscher}{Poetscher} and Preinerstorfer 
(2021), depending on whether \eqn{q = 1} or \eqn{q > 1}. Given a user-supplied 
critical value, the respective algorithm can be used to determine 
the size of a test based on one of the test statistics \eqn{T_{uc}}, 
\eqn{T_{Het}} (with HC0-HC4 weights), \eqn{\tilde{T}_{uc}}, or 
\eqn{\tilde{T}_{Het}} (with HC0R-HC4R weights).

We refer the user to the description of the three functions below, 
and to \enc{Pötscher}{Poetscher} and Preinerstorfer (2021) for details 
concerning the framework, the test statistics, the 
algorithms, and the underlying theoretical results.
}
}
\references{
\enc{Pötscher}{Poetscher}, B. M. and Preinerstorfer, D. (2021). Valid Heteroskedasticity Robust Testing. <arXiv:2104.12597>
}