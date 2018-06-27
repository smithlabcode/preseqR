\name{preseqR.rSAC}
\alias{preseqR.rSAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Best practice for \eqn{r}-SAC -- a fast version
}
\description{
\code{preseqR.rSAC} predicts the expected number of species represented at least 
\eqn{r} times in a random sample based on the initial sample.
}
\usage{
preseqR.rSAC(n, r=1, mt=20, size=SIZE.INIT, mu=MU.INIT)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{N_j}, the number of species with each species represented exactly \eqn{j}
    times in the initial sample. The first column must be sorted in an
    ascending order.
  }
  \item{mt}{
    A positive integer constraining possible rational function
    approximations. Default is 20.
  }
  \item{r}{
    A positive integer. Default is 1.
  }
  \item{size}{
    A positive double, the initial value of the parameter \code{size} in
    the negative binomial distribution for the EM algorithm. Default value is 1.
  }
  \item{mu}{
    A positive double, the initial value of the parameter \code{mu} in the
    negative binomial distribution for the EM algorithm. Default value is 0.5.
  }
}
\details{
  \code{preseqR.rSAC} combines the nonparametric approach using the rational
  function approximation and the parametric approach using the
  zero-truncated negative binomial (ZTNB). For a given initial sample, if the sample
  is from a heterogeneous population, the function calls
  \code{\link{ds.rSAC}}; otherwise it calls \code{\link{ztnb.rSAC}}. The degree
  of heterogeneity is measured by the coefficient of variation, which is
  estimated by the ZTNB approach.

  \code{preseqR.rSAC} is the fast version of \code{\link{preseqR.rSAC.bootstrap}}.
  The function does not provide the confidence interval. To obtain the
  confidence interval along with the estimates, one should use the function
  \code{\link{preseqR.rSAC.bootstrap}}.
}
\value{
  The estimator for the \eqn{r}-SAC. The input of the estimator is a vector of
  sampling efforts \eqn{t}, i.e., the relative sample sizes comparing with the initial
  sample. For example, \eqn{t = 2} means a random sample that is twice the size of
  the initial sample.
}
\references{
Deng, C., Daley, T., Calabrese, P., Ren, J., & Smith, A.D. (2016). Estimating
the number of species to attain sufficient representation in a random sample.
arXiv preprint arXiv:1607.02804v3.
}
\author{
  Chao Deng
}

\examples{
## load library
library(preseqR)

## import data
data(FisherButterfly)

## construct the estimator for SAC
estimator1 <- preseqR.rSAC(FisherButterfly, r=1)
## The number of species represented at least once in a sample, 
## when the sample size is 10 or 20 times of the initial sample
estimator1(c(10, 20))

## construct the estimator for r-SAC
estimator2 <- preseqR.rSAC(FisherButterfly, r=2)
## The number of species represented at least twice in a sample, 
## when the sample size is 50 or 100 times of the initial sample
estimator2(c(50, 100))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ estimator, r-SAC, RFA, ZTNB, best practice}
