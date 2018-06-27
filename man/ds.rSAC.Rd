\name{ds.rSAC}
\alias{ds.rSAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  RFA estimator
}
\description{
\code{ds.rSAC} predicts the expected number of species represented at least 
\eqn{r} times in a random sample, based on the initial sample. 
}
\usage{
ds.rSAC(n, r=1, mt=20)
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
    An positive integer constraining possible rational function
    approximations. Default is 20.
  }
  \item{r}{
    A positive integer. Default is 1.
  }
}
\details{
  The estimator is based on an empirical Bayes approach using rational
  function approximation (RFA), as described in the paper in the 
  references section.

  \code{ds.rSAC} is the fast version of \code{\link{ds.rSAC.bootstrap}}.
  The function does not provide the confidence interval. To obtain the
  confidence interval along with the estimates, one should use the function
  \code{\link{ds.rSAC.bootstrap}}.
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
ds1 <- ds.rSAC(FisherButterfly, r=1)
## The number of species represented at least once in a sample, 
## when the sample size is 10 or 20 times of the initial sample
ds1(c(10, 20))

## construct the estimator for r-SAC
ds2 <- ds.rSAC(FisherButterfly, r=2)
## The number of species represented at least twice in a sample, 
## when the sample size is 50 or 100 times of the initial sample
ds2(c(50, 100))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ estimator, r-SAC, RFA, nonparametric}
