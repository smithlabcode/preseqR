\name{preseqR.sample.cov}
\alias{preseqR.sample.cov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the probability of observing a species represented at least r times
in a random sample
}
\description{
The function estimates the probability of observing a species represented at
least r times in a random sample
}
\usage{
  preseqR.sample.cov(n, r=1, mt=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with each species represented \eqn{j}
    times in the initial sample. The first column must be sorted in an
    ascending order.
  }
  \item{r}{
    A positive integer. Default is 1.
  }
  \item{mt}{
    An positive integer constraining possible rational function
    approximations. Default is 20.
  }
}
\details{
  Suppose a sample is given and one more individual is randomly drawn from the
  population. The function estimates the probability of the species, which 
  represents the individual, has been observed at least \eqn{r} times in the
  sample. When \eqn{r = 1}, the probability is called the sample coverage.

  Let \eqn{N_j} be the number of species represented exactly \eqn{j} times in 
  a sample. The probability of observing a species represented at
  least \eqn{r} times in the sample is estimated as 
  \eqn{\sum_{j=r+1}^\infty jN_j / \sum_{j=1}^\infty jN_j}. The theory is
  described in Mao and Lindsay (2002). For a random sample
  where \eqn{N_j} is unknown, a modified rational function approximation is
  first used to predict the value of \eqn{N_j}. Then the estimates are
  substituted to obtain an estimator for the probability of observing a species
  represented at least \eqn{r} times in the sample.

  This function is the fast version of \code{\link{preseqR.sample.cov.bootstrap}}.
  The function does not provide an estimate for the confidence interval.
  If one needs the confidence interval, please use 
  \code{\link{preseqR.sample.cov.bootstrap}} instead.
}
\value{
  The estimator for the probability of observing a species represented at least
  r times in a sample as a function of the sample size.
  The input of the estimator is a vector of sampling 
  efforts t, i.e. the relative sample sizes comparing with the initial sample.
  For example, t = 2 means a random sample that is twice the size of the 
  initial sample.
}
\references{
Good, I. J. (1953). The population frequencies of species and the estimation of
population parameters. Biometrika, 40(3-4), 237-264.

Mao, C. X. and Lindsay, B. G. (2002). A Poisson model for the coverage problem
with a genomic application. Biometrika, 89(3), 669-682.

Deng, C and Smith, A. D. (2016). Estimating the number of species to attain 
sufficient representation in a random sample. arXiv preprint arXiv:1607.02804
}
\author{
  Chao Deng
}

\examples{
## load library
library(preseqR)

## import data
data(FisherButterflyHist)

## construct the estimator for the sample coverage
## in a random sample
estimator1 <- preseqR.sample.cov(FisherButterflyHist, r=1) 
## Given a sample that is 10 times or 20 times the size of an initial samples,
## suppose one randomly draws one more individual from the population. The
## value of the function is the probability that the representing species 
## has been observed in the sample
estimator1(c(10, 20))

## construct the estimator for the probability of observing a species represented
## at least r times in a random sample for r = 2
estimator2 <- preseqR.sample.cov(FisherButterflyHist, r=2)
## the probability when the sample size is 50 times or 100 times of the initial
## sample
estimator2(c(50, 100))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Estimator, RFA, Sample coverage}
