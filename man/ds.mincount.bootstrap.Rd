\name{ds.mincount.bootstrap}
\alias{ds.mincount.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the number of species represented r or more times
}
\description{
    The function estimates the expected number of species represented at least r
    times in a random sample based on the initial sample. The initial sample is
    bootstrapped to improve the stability of estimates.
}
\usage{
ds.mincount.bootstrap(n, r=1, mt=20, times=100)
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
  \item{mt}{
    An positive integer constraining possible rational function
    approximations. Default is 20.
}
  \item{r}{
    A vector of positive integers. Default is 1.
  }
  \item{times}{
    An positive integer representing the minimum required number of successful
    estimation. Default is 100. See detail below.
}
}
\details{
    Under a mixture of Poisson models, the expected number of
    species represented at least r times in a random sample 
    can be expressed as higher derivatives
    of the expected number of species represented at least once. 
    We first use rational function approximations to the 
    modified Good and Toulmin's (1956) non-parametric empirical Bayes power 
    series to estimate the average discovery rate. By differentiating the
    rational function approximation, we obtain an estimator 
    for the number of species represented at least r times in a random
    sample.
}
\value{
  \item{FUN.nobootstrap}{
    The estimator constructed based on the initial sample by the function. No
    bootstrap procesure is involved.
  }
  \item{FUN.bootstrap}{
    The bootstrap samples from an initial sample are used to construct estimators.
    The median value of these estimators are estimates of the number of species
    represented at least r times in a sample. 
  }
  \item{var}{
    The estimated variance for the estimator \code{FUN.nobootstrap} by bootstrap.
  }
}

\references{
Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.

Kalinin V (1965). Functionals related to the poisson distribution and 
statistical structure of a text. Articles on Mathematical Statistics and 
the Theory of Probability pp. 202-220.

Daley, T., & Smith, A. D. (2013). Predicting the molecular complexity of
sequencing libraries. Nature methods, 10(4), 325-327.
}

\author{
  Chao Deng
}

\examples{
## load library
#library(preseqR)

## import data
#data(FisherButterflyHist)

## estimate the number of species captured at least once, twice or 20 times 
## as a function of the number of individuals

# result = ds.mincount.bootstrap(FisherButterflyHist, r=c(1,2, 20), times=10)

## estimates of the number of unique words appeared at least once, twice or three
## times when the sample size 10 times the size of the initial sample

## estimates by the function ds.mincount
# result$FUN.nobootstrap$FUN(10)

## estimates by the bootstrapped estimator
# result$FUN.bootstrap(10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Estimator, At least r times, Bootstrap, RFA}
