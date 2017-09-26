\name{ds.rSAC.bootstrap}
\alias{ds.rSAC.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the number of species represented r or more times
}
\description{
The function estimates the expected number of species represented at least r
times in a random sample based on the initial sample. The initial sample is
then bootstrapped to improve the stability of estimates.
}
\usage{
preseqR.rSAC.bootstrap(n, r=1, mt=20, times=100)
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
    An positive integer to truncate the power series of the average discovery
    rate. Default is 20.
  }
  \item{r}{
    A positive integer. Default is 1.
  }
  \item{times}{
    The number of bootstrap samples.
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

## construct the estimator for the number of species represented at least once
## in a random sample
# estimator1 = preseqR.rSAC.bootstrap(FisherButterflyHist, r=1)
## The number of species represented at least once, when the sample size is 
## 10 or 20 times of the initial sample
# estimator1$FUN.nobootstrap$FUN(c(10, 20))

## construct the estimator for the number of species represented at least twice
## in a random sample
#estimator2 <- preseqR.rSAC.bootstrap(FisherButterflyHist, r=2)
## The number of species represented at least twice, when the sample size is 
## 50 or 100 times of the initial sample
# estimator2(c(50, 100))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Estimator, r-SAC, Bootstrap, Nonparametric}
