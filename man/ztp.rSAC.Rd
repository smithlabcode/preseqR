\name{ztp.rSAC}
\alias{ztp.rSAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ZTP estimator
}
\description{
  \code{ztp.rSAC} predicts the expected number of species represented at least 
\eqn{r} times in a random sample, based on the initial sample. 
}
\usage{
ztp.rSAC(n, r=1)
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
  \item{r}{
    A positive integer. Default is 1.
  }
}
\details{
  The statistical assumption is that for each species the number of individuals
  in a sample follows a Poisson distribution. The Poisson rate \code{lambda}
  is the same among all species. So the random variable \code{X}, which is 
  the number of species represented x (x > 0) times, follows a zero-truncated 
  Poisson distribution. The unknown parameters are estimated by 
  Cohen (1960). Based on the estimated distribution,
  we calculate the expected number of species in a random sample. 
}
\value{
  The estimator for the \eqn{r}-SAC. The input of the estimator is a vector of
  sampling efforts \eqn{t}, i.e., the relative sample sizes comparing with the initial
  sample. For example, \eqn{t = 2} means a random sample that is twice the size of
  the initial sample.
}
\author{
  Chao Deng
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{
Cohen, A. Clifford. "Estimating the parameter in a conditional Poisson 
distribution." Biometrics 16, no. 2 (1960): 203-211.
}

\examples{
## load library
library(preseqR)

## import data
data(FisherButterfly)

## construct the estimator for SAC
ztp1 <- ztp.rSAC(FisherButterfly, r=1)
## The number of species represented at least once in a sample, 
## when the sample size is 10 or 20 times of the initial sample
ztp1(c(10, 20))

## construct the estimator for r-SAC
ztp2 <- ztp.rSAC(FisherButterfly, r=2)
## The number of species represented at least once in a sample, 
## when the sample size is 10 or 20 times of the initial sample
ztp2(c(50, 100))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ estimator, r-SAC, parametric, Poisson}
