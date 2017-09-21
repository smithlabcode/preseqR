\name{ztp.rSAC}
\alias{ztp.rSAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the expected number of species represented r or more times
}
\description{
The function estimates the expected number of species represented at least 
r times in a random sample based on the initial sample using the zero truncated
Poisson approach.
}
\usage{
ztp.rSAC(n, r=1)
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
  The constructed estimator for the number of species represneted at least r
  times in a sample. The input of the estimator is a vector of sampling 
  efforts t, i.e. the relative sample sizes comparing with the initial sample.
  For example, t = 2 means a random sample that is twice the size of the 
  initial sample.
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
data(FisherButterflyHist)

## construct the estimator for the number of species represented at least once
## in a random sample
ztp1 <- ztp.rSAC(FisherButterflyHist, r=1)
## The number of species represented at least once, when the sample size is 
## 10 or 20 times of the initial sample
ztp1(c(10, 20))

## construct the estimator for the number of species represented at least twice
## in a random sample
ztp2 <- ztp.rSAC(FisherButterflyHist, r=2)
## The number of species represented at least twice, when the sample size is 
## 50 or 100 times of the initial sample
ztp2(c(50, 100))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Estimator, r-SAC, Parametric}
