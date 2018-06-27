\name{ztnb.rSAC}
\alias{ztnb.rSAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ZTNB estimator
}
\description{
\code{ztnb.rSAC} predicts the expected number of species represented at least 
\eqn{r} times in a random sample, based on the initial sample. 
}
\usage{
ztnb.rSAC(n, r=1, size=SIZE.INIT, mu=MU.INIT)
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
  The statistical assumption is that for each species the number of individuals
  in a sample follows a Poisson distribution. The Poisson rate \code{lambda}
  are numbers generated from a gamma distribution. So the random variable
  \code{X}, which is the number of species represented x (x > 0) times in the
  sample, follows a zero-truncated negative binomial distribution. The
  unknown parameters are estimated by the function \code{preseqR.ztnb.em} based
  on the initial sample. Using the estimated distribution, we calculate the 
  expected number of species represented at least r times in a random sample.
  Details of the estimation procedure can be found in the supplement of 
  Daley T. and Smith AD. (2013).
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

\seealso{
  \code{\link{preseqR.ztnb.em}}
}

\references{
Daley, T., & Smith, A. D. (2013). Predicting the molecular complexity of
sequencing libraries. Nature methods, 10(4), 325-327.

Deng C, Daley T & Smith AD (2015). Applications of species accumulation curves
in large-scale biological data analysis. Quantitative Biology, 3(3), 135-144. 
}

\examples{
## load library
library(preseqR)

## import data
data(FisherButterfly)

## construct the estimator for SAC
ztnb1 <- ztnb.rSAC(FisherButterfly, r=1)
## The number of species represented at least once in a sample, 
## when the sample size is 10 or 20 times of the initial sample
ztnb1(c(10, 20))

## construct the estimator for r-SAC
ztnb2 <- ztnb.rSAC(FisherButterfly, r=2)
## The number of species represented at least twice in a sample, 
## when the sample size is 50 or 100 times of the initial sample
ztnb2(c(50, 100))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ estimator, r-SAC, parametric, negative binomial}
