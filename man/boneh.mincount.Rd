\name{boneh.mincount}
\alias{boneh.mincount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the expected number of species represented r or more times
}
\description{
The function estimates the expected number of species represented at least 
r times in a random sample based on the initial sample using a nonparametric
approach by Boneh et al. (1998).
}
\usage{
boneh.mincount(n, r=1)
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
    A vector of positive integers. Default is 1.
  }
}
\value{
    The constructed estimator for the number of species represneted at least r
    times in a sample. The input of the estimator is a vector of sampling efforts
    t, i.e. the relative sample sizes comparing with the initial sample. For example,
    t = 2 means the sample is twice the size of the initial sample.
}
\author{
  Chao Deng
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{
Boneh, S., Boneh, A., & Caron, R. J. (1998). Estimating the prediction function
and the number of unseen species in sampling with replacement.
Journal of the American Statistical Association, 93(441), 372-379.
}

\examples{
## load library
library(preseqR)

## import data
data(FisherButterflyHist)

## construct the estimator for the number of species
## represented at least once, twice or three times in a sample
boneh.estimator <- boneh.mincount(FisherButterflyHist, r=1:3)

## The number of species represented at least once, twice or three times
## when the sample size is 10 or 20 times of the initial sample
boneh.estimator(c(10, 20))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Estimator, At least r times, Nonparametric }
