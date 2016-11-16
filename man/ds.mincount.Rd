\name{ds.mincount}
\alias{ds.mincount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the expected number of species represented r or more times
}
\description{
The function estimates the expected number of species represented at least 
r times in a random sample based on the initial sample.
}
\usage{
ds.mincount(n, r=1, mt=100)
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
    approximations. Default is 100.
}
  \item{r}{
    A vector of positive integers. Default is 1.
  }
}
\details{
The difference between this function and 
\code{\link{ds.mincount.bootstrap}} is that no bootstrapping for the 
initial sample. Therefore the function could be less stable than estimates by
bootstrap.
However, this function is much faster. In general, we recommand
\code{\link{ds.mincount.bootstrap}} for estimating the expected number
of species represented at least r times in a sample. 

See \code{\link{ds.mincount.bootstrap}} for more information.
}
\value{
  \item{FUN}{
    The constructed estimator for the number of species represneted at least r
    times in a sample. The input of the estimator is a vector of sampling efforts
    t, i.e. the relative sample sizes comparing with the initial sample. For example,
    t = 2 means the sample is twice the size of the initial sample.
  }
  \item{FUN.elements}{
    A list of two components for the estimator. The estimator can be
    expressed as \deqn{\hat{E}(S_r(t)) = \sum_{i=1}^l c_i \left(\frac{t}{t - x_i}\right)^r.}
    PF.elements contains both coefficients \eqn{c_i} and roots \eqn{x_i}.
  }
  \item{m}{
    The number of terms used when applying rational function approximation to 
    the power series of the average discovery rate.
  }
  \item{m.adjust}{
    The number of terms in the estiamtor, equal to l
  }
}
\references{
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
library(preseqR)

## import data
data(ShakespeareWordHist)

## construct the estimator for the number of unique word
## represented at least once, twice or twenty times in a sample
estimator = ds.mincount(ShakespeareWordHist, r=c(1,2,20))

## print the elements of the estimator
estimator$FUN.elements

## The number of unique words represented at least once, twice or twenty times
## when the sample size is 10 or 20 times of the initial sample
estimator$FUN(c(10, 20))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Estimator, At least r times, RFA}
