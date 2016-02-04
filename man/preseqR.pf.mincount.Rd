\name{preseqR.pf.mincount}
\alias{preseqR.pf.mincount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the expected number of species represented r or more times}
\description{
The function estimates the expected number of species represented at least 
r times in a random sample based on the initial sample.
}
\usage{
preseqR.pf.mincount(n, mt = 100, ss = NULL, max.extrapolation = NULL, r=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column 
    is \eqn{n_j}, the number of species with each species represented by \eqn{j} 
    individuals in the initial sample. The first column must be sorted in an
    ascending order.
}
  \item{mt}{
    An positive integer equal to the maximum degree allowed in the continued
    fraction approximation. Default is 100.
}
  \item{ss}{
    An positive double equal to the step size between samples. Default value
    is the size of the initial sample.
}
  \item{max.extrapolation}{
    A positive double equal to the maximum possible size of a random sample. 
    Default value is 100 times the size of the initial sample.
}
  \item{r}{
    A vector of positive integers. Default is 1.
  }
}
\details{
The difference between this function and 
\code{\link{preseqR.pf.mincount.bootstrap}} is that no bootstrapping for the 
initial sample. Therefore the function can not construct confidence intervals.
However, the function is much faster. In general, we recommand
\code{\link{preseqR.pf.mincount.bootstrap}} for estimating the expected number
of species represented at least r times in a sample. 

See \code{\link{preseqR.pf.mincount.bootstrap}} for more.
}
\value{
  \item{PF.elements}{
    A list of two components for the estimator. The estimator can be
    expressed as \deqn{\hat{E}(S_r(t)) = \sum_{i=1}^M c_i \left(\frac{t}{t - x_i}\right)^r.}
    PF.elements contains both coefficients \eqn{c_i} and roots \eqn{x_i}.
  }
  \item{yield.estimates}{
    A list of two-column matrices representing estimates of the expected number
    of species represented r or more times in a random sample. Each matrix 
    in the list corresponds to a specified value of r. The first column of a matrix is 
    the size of a random sample; the second column is the estimate of the
    number of species represented r or more times in the sample. 
  }
  \item{NULL}{
    NULL if the method fails to construct an estimator.
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
# library(preseqR)

## import data
# data(ShakespeareWordHist)

## construct the estimator and predicting estimates of the number of unique word
## represented at least once, twice or twenty times in larger samples
# result = preseqR.pf.mincount(ShakespeareWordHist, r=c(1,2,20))

## print the elements of the estimator
# result$PF.elements

## print the number of word types observed at least once as a function of words
# result$yield.estimates[[1]]
## print the number of word types observed at least twice
# result$yield.estimates[[2]]
## print the number of word types observed at least twenty times
# result$yield.estimates[[3]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sufficient representation, Partial fraction }
