\name{preseqR.pf.mincount}
\alias{preseqR.pf.mincount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the number of species represented r or more times in a random sample
}
\description{
The function estimates the expected number of species represented at least 
r times in a random sample based on the initial sample. The only difference
between this function and \code{\link{preseqR.pf.mincount.bootstrap}} is that no
bootstrapping for the initial sample. However, the function is much faster. In
general, we recommand \code{\link{preseqR.pf.mincount.bootstrap}} for estimating 
the expected number of species represented r or more times in a random sample.
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
    individuals in the sample. The first column must be sorted in an ascending order.
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
    A vector of positive integers equal to minimum counts.
  }
}
\details{
  See \code{\link{preseqR.pf.mincount}} for the method.
}
\value{
  \item{PF.elements}{
    A list contains two components of the estimator for predicting the number of
    species observed at least r time. The estimator is expressed as 
    \deqn{\hat{E}(S_r(t)) = \sum_{i=1}^M c_i \left(\frac{t}{t - x_i}\right)^r}.
    PF.elements contains the coefficients \eqn{c_i} and the roots \eqn{x_i}.
  }
  \item{yield.estimates}{
    A two-column matrix representing for estimates of the expected number of 
    species represented in a random sample. The first column is the size of the
    random sample; the second column is estimates of the expected number of species
    represented in the sample. 
  }
  \item{NULL}{
    NULL if the method fails to construct an estimator.
}
}
\references{
}
\author{
  Chao Deng
}


\seealso{
  \code{\link{preseqR.pf.mincount.bootstrap}}
}

\examples{
## load library
library(preseqR)

## import data
data(ShakespeareWordHist)

## construct the estimator and predicting estimates of the number of unique word
## represented at least once and twice in larger samples
result = preseqR.pf.mincount(ShakespeareWordHist, r=c(1,2))

## print the elements of the estimator
result$PF.elements

## print the number of word types observed once and twice as a function of the number of words
## in Shakespeare's works.
result$yield.estimates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Minimum count }
