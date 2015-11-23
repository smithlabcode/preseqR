\name{preseqR.pf.mincount}
\alias{preseqR.pf.mincount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the number of species observed r or more times as a function of captures
}
\description{
This function estimates the expected number of species observed at least r times 
as a function of the number of total individuals captured.}
\usage{
preseqR.pf.mincount(n, mt = 100, ss = NULL, max.extrapolation = NULL, r=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    Frequencies of the number of individuals of each species captured. 
    The data must be a two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with \eqn{j} individuals observed in the
    sample. The first column must be sorted in an ascending order.
}
  \item{mt}{
    An positive integer equal to the maximum degree allowed in the continued
    fraction approximation. Default is 100.
}
  \item{ss}{
    An positive double equal to the step size between samples. Default value
    is the size of the sample in the initial survey.
}
  \item{max.extrapolation}{
    A positive double equal to the maximum possible size of a sample in a
    survey. Default value is 100 times the size of the initial
    survey.
}
  \item{r}{
    A vector of positive integers equal to minimum counts.
  }
}
\details{
    This function estimates the expected number of species observed at least r 
    times through rational function approximation without bootstrapping. 
}
\value{
  \item{PF.elements}{
    A list contains two components of the estimator for predicting the number of
    species observed at least r times:
    ps.coeffs is the coeffiencients of the Good and Toulmin power series; 
    coef is the coeffiencients of the estimator;
    denom.roots is the roots of the estimator. 
  }
  \item{yield.estimates}{
    A two-column table representing estimated number of species given samples of
    different sizes.
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

## construct the estimator and predicting estimates of the number of word types
## observed at least once and twice
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
