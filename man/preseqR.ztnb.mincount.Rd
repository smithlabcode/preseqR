\name{preseqR.ztnb.mincount}
\alias{preseqR.ztnb.mincount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the number of species represented at least r times in a random sample
}
\description{
  The function estimates the expected number of species represented at least r
  times in a random sample by fitting a zero-truncated negative binomial (ZTNB)
  distribution to the initial sample.
}
\usage{
preseqR.ztnb.mincount(n, ss = NULL, max.extrapolation = NULL, r=1, 
                      size=SIZE.INIT, mu=MU.INIT)
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
  \item{ss}{
    An positive double equal to the step size between samples. Default value
    is the size of the initial sample.
}
  \item{max.extrapolation}{
    A positive double equal to the maximum possible size of a random sample.
    Default value is the 100 times the size of the initial sample.
}
  \item{size}{
    A positive double setting the initial value of the parameter \code{size} in
    a negative binomial distribution for the EM algorithm. Default value is 1.
}
  \item{mu}{
    A positive double setting the initial value of the parameter \code{mu} in a
    negative binomial distribution for the EM algorithm. Default value is 0.5.
}
  \item{r}{
    A vector of positive integers. Default is 1.
  }
}
\details{
  See \code{\link{preseqR.ztnb.species.accum.curve}}.
}
\value{
    A list of two-column matrices representing estimates of the expected number
    of species represented r or more times in a random sample. Each matrix 
    in the list corresponds to a specified value of r. The first column of a
    matrix is the size of a random sample; the second column is the estimate of
    the number of species represented r or more times in the sample. 
}
\author{
  Chao Deng
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## load library
library(preseqR)

## import data
data(FisherButterflyHist)

## extrapolation
result = preseqR.ztnb.mincount(FisherButterflyHist, r=c(1,2))
## estimates of the number of butterflies represented at least once in a random sample
result[[1]]
## estimates of the number of butterflies represented at least twice
result[[2]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Zero truncated negative binomial, Sufficient representation }
