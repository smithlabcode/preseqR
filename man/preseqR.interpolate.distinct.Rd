\name{preseqR.interpolate.distinct}
\alias{preseqR.interpolate.distinct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Interpolating a species accumulation curve
}
\description{
  Interpolating the expected number of species represented in a random sample
  based on an initial sample. 
}
\usage{
  preseqR.interpolate.distinct(ss, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{
    An positive double equal to the step size between samples.
}
  \item{n}{
    A two-column matrix.
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column 
    is \eqn{n_j}, the number of species with each species represented by \eqn{j} 
    individuals in the initial sample. The first column must be sorted in an
    ascending order.
}
}
\details{
    Assume that a random sample (subsample) follows a multivariant hypergeometric 
    distribution given an initial sample. The expected number of unique species
    represented in the subsample is then calculated by the formula (1) in 
    Heck Jr, K. L., van Belle, G., & Simberloff, D. (1975).   
}
\value{
    A two-column matrix for the expected number of species represented in a random sample.
    The first column is the size of the random sample; the second column is the 
    expected number of species represented in the sample.  

    NULL if failed.
}
\references{
Heck Jr, K. L., van Belle, G., & Simberloff, D. (1975). 
Explicit calculation of the rarefaction diversity measurement and the
determination of sufficient sample size. Ecology, 1459-1461.
}
\author{
  Chao Deng
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## load library
library(preseqR)

## import data
data(ShakespeareWordHist)

## interpolating an accumulation curve; the step size is 100,000
preseqR.interpolate.distinct(100000, ShakespeareWordHist)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Interpolation }
