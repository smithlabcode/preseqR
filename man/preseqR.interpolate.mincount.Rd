\name{preseqR.interpolate.mincount}
\alias{preseqR.interpolate.mincount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Interpolating the number of species represented r or more times
}
\description{
    Interpolating the expected number of species represented at least r times
    in a random sample based on an initial sample. 
}
\usage{
  preseqR.interpolate.mincount(ss, n, r=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{
    An positive double equal to the step size between samples.
}
  \item{n}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with each species represented \eqn{j}
    times in the initial sample. The first column must be sorted in an
    ascending order.
}
  \item{r}{
    A positive integer.
  }
}
\details{
  Assume that a random sample (subsample) follows a multivariant hypergeometric
  distribution given an initial sample. The expected number of unique species
  represented at least r times in the subsample is then calculated by an 
  expended version of the formula in 
  Heck Jr, KL. et al. (1975).
}
\value{
    A two-column matrix for the number of species represented at least r times
    in a random sample. The first column is the size of the random sample; 
    the second column is the expected number of species represented at least r
    times in the sample.

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

## The expected number of species represented twice or more in a random sample 
## The step size is 1e5; the initial sample is "ShakespeareWordHist"
preseqR.interpolate.mincount(n=ShakespeareWordHist, ss=1e5, r=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Interpolation, At least r times}
