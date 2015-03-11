\name{preseqR.interpolate.distinct}
\alias{preseqR.interpolate.distinct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Interpolating a species accumulation curve
}
\description{
    Interpolating a species accumulation curve given a survey
}
\usage{
  preseqR.interpolate.distinct(ss, hist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{
    An positive double equal to the step size between samples.
}
  \item{hist}{
    Frequencies of the number of individuals of each species captured. 
    The data must be a two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with \eqn{j} individuals observed in the
    sample. The first column must be sorted in an ascending order.
}
}
\details{
  Assume that the sub sample follows a hypergeometric distribution for a
  given survey. So the expected number of unique species for a fixed size of 
  the sub sample is calculated based on the formula in
  Heck Jr, K. L., van Belle, G., & Simberloff, D. (1975).

  ss should be at least one. 
}
\value{
    A two-column matrix representing a interpolated accumulation curve.
    The first column is the size of the survey; the second column is the 
    expected number of unique species being observed.

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

## interpolating an accumulation curve; the step size is 1e5
preseqR.interpolate.distinct(1e5, ShakespeareWordHist)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Interpolation }
