\name{preseqR.interpolate.mincount}
\alias{preseqR.interpolate.mincount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Interpolating the number of species observe at least r times
}
\description{
    Interpolating the number of species observe at least r times in a random sample
    based on a survey
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
    Frequencies of the number of individuals of each species captured. 
    The data must be a two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with \eqn{j} individuals observed in the
    sample. The first column must be sorted in an ascending order.
}
  \item{r}{
    A vector of positive integers equal to minimum counts.
  }
}
\details{
  Assume that the sub sample follows a hypergeometric distribution for a
  given survey. So the expected number of unique species observed at least r times for a fixed size of 
  the sub sample is calculated based on an extension of the formula in
  Heck Jr, K. L., van Belle, G., & Simberloff, D. (1975).

  ss should be at least one. 
}
\value{
    A list of two-column matrices representing interpolated species observed at least r times.
    In each matrix, the first column is the size of the survey; the second column is the 
    expected number of unique species being observed at least r times. Matrices correspond 
    to different minimum counts. 

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

## interpolating for the number of species observed once and twice in a sample 
## smaller than the size of the initial sample; the step size is 1e5
preseqR.interpolate.mincount(1e5, ShakespeareWordHist, r=c(1,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Interpolation }
