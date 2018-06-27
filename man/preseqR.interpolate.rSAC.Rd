\name{preseqR.interpolate.rSAC}
\alias{preseqR.interpolate.rSAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Interpolation
}
\description{
Interpolating the number of species represented at least \eqn{r} times
in a subsample given an initial sample
}
\usage{
preseqR.interpolate.rSAC(n, ss, r=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{N_j}, the number of species with each species represented exactly \eqn{j}
    times in the initial sample. The first column must be sorted in an
    ascending order.
  }
  \item{ss}{
    A positive double equal to the step size between subsamples.
  }
  \item{r}{
    A positive integer. Default is 1
  }
}
\details{
  The expected number of species represented at least \eqn{r}
  times in the subsample is estimated based on an expended version of the
  formula by Heck Jr, KL. et al. (1975).
}
\value{
  A two-column matrix for the number of species represented at least \eqn{r}
  times in a random sample. The first column is the size of the random sample; 
  the second column is the expected number of species represented at least \eqn{r}
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
data(Shakespeare)

## The expected number of distinct words represented twice or more in the
## subsample 
preseqR.interpolate.rSAC(n=Shakespeare, ss=1e5, r=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ interpolation, r-SAC}
