\name{preseqR.ztnb.species.accum.curve}
\alias{preseqR.ztnb.species.accum.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the number of species through a ZTNB
}
\description{
  The function estimates the expected number of species represented in a random 
  sample by fitting a negative binomial distribution to the initial sample.
}
\usage{
preseqR.ztnb.species.accum.curve(n, ss = NULL, max.extrapolation = NULL,
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
}
\details{
    The statistical assumption is that the number of individuals in a species
    follows a negative binomial distribution. The parameters are estimated by
    the function \code{preseqR.ztnb.em}. Based on the estimated distribution,
    we calculate the expected number of species in a random sample. Details of
    the method see supplementary of Daley, T., & Smith, A. D. (2013).
}
\value{
    A two-column matrix representing for estimates of the expected number of        
    species represented in a random sample. The first column is the size of the  
    random sample; the second column is estimates of the expected number of species
    represented in the sample. 
}
\author{
  Chao Deng
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{preseqR.ztnb.em}}
}

\examples{
## load library
library(preseqR)

## import data
data(FisherButterflyHist)

## extrapolation
preseqR.ztnb.species.accum.curve(FisherButterflyHist)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Zero truncated negative binomial }
