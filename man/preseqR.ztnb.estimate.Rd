\name{preseqR.ztnb.estimate}
\alias{preseqR.ztnb.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the number of new species through a ZTNB
}
\description{
This function estimates the expected number of new species in a random sample
that is t times size of the initial sample based on a zero truncated negative
binomial.
}
\usage{
preseqR.ztnb.estimate(n, t, size=SIZE.INIT, mu=MU.INIT)
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
  \item{t}{
    t is the relative size of a random sample to the initial sample. 
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
    A zero-truncated negative binomial distribution is fitted to the initial
    sample. The function \cr \code{preseqR.ztnb.em} is used to estimate the parameters. 
    Based on the estimated distribution, we calculate the expected gain in new
    species in a random sample. Details of the method see supplementary
    of Daley, T., & Smith, A. D. (2013).
}

\value{
    a double as the function value. 
}
\author{
    Chao Deng
}

\seealso{
    \code{\link{preseqR.ztnb.em}}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
## load library
library(preseqR)

## import data
data(FisherButterflyHist)

## estimate the number of new butterfly species when the size of 
## a random sample is twice the size of the initial sample
preseqR.ztnb.estimate(FisherButterflyHist, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Zero truncated negative binomial }
