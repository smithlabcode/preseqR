\name{preseqR.ztnb.em}
\alias{preseqR.ztnb.em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a zero-truncated negative binomial distribution
}
\description{
The function fits a zero-truncated negative binomial (ZTNB) distribution
to the frequencies of count data. Since the zero count missed in the survey, an
EM algorithm is used to estimate the unobserved zero count and the associated
parameters of a negative binomial distribution.
}
\usage{
preseqR.ztnb.em(hist, size = SIZE.INIT, mu = MU.INIT, header = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
    Frequencies of count data in an initial survey. \code{hist} can be either a 
    file name of the data or a variable of a matrix containing the data. The 
    data must be a two-column table. Each row (j, n_j) in the table represents
    there were n_j species with j individuals observed in each single species. 
    Both j and n_j are positive integers and the first column must be sorted in
    an ascending order.
}
  \item{size}{
    A positive double setting the initial value of the parameter \code{size} in
    a negative binomial distribution. Default value is 1. See details.
}
  \item{mu}{
    A positive double setting the initial value of the parameter \code{mu} in a
    negative binomial distribution. Default value is 0.5. See details.
}
  \item{header}{
    A boolean variable number representing whether or not the frequencies of 
    count data file contains a header.
}
}
\details{
    See supplementary of Daley, T., & Smith, A. D. (2013).
}

\value{
  \item{size}{
    Estimated value of the parameter \code{size} in the negative binomial.
}
  \item{mu}{
    Estimated value of the parameter \code{mu} in the negative binomial.
}
  \item{loglik}{
    Log-likelihood under estimated ZTNB.
}
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

## print the parameters of a fitting negative binomial distribution
preseqR.ztnb.em(FisherButterflyHist)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ZTNB }
