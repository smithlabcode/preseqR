\name{preseqR.ztnb.em}
\alias{preseqR.ztnb.em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting a zero-truncated negative binomial distribution
}
\description{
This function fits a zero-truncated negative binomial (ZTNB) distribution
to the initial sample.
Since the species with zero observations are missed in the sample, an
EM algorithm is used to estimate the parameters assuming the number of
individuals for each species follows a Negative Binomial distribution
with the zero counts as a missing latent data.
}
\usage{
preseqR.ztnb.em(n, size = SIZE.INIT, mu = MU.INIT)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with each species represented \eqn{j}
    times in the initial sample. The first column must be sorted in an
    ascending order.
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
    See the supplement of Daley and Smith (2013).
}

\value{
  \item{size}{
    The estimate of the parameter \code{size} in the negative binomial.
}
  \item{mu}{
    The estimate of the parameter \code{mu} in the negative binomial.
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
\keyword{ Zero truncated negative binomial, EM}
