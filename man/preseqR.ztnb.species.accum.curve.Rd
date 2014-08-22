\name{preseqR.ztnb.species.accum.curve}
\alias{preseqR.ztnb.species.accum.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict the number of species through a ZTNB
}
\description{
    This function extrapolates the number of species as a function of the number
    of captures. It assumes that the number of observations in an individual
    species follows a negative binomial distribution.
}
\usage{
preseqR.ztnb.species.accum.curve(hist, ss = NULL,
max.extrapolation = NULL, header = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
    Frequencies of the number of individuals of each species captured. 
    It can be either a file name of the data or a variable of a matrix 
    containing the data. The data must be a two-column table.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with \eqn{j} individuals observed in the
    sample. The first column must be sorted in an ascending order.
}
  \item{header}{
    A boolean variable number representing whether or not the frequencies of 
    count data file contains a header.
}
  \item{ss}{
    An positive double equal to the step size between samples. Default value
    is the size of the sample in the initial survey.
}
  \item{max.extrapolation}{
    A positive double equal to the maximum possible size of a sample in a
    survey. Default value is the 100 times the size of a sample in the initial 
    survey.
}
}
\details{
    The function applies \code{preseqR.ztnb.em} for parameter estimation.
    Then it extrapolates the number of species as a function of the number of
    captures through estimated negative binomial distribution. See supplementary
    of Daley, T., & Smith, A. D. (2013)
}
\value{
    A two-columns table represents estimated number of species given samples
    with different sizes.
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
\keyword{ ZTNB }
