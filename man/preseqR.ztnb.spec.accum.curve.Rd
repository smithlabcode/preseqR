\name{preseqR.ztnb.spec.accum.curve}
\alias{preseqR.ztnb.spec.accum.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict the number of species through a ZTNB
}
\description{
    The function extrapolates the number of species as a function of the number
    of captures. It assumes that the number of observations in an individual
    species follows a zero-truncated negative binomial (ZTNB) distribution.
}
\usage{
preseqR.ztnb.spec.accum.curve(hist, ss = NULL,
max.extrapolation = NULL, header = FALSE)
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
  \item{header}{
    A boolean variable number representing whether or not the frequencies of 
    count data file contains a header.
}
  \item{ss}{
    An positive double representing the step size between samples. Default value
    is the size of the sample in the initial survey.
}
  \item{max.extrapolation}{
    A positive double representing the maximum possible size of a sample in a
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
data(FisherButterflyCaptures)

## extrapolation
preseqR.ztnb.spec.accum.curve(FisherButterflyCaptures)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ZTNB }
