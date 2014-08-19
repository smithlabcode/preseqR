\name{preseqR.ztnb.estimate}
\alias{preseqR.ztnb.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict the number of species through a ZTNB
}
\description{
The function extrapolates the number of species given the size of a sample in a
survey. It assumes that the number of observed individual species follows a
zero-truncated negative binomial (ZTNB) distribution.
}
\usage{
preseqR.ztnb.estimate(hist, n, header = FALSE)
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
  \item{n}{
    A positive double representing the size of a sample in the survey.
}
  \item{header}{
    A boolean variable number representing whether or not the frequencies of 
    count data file contains a header.
}
}
\details{
    The frequencies of count data is fitted into a zero-truncated negative 
    binomial distribution. We use \code{preseqR.ztnb.em} to estimate the
    parameters. 
	
	For extrapolation the number of species through a negative binomial 
    distribution, See supplementary of Daley, T., & Smith, A. D. (2013).
}

\value{
    The expected number of species when a sample with \code{n} individuals
    observed in the survey.
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
data(FisherButterflyCaptures)

## estimate the number of distinct butterflies if 10000 butterflies are captured
preseqR.ztnb.estimate(FisherButterflyCaptures, 1e4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ZTNB }
