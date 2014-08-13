\name{preseqR.ztnb.complexity.curve}
\alias{preseqR.ztnb.complexity.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict the number of species through a ZTNB
}
\description{
    The function estimates the number of species given samples with different 
    sizes. It assumes that the observed count of species follows a 
    zero-truncated negative binomial (ZTNB) distribution.
}
\usage{
preseqR.ztnb.complexity.curve(hist, ss = NULL, max.extrapolation = NULL, header = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
    Frequencies of count data in an inital survey. \code{hist} can be either a 
    file name of the data or a variable of a matrix containing the data. The 
    data must be a two-column table. Each row (x, y) in the table represents
    there were y species with x individuals observed in each single species. 
    Both x, y are positive integers and the first column must be sorted in an 
    ascending order.
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
    survey. Default value is the \code{ss} times 100.
}
}
\details{
  see the supplymentary of Daley, T., & Smith, A. D. (2013).
}
\value{
    A two-columns table represents estimated number of species given samples
    with different sizes.
}
\author{
  Chao Deng
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## load library
library(preseqR)

## import data
data(FisherButterflyCaptures)

## estimate a complexity curve
preseqR.ztnb.complexity.curve(FisherButterflyCaptures)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ZTNB }
