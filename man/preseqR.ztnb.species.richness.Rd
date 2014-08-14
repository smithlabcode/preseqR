\name{preseqR.ztnb.species.richness}
\alias{preseqR.ztnb.species.richness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict the number of species through a ZTNB
}
\description{
    The function extrapolates the number of species given samples with different
    sizes. It assumes that the observed count of species follows a 
    zero-truncated negative binomial (ZTNB) distribution.
}
\usage{
preseqR.ztnb.species.richness(hist, ss = NULL, max.extrapolation = NULL, header = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
    Frequencies of count data in an inital survey. \code{hist} can be either a 
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
    survey. Default value is the \code{ss} times 100.
}
}
\details{
    This is multiple samples version of \code{preseqR.ztnb.estimate}. Instead of
    extrapolating the number of species given the size of a sample, the function
    extrapolate number of species given various sizes of samples.
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
  \code{\link{preseqR.ztnb.estimate}}
}

\examples{
## load library
library(preseqR)

## import data
data(FisherButterflyCaptures)

## extrapolation
preseqR.ztnb.species.richness(FisherButterflyCaptures)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ZTNB }
