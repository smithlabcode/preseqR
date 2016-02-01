\name{preseqR.nonreplace.sampling}
\alias{preseqR.nonreplace.sampling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Sampling without replacement
}
\description{
    Generating a histogram by subsampling without replacement a given sample.
}
\usage{
  preseqR.nonreplace.sampling(size, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
    An positive integer representing the size of the subsample.
  }
  \item{n}{
    A two-column matrix.
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with each species represented by \eqn{j}
    individuals in the initial sample. The first column must be sorted in an
    ascending order.
  }
}
\details{
  The function sample() in R is used to implement the function. We
  wrap the sample() function in a way that both input and output are histograms. 
}
\value{
    A two-column matrix as a subsample.
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with each species represented by \eqn{j}
    individuals in the subsample. 
}
\references{
\url{https://stat.ethz.ch/R-manual/R-patched/library/base/html/sample.html}
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

## generate a sub sample of size 1000.
preseqR.nonreplace.sampling(size=1000, FisherButterflyHist)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Sampling}
