\name{preseqR.nonreplace.sampling}
\alias{preseqR.nonreplace.sampling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Sampling without replacement
}
\description{
    Generating a histogram by subsampling without replacement.
}
\usage{
  preseqR.nonreplace.sampling(n, size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species represented \eqn{j} times in the initial
    sample. The first column must be sorted in an ascending order.
  }
  \item{size}{
    An positive integer representing the size of the subsample.
  }
}
\details{
  The function sample() in R is used to implement the function. We wrap up
  the function sample() in a way such that both input and output are histograms. 
}
\value{
  A two-column matrix as a subsample.
  The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
  is \eqn{n_j}, the number of species represented \eqn{j} times in the subsample. 
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
## generate a subsample of size 1000.
preseqR.nonreplace.sampling(n=FisherButterflyHist, size=1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sampling, histogram}
