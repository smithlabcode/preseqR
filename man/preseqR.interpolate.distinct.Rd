\name{preseqR.rfa.species.accum.curve}
\alias{preseqR.rfa.species.accum.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Interpolating a species accumulation curve
}
\description{
    interpolating a species accumulation curve given a survey
}
\usage{
  preseqR.interpolate.distinct <- function(hist.count, ss)
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
  \item{ss}{
    An positive double equal to the step size between samples. Default value
    is the size of the sample in the initial survey.
}
}
\details{

}
\value{
    A two-column table representing a interpolated accumulation curve.
    The first column is the size of the survey; the second column is the expected number of unique species.

    NULL if failed.
}
\references{
  Heck 1975
}
\author{
  Chao Deng
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## load library
library(preseqR)

## import data
data(ShakespeareWordHist)

## interpolating an accumulation curve; the step size is 1e5
preseqR.rfa.species.accum.curve(ShakespeareWordHist, 1e5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rational Function Approximation }
