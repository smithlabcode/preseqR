\name{preseqR.rfa.estimate}
\alias{preseqR.rfa.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the value of the continued fraction
}
\description{
This function calculates the value of a continued fraction at the
value t
}
\usage{
preseqR.rfa.estimate(CF, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CF}{
    A truncated continued fraction.
  }
  \item{t}{
    A double as the argument.
  }
}
\details{
  t is the ratio of the number of species observed between the new survey
and the initial survey. The function value represents the number of unique 
new species captured that are unobserved in the initial survey. The t must
be a positive real number no less than 1. 
}
\value{
    A double as the function value.
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

## the rational fraction approximation
result = preseqR.rfa.curve(ShakespeareWordHist)

## record the constructed continued fraction approximation
CF = result$continued.fraction

## calculate values when the argument equals to two
preseqR.rfa.estimate(CF, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rational Function Approximation }
