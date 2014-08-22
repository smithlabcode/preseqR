\name{preseqR.rfa.estimate}
\alias{preseqR.rfa.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the value of the continued fraction
}
\description{
This function calculates the value of a continued fraction at the
value t which will correspond to the relative number of additonal
individuals to be captured.
}
\usage{
preseqR.rfa.estimate(CF, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CF}{
    A continued fraction.
  }
  \item{t}{
    A double as the argument.
  }
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

## calculate values when the argument equals to one
preseqR.rfa.estimate(CF, 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rational Function Approximation }
