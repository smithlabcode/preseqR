\name{preseqR.calculate.continued.fraction.approximation}
\alias{preseqR.calculate.continued.fraction.approximation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the value of a continued fraction
}
\description{
The function calculates the value of a continued fraction.
}
\usage{
preseqR.calculate.continued.fraction.approximation(CF, t)
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
data(ShakespeareWordCounts)

## the rational fraction approximation
result = preseqR.rational.function.approximation.estimate(ShakespeareWordCounts)

## record the constructed continued fraction approximation
CF = result$continued.fraction

## calculate values when the argument equals to one
preseqR.calculate.continued.fraction.approximation(CF, 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rational Function Approximation }
