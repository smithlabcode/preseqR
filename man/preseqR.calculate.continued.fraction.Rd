\name{preseqR.calculate.continued.fraction}
\alias{preseqR.calculate.continued.fraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the value of a continued fraction
}
\description{
The function calculates the value of a continued fraction given its coordinates.
}
\usage{
preseqR.calculate.continued.fraction(CF, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CF}{
    A continued fraction.
  }
  \item{x}{
    A real number as the coordinate.
  }
}

\value{
    A real number is the function value.
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

## estimate the continued fraction
result = preseqR.continued.fraction.estimate(ShakespeareWordCounts)
## record the constructed continued fraction
CF = result$continued.fraction
## calculate values when the coordinate equals to one
preseqR.calculate.continued.fraction(CF, 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ calculate }
\keyword{ continued }
\keyword{ fraction }
