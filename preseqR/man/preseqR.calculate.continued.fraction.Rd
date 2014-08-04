\name{preseqR.calculate.continued.fraction}
\alias{preseqR.calculate.continued.fraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to calculate the function value given its coordinate
}
\description{
preseqR.calculate.continued.fraction calculates the value of a continued 
fraction given its coordinates.
}
\usage{
preseqR.calculate.continued.fraction(CF, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CF}{
	  A continued fraction with "continued.fraction" attribute.
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

## convert a histogram into a count vector 
hist.count = vector(mode="numeric", length=max(ShakespeareWordCounts$frequencies))
hist.count[ShakespeareWordCounts$frequencies] = ShakespeareWordCounts$word.types

## estimate the continued fraction
result = preseqR.continued.fraction.estimate(hist.count)
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
