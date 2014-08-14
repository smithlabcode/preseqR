\name{preseqR.print2file}
\alias{preseqR.print2file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write results down to a file
}
\description{
The function writes down a continued fraction approximation or a matrix into a
file
} 
\usage{
preseqR.print2file(X, filename = NULL, digit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    A variable representing either a continued fraction with "RFA" attribute or 
	a matrix with "matrix" attribute. See details for more information.
}
  \item{filename}{
    A string that is the name of the output file.
}
  \item{digit}{
    A nonnegative number controlling the format of output file. It determines
    the number of digits preserved after the decimal point. Default value is 4
    for a continued fraction class and 2 for a matrix class.
}
}

\value{
  It returns either 0 or 1. 0 means normal. 1 means abnormal. If the returning
  value is 1, an error message will print on stderr.
}

\author{
  Chao Deng
}

\note{
  The function uses \code{formatC} in R to control the output format. If users
  do not concern the format, you can use \code{write.table} to write down a
  matrix into a file.
}

\examples{
## load library
library(preseqR)

## import data
data(ShakespeareWordCounts)

## construct a continued fraction approximation and estimate the number of
## observed word types
## if different quantities of Shakespeare's works are discovered.
result = preseqR.rational.function.approximation.estimate(ShakespeareWordCounts)
preseqR.print2file(result$continued.fraction, "continued_fraction_approximation.txt")
preseqR.print2file(result$yield.estimates, "observed_wordtypes.txt")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Print }
