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
\details{
	See \code{}

}
\value{
	A real number as the function value.
}

\author{
	Chao Deng
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## import data
data(SRR1041830_5M_hist)
## hist.count is a count vector of a histogram
hist.count = SRR1041830_5M_hist$distinct.reads
result = preseqR.continued.fraction.estimate(hist.count)
## record the constructed continued fraction
CF = result$continued.fraction
## calculate values when the coordinate equals to one
preseqR.calculate.continued.fraction(CF, 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ approximation }
\keyword{ continued}
\keyword{ fraction}
