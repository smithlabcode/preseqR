\name{preseqR.ztnb.complexity.curve}
\alias{preseqR.ztnb.complexity.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	A function to estimate a complexity curve
}
\description{
	The function estimates the number of distinct items given experiments with 
	various sample sizes. It assumes that the count of each item follows a
	negative binomial distribution. The number of distinct items is estimated by
	its expectation.
}
\usage{
preseqR.ztnb.complexity.curve(hist, ss = NULL, max.extrapolation = NULL, header = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
	  The histogram of an experiment. It can be either a file name or a 
	  variable of a histogram. Histogram should have two columns. Each row (x, y)
	  represents there are y distinct items, each of which occurs x times in the
	  experiment. Both x, y are positive integers and the first column must be 
	  an increasing order. 
}
  \item{header}{
	  A logic number represents whether or not histogram file contains a header.
}
  \item{ss}{
	  A positive real number represents the step size. It defines the size of  
	  an additional experiment.
}
   \item{max.extrapolation}{
	  A positive real number represents the upperbound for extrapolation. In
	  other word, it defines the maximum possible size of experiment, in which 
	  the constructed continued fraction would be called to estimate the total
	  number of distinct items.   
}
}
\details{
	see \code{\link{preseqR.ztnb.estimate}} for methods.
}
\value{
	  A two-columns table represents estimated number of distinct items given 
	  experiments with different sample sizes.
}
\author{
	Chao Deng
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## load library
library(preseqR)

## import data
data(FisherButterflyCaptures)

## estimate a complexity curve
preseqR.ztnb.complexity.curve(FisherButterflyCaptures)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ complexity }
\keyword{ negative }
\keyword{ binomial }
\keyword{ EM }
