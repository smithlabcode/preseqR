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
preseqR.ztnb.complexity.curve <- function(hist, ss = NULL, 
max.extrapolation = NULL, header = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
	  The histogram of an experiment. It can be either the file name of a 
	  histogram or its count vector. For a histogram file, it should have two 
	  columns. Each row (x, y) represents there are y distinct items, each of
	  which occurs x times in the experiment. Both x, y are positive integers
	  and The first columns must be an increasing order. 

	  A count vector is a non-negative integral vector. The ith coordinate x in
	  the vector represents that x distinct items appear in the experiment, each
	  of which occurs i times.
}
  \item{header}{
	  An integer represents diagonal value for a constructed continued fraction.
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
	The function assumes that the count of each item follows a negative binomial
	distribution. Given the histogram of the initial experiment, an EM algorithm
	is used to learn the parameters (size = s, mu = m) in the negative binomial 
	distribution. According to theory, for experiments with different sample 
	sizes, they also follow a negative binomial distribution with parameters 
	(size = s, mu = m * t). t is the proportion of sample size to the size of
	the inital experiment. 
	
	Given the sample size of an experiment, preseqR.ztnb.complexity.curve will
	calculate the parameters (size, mu) of a negaative binomial distribution
	corresponding to the experiment. Then the expected number of distinct items
   	is used as estimation for the experiment.
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
data(FisherButterlyCaptures)

## convert a histogram into a count vector
hist.count = vector(mode = 'numeric', length = max(FisherButterlyCaptures$n))
hist.count[FisherButterlyCaptures$n] = FisherButterlyCaptures$S

## estimate a complexity curve
preseqR.ztnb.complexity.curve(hist.count)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ complexity }
\keyword{ negative }
\keyword{ binomial }
\keyword{ EM }
