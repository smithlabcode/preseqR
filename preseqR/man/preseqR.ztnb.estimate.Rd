\name{preseqR.ztnb.estimate}
\alias{preseqR.ztnb.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to predict the number of distinct items
}
\description{
The function estimates the number of distinct items if an experiment with
different sample size conducts. It assumes that the count of each item follows a
negative binomial distribution. The number of distinct items is estimated by its
expectation.
}
\usage{
preseqR.ztnb.estimate(hist, n, header = FALSE)
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
  \item{n}{
	  A positive real number represents the sample size of an experiment.
}
  \item{header}{
	  A logic number represents whether or not histogram file contains a header.
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
	The expected number of distinct items in the new experiment with sample size
	n.
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

## estimate the number of distinct butterflies if 10000 butterflies are trapped
preseqR.ztnb.estimate(hist.count, 1e4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ estimate }
\keyword{ negative }
\keyword{ binomial }
\keyword{ EM }
