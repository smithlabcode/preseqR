\name{preseqR.extrapolate.distinct}
\alias{preseqR.extrapolate.distinct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to predict the number of distinct items
}
\description{
preseqR.extrapolate.distinct predicts the number of distinct items when 
experiments with large sample sizes conduct.
}
\usage{
preseqR.extrapolate.distinct(hist.count, CF, start.size = NULL, 
step.size = NULL, max.size = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist.count}{
	  The count vector of a histogram. A count vector is a non-negative integral
	  vector. The ith coordinate x in the vector represents that x distinct 
	  items appear in the experiment, each of which occurs i times.
  }
  \item{CF}{
	  A continued fraction with "continued.fraction" attribute.
  }
  \item{start.size}{
	  A positive real number represent a starting size of experiment for 
	  extrapolation.
  }
  \item{step.size}{
	  A positive real number represents the step size. It defines the size of 
	  an additional experiment.
  }
  \item{max.size}{
	  A positive real number represents the upperbound for extrapolation. In
	  other word, it defines the maximum possible size of experiment, in which 
	  the constructed continued fraction would be called to estimate the total
	  number of distinct items.   
  }
}
\details{
	The function is a part of preseqR.continued.fraction.estimate used to 
	estimate the number of distinct items when the sample size of an experiment
	is beyond the initial experiment. The constructed continued fraction is used
   	to extrapolate.
}
\value{
  \item{yield.estimates}{
	  A two-columns table represents estimated number of distinct items given 
	  various sizes of experiments.
  }
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

## construct a continued fraction
CF = preseqR.continued.fraction.estimate(hist.count)$continued.fraction

## calculate the total number of reads
total.reads = 0
for (i in 1:length(hist.count)) total.reads <- total.reads + i * hist.count[i]

## extrapolate
step.size = total.reads;
preseqR.extrapolate.distinct(hist.count, CF, start.size = total.reads + step.size, 
step.size = step.size, max.size = 100 * step.size)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ extrapolate }
\keyword{ continued }
\keyword{ fraction }
