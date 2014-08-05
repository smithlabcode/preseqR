\name{preseqR.interpolate.distinct}
\alias{preseqR.interpolate.distinct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to estimate the number of distinct items
}
\description{
The function estimates the number of distinct items when experiments with small
sample sizes conduct.
}
\usage{
preseqR.interpolate.distinct(hist.count, ss)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist.count}{
	  The count vector of a histogram. A count vector is a non-negative integral
	  vector. The ith coordinate x in the vector represents that x distinct 
	  items appear in the experiment, each of which occurs i times.
}
  \item{ss}{
	  A positive real number represents the step size. It defines the sample 
	  size of an additional experiment.
}
}

\details{
	The function is a part of preseqR.continued.fraction.estimate used to 
	estimate the number of distinct items when the sample size of an experiment
	is less than the initial experiment. Sampling without replacement is used
	to interpolate.
}
\value{
    A two-columns table represents estimated number of distinct items given 
  	experiments with various sample sizes.
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

## calculate the total number of reads
total.reads = 0
for (i in 1:length(hist.count)) total.reads <- total.reads + i * hist.count[i]

## interpolate with step size around fifth of total.reads
preseqR.interpolate.distinct(hist.count, ss = total.reads / 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ interpolate }
\keyword{ nonreplacemnt }
\keyword{ sampling }
