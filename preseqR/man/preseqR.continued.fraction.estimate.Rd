\name{preseqR.continued.fraction.estimate}
\alias{preseqR.continued.fraction.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to predict the number of distinct items
}
\description{
The function constructs a continued fraction to estimate the number of distinct
items, which did not appear in an intial experiment, should be covered if an 
additional experiment conducts.
}
\usage{
preseqR.continued.fraction.estimate(hist, di = 0, mt = 100, ss = NULL, 
max.extrapolation = NULL, step.adjust = TRUE, header = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
	  The histogram of an experiment. It can be either the file name of a 
	  histogram or its count vector. For a histogram file, it should have two 
	  columns. Each row (x, y) represents there are y distinct items, each of
	  which appears x times in the experiment. Both x, y are positive integers
	  and The first columns must be an increasing order. 
	  
	  A count vector is a non-negative integral vector. The ith coordinate x in
	  the vector represents that there are x distinct items appear in the 
	  experiment, each of which appears i times.
}
  \item{header}{
	  A logic number represents whether or not histogram file contains a header.
}
  \item{di}{
	  An integer represents diagonal value for a constructed continued fraction.
}
  \item{mt}{
	  An integer represents maximum number of parameters in a continued fraction.
}
  \item{ss}{
	  A positive real number represents the step size. It defines the size of 
	  an additional experiment.
}
  \item{max.extrapolation}{
	  A positive real number represents the upperbound for extrapolation. In
	  other word, it defines the maximum possible sample size of experiment, in
	  which a constructed continued fraction would be called to estimate the
	  total number of distinct items.   
}
  \item{step.adjust}{
	  A logic value represents whether or not the step size could change. If
	  step.adjust equals TRUE, the value of step size (ss) could adjust in 
	  respond to the size of an initial experiment.
  }
}
\details{
Items could be anything depending on the experiment. For example, items could be
reads or molecules if it is an DNA or RNA sequencing experiment. Items could be
species if it is a population experiment. Generally, the function applies to
all experiments, in which the count of each item follows poisson distribution.

Histogram records the observed results in an initial experiment. It is used to
train a constructed continued fraction. The continued fraction estimates the
number of additional distinct items discovered when another experiment conducts
under the similar conditions as the inital experiment. Then by adding the value
to the number of distinct items in the inital experiment, we could estimate the
total number of distinct items in a new experiment with bigger sample size.
}
\value{
  \item{continued.fraction}{
	  A list contains five components of a continued fraction: ps.coeffs is the 
	  coeffiencies of a power series, which the continued fraction tries to 
	  estimate; cf.coeffs is coeffiencies of the continued fraction; 
	  offset.coeffs is offset coeffiencies of the continued fraction; di is the 
	  diagonal value of the continued fraction; de is the degree of the 
	  continued fracton.
  }
  \item{yield.estimates}{
	  A two-column table represents estimated number of distinct items given 
	  experiments with different sample sizes.
  }
}
\references{
	\url{http://smithlabresearch.org/software/preseq/}
}
\author{
	Chao Deng
}
\note{
	  The format of a constructed continued fraction is totally defined by its 
	  coefficients, its diagonal value and offset coefficients. For example,
	  when the diagonal value is zero, there are no offset coefficients. The 
	  format of a continued fraction should like this:
	  \deqn{\frac{a_0t}{1 + \frac{a_1t}{1 + \frac{a_2t}{\ldots}}}}
	  See Baker, G. & Graves-Morris, P. \emph{Pade Approximants} for the detail
	  of the continued fraction.
}

\examples{
## load library
library(preseqR)

## import data
data(ShakespeareWordCounts)

## convert a histogram into a count vector 
hist.count = vector(mode="numeric", length=max(ShakespeareWordCounts$frequencies))
hist.count[ShakespeareWordCounts$frequencies] = ShakespeareWordCounts$word.types

## construct a continued fraction and estimate the number of observed word types
## if different quantities of Shakespeare's works are discovered.
result = preseqR.continued.fraction.estimate(hist.count)

## print the constructed continued fraction
result$continued.fraction

## print estimated the number of word types under various quantities of 
## Shakespeare's works
result$yield.estimates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ approximation }
\keyword{ continued }
\keyword{ fraction }
