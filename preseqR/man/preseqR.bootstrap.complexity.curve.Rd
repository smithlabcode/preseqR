\name{presesqR.bootstrap.complexity.curve}
\alias{preseqR.bootstrap.complexity.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	A function to bootstrap a complexity curve
}
\description{
	The function estimates number of distinct items given experiments with 
	various sample size. All experiments should be under the similar condition
   	as it was in the inital experiment but with different sizes. Bootstrap 
	method is used to improve estimation and make a confident interval.	
}
\usage{
preseqR.bootstrap.complexity.curve(hist, bootstrap.times = 100, di = 0, mt = 100, 
ss = NULL, max.extrapolation = NULL, step.adjust = TRUE, header = FALSE)
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
	  A logic number represents whether or not histogram file contains a header.
}
  \item{bootstrap.times}{
	  An integer represents the number of successful estimation times for
	  bootstrapping.
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
	  A positive real number represents the maximum total size of an experiment
	  that the constructed continued fraction can apply. The total size of an  
	  experiment is calculated by the size of an initial experiment plus the 
	  size of additional experiment.
}
  \item{step.adjust}{
	  A logic value represents whether or not the step size could change. If
	  step.adjust equals TRUE, the value of step size (ss) could adjust in
	  respond to the size of an initial experiment.
}
}
\details{
	The function resamples a histogram as the histogram of an intial experiment
	and applies preseqR.continued.fraction.estimate, which could return
	a valid estimation or NULL. When times of valid estimation are beyond the 
	bootstrap.times, it finish resampling and returns the mean of results from
    each estimation	as final results. When total resampling times are beyond an
    upperbound, it stops and returns failture as the inital experiment
	is poor for bootstrapping. 
}
\value{
	  A four-column table represents estimated number of distinct items given
	  experiments with various sample sizes. The first column is the sample size.
	  The second column is estimated number of distinct items. The third column
	  is the lower bound of the 95\% confident interval. The fourth column is the
	  upper bound of the 95\% confident interval.
}
\references{
\url{http://smithlabresearch.org/software/preseq/}
}
\author{
	Chao Deng
}
\note{
	A global variable MULTINOMIAL.SAMPLE.TIMES defines the number of histograms
	generated through sampling in each iteration. The default value is 19. Each
   	resampled histogram will be used as the histogram of an inital experiment 
	and make an estimation. It applies R's vector implementation feature to 
	accelerate the function. Proper MULTINOMIAL.SAMPLE.TIMES could make the
	function faster.

	A global variable BOOTSTRAP.factor defines maximum resampling times allowed
	for bootstrapping. THe default value is 0.1. When resampling times are 
	beyond bootstrap.times / BOOTSTRAP.factor, the function will terminate as
	the experiment is too poor to bootstrap.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{preseqR.continued.fraction.estimate}}
}
\examples{
## load library
library(preseqR)

## import data
data(FisherButterlyCaptures)

## convert a histogram into a count vector
hist.count = vector(mode = 'numeric', length = max(FisherButterlyCaptures$n))
hist.count[FisherButterlyCaptures$n] = FisherButterlyCaptures$S

## estimate the number of distinct butterflies if various sizes of butterflies
## are trapped
preseqR.bootstrap.complexity.curve(hist.count)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ complexity }
\keyword{ bootstrap }
\keyword{ confident }
\keyword{ interval }
