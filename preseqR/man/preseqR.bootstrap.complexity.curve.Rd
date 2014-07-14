\name{presesqR.bootstrap.complexity.curve}
\alias{preseqR.bootstrap.complexity.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	A function to bootstrap a complexity curve
}
\description{
	preseqR.bootstrap.complexity.curve estimates number of distinct items given
	various experiments. All experiments should be under the same condition as
	it is in the inital experiment but with different sizes. Bootstrap method is
	used to improve estimation and make a confident interval.	
}
\usage{
preseqR.bootstrap.complexity.curve(hist, bootstrap.times = 100, di = 0, mt = 100, 
ss = 1e+6, mv = 1e+10, max.extrapolation = 1e+10, step.adjust = TRUE)
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
  \item{mv}{
	  A positive real number represents the maximum total size of an experiment
	  that the constructed continued fraction can apply. The total size of an 
	  experiment is calculated by the size of an initial experiment plus the 
	  size of additional experiment.
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
	and applies \code{preseqR.continued.fraction.estimate}, which could return
	a valid estimation or NULL. When times of valid estimation are beyond the 
	bootstrap.times, it finish resampling and returns the mean of results from
    each estimation	as final estimation. When times of total resampling are
	beyond a upperbound, it stops and returns failture as the inital experiment
	is poor for bootstrappling. 
}
\value{
  \item{yield.estimates}{
	  A two-columns table represents estimated number of distinct items given
	  various sizes of experiments.}
  \item{LOWER_0.95CI}{Lower bound for a 95\% confident interval}
  \item{UPPER_0.95CI}{Upper bound for a 96\% confident interval}
}
\references{
\url{http://smithlabresearch.org/software/preseq/}
}
\author{
	Chao Deng
}
\note{
	A global variable MULTINOMIAL.SAMPLE.TIMES defines the number of histograms
	resampled in each iteration. The default value is 19. Each resampled 
	histogram will be used as the histogram of an inital experiment and make an
	estimation. Thus there are total MULTINOMIAL.SAMPLE.TIMES for each iteration.
	Proper MULTINOMIAL.SAMPLE.TIMES could accelerate the function.

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
## import data
data(SRR1041830_5M_hist)
## hist.count is a count vector of a histogram
hist.count = SRR1041830_5M_hist$distinct.reads
result = preseqR.bootstrap.complexity.curve(hist.count)
## print the estimated yields given various sequencing sample size
print(c("sample.size", "yield estimates", "LOWER_0.95CI", "UPPER_0.95CI"))
yields = c(result$yield.estimates$sample.size, result$yield.estimates$yields, 
           result$LOWER_0.95CI, result$UPPER_0.95CI)
matrix(yields, ncol = 4, byrow = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ complexity }
\keyword{ bootstrap }
