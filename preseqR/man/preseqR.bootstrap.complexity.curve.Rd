\name{presesqR.bootstrap.complexity.curve}
\alias{bootstrap.complexity.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	bootstrapping Complexity curve
}
\description{
	The function estimates the yields of distinct reads if sequencing experiments 
	occur under the same condition as the initial experiment but with different 
	sample sizes. The histogram of the initial experiment is provided as input. 
	It uses bootstrap to improve yields estimation and make a confidence interval.
}
\usage{
preseqR.bootstrap.complexity.curve(hist, bootstrap.times = 100, di = 0, mt = 100, 
ss = 1e+6, mv = 1e+10, max.extrapolation = 1e+10, step.adjust = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
	  A histogram. It can be either a file name or a count vector of the 
	  histogram. For histogram file, it contains two columns. The first column
	  is frequencies of sequenced reads. Values of frequency should be at least
	  one.  For each given frequency, the second column is the number of distinct 
	  reads sequenced with that frequency. A count vector is a non-negative
	  integral vector. The value x in row i represents x distinct reads are
	  sequenced i times. 
}
  \item{bootstrap.times}{
     The number of resampling times for bootstrapping.
}
  \item{di}{
	  Diagonal value for a constructed continued fraction.
}
  \item{mt}{
	  Maximum number of parameters in a continued fraction.
}
  \item{ss}{
	  Step size of sampling points along a sequenced library.
}
  \item{mv}{
	  The maximum value to train a continued fraction.
}
  \item{max.extrapolation}{
	  The maximum possible value to extrapolat.
}
  \item{step.adjust}{
	  Step size of sampling points could adapt in order to correspond to the size
	  of an inital experiment.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{yield.estimates}{Estimated yields of distinct reads if a experiment with 
	  sample size from the first column occurs}
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
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Complexity }
\keyword{ Library }% __ONLY ONE__ keyword per line
