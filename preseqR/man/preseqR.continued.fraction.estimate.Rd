\name{preseqR.continued.fraction.estimate}
\alias{preseqR.continued.fraction.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to predict library complexity
}
\description{
preseqR.continued.fraction.estimate creates a continued fraction to estimate 
the number of distinct molecules, reads, or species given its histogram. It also
provides a complexity curve to describe the complexity of the capture-recapture
data.
}
\usage{
preseqR.continued.fraction.estimate(hist, di = 0, mt = 100, ss = 1e+06, 
		     mv = 1e+10, max.extrapolation = 1e+10, step.adjust = TRUE)
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
  \item{di}{
	  Diagonal value for a constructed continued fraction.
}
  \item{mt}{
	  Maximum number of parameters in a continued fraction.
}
  \item{ss}{
	  Step size of sampling points along a library.
}
  \item{mv}{
	  The maximum value to train a continued fraction.
}
  \item{max.extrapolation}{
	  The maximum possible value to extrapolate
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
  \item{CF}{All components of a continued fraction. ps.coeffs is the 
	  coeffiencies of a power series, which the continued fraction estimates;
	  co.coeffs is the coeffiencies of the continued fraction; offset.coeffs is
  	  offset coeffiencies of the contined fraction; di is the diagonal value of
      the continued fraction; de is the degree of the contined fracton}
  \item{yield.estimates}{Estimated yields of distinct reads if a experiment with 
        sample size from the first column occurs}
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
\keyword{ approximation }
