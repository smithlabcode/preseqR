\name{preseqR.zerotruncated.complexity.curve}
\alias{preseqR.zerotruncated.complexity.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	A function to estimate a complexity curve
}
\description{
	preseqR.ztnb.complexity.curve estimates number of distinct items given
	various experiments. All experiments should be under the same condition as
	it is in the inital experiment but with different sizes. Observed read count
	is fitted into a negative binomial and the estimated number of distinct
	items for each experiment is its expectation. An EM algorithm 
	preseqR.nbinom.em in preseqR is used to estiamte the parameters in the
   	negative binomial distributin.
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
	  A positive real number represents the maximum total size of an experiment
	  that the constructed continued fraction can apply. The total size of an  
	  experiment is calculated by the size of an initial experiment plus the 
	  size of additional experiment.
}
}
\details{
	See \code{\link{preseqR.zerotruncated.estimate}} for the estimation method.
}
\value{
  \item{sample.size}{
	  Various sizes of experiments.
  }
  \item{yield.estimates}{
	  Estimated number of distinct items given various sizes of experiments.
  }
}
\author{
	Chao Deng
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## import data
data(SRR1041830_5M_hist)
## hist.count is a count vector of a histogram
hist.count = vector(mode = 'numeric', length = max(SRR1041830_5M_hist$frequencies))
hist.count[SRR1041830_5M_hist$frequencies] = SRR1041830_5M_hist$distinct.reads
## estimate a complexity curve
res <- preseqR.ztnb.complexity.curve(hist.count)
## print the estimated yields given various sequencing sample size
print(c("sample.size", "yield estimates"))
yields = c(res$sample.size, res$yield.estimates)
matrix(yields, ncol = 2, byrow = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ complexity }
\keyword{ negative }
\keyword{ binomial }
