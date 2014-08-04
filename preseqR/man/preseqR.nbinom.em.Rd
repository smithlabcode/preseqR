\name{preseqR.zerotruncated.mle}
\alias{preseqR.zerotruncated.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to fit a negative binomial distribution through MLE
}
\description{
preseqR.zerotruncated.mle fits the counts of a item through a negative binomial 
distribution.}
\usage{
preseqR.zerotruncated.mle(hist, size = SIZE.INIT, mu = MU.INIT)
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
  \item{size}{
	  a positive real value to set the inital value of number of successful trials, or 
	  dispersion parameter (the shape parameter of the gamma mixing distribution).
	  Must be strictly positive, need not be integer. The default value 
	  SIZE.INIT = 1
  }
  \item{mu}{
	  a positive real value to set the initial mean value of a negative binomial
	  distribution. The default value MU.INIT = 0.5
  }
}
\details{
	  The function assumes that the read count follows a negative binomial 
	  distribution. Under the assumption, given a histogram, each observed read 
	  count should follow a zero truncated negative binomial distributon. Thus a
	  conditional maximum likelyhood estimation is used to estimate two 
	  parameters size and mu.

	  Both size and mu are parameters for the negative binomial distribution.
	  See \url{http://stat.ethz.ch/R-manual/R-patched/library/stats/html/NegBinomial.html}
	  for more information on definition of parameters. 
	  
	  The function uses optim() in R to find the optimal parameters. It returns
	  results from optim(). See Value in 
	  \url{http://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html}

}
\value{
  \item{par}{
	  The best set of parameters found.
	   }
  \item{value}{
	  The value of negative loglikelyhood  corresponding to par
  }
  \item{counts}{
	  A two-element integer vector giving the number of calls to fn and gr 
	  respectively. This excludes those calls needed to compute the Hessian, 
	  if requested, and any calls to fn to compute a finite-difference 
	  approximation to the gradient.
  }
  \item{convergence}{
	  An integer code.
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
hist.count = vector(mode="numeric", length=max(SRR1041830_5M_hist$frequencies))
hist.count[SRR1041830_5M_hist$frequencies] = SRR1041830_5M_hist$distinct.reads
## print the parameters of a fitting negative binomial distribution
preseqR.zerotruncated.mle(hist.count)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ zerotruncated }
\keyword{ negative }
\keyword{ binomial }
\keyword{ MLE }
