\name{preseqR.nbinom.em}
\alias{preseqR.nbinom.em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to fit a negative binomial distribution through MLE
}
\description{
preseqR.nbinom.em fits a histogram of observations with a negative binomial 
distribution.}
\usage{
preseqR.nbinom.em(hist, size = SIZE.INIT, mu = MU.INIT, header = FALSE)
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
  \item{header}{
	  A logic number represents whether or not histogram file contains a header.
  }
}
\details{
	  The function assumes that the read count follows a negative binomial 
	  distribution. The histogram gave frequencies of observed items. The number
	  of unobserved items is the missing data. The function applies an EM 
	  algorithm to estimate parameters in the negative binomial distribution. In
	  E-step, the number of unobserved items is estimated by its expected value,
	  under current parameters. Then in M-step, preseqR.nbinom.em uses 
	  conditional maximum loglikelyhood given the number of observations to
	  update parameters based on a complete data set. It stops when the
	  conditional loglikelyhood stops increasing.
}
\value{
  \item{size}{
	  Estimated value of the parameter size in a negative binomial distribution.
	   }
  \item{mu}{
	  Estimated value of the parameter mu in a negative binomial distribution.
  }
  \item{loglik}{
	  Conditional loglikelyhood given the number of observations under an 
	  estimated negative binomial distribution.
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
preseqR.nbinom.em(hist.count)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ zerotruncated }
\keyword{ negative }
\keyword{ binomial }
\keyword{ MLE }
