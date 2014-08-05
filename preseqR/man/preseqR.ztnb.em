\name{preseqR.nbinom.em}
\alias{preseqR.nbinom.em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to fit a negative binomial distribution through MLE
}
\description{
The function fits a histogram of observations with a negative binomial 
distribution. Since the number of unobserved items missed in the histogram,
an EM algorithm is used to estimate the number of unobserved items and 
parameters (size, mu) in a negative binomial distribution.}
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
	  A positive real value sets the inital value of the parameter size. The
	  parameters of the negative binomial distribution should be initialized
	  when an em algorithm applies.
  }
  \item{mu}{
	  A positive real value sets the initial of the parameter mu. The
	  parameters of the negative binomial distribution should be initialized
	  when an em algorithm applies.
  }
  \item{header}{
	  A logic number represents whether or not histogram file contains a header.
  }
}
\details{
	  The function assumes that the read count follows a negative binomial 
	  distribution. The histogram gave frequencies of observed items. The number
	  of unobserved items is a missing data. The function applies an EM 
	  algorithm to estimate parameters in the negative binomial distribution. In
	  E-step, the number of unobserved items is estimated by its expected value 
	  under current parameters settings. Then in M-step, preseqR.nbinom.em uses 
	  maximum likelyhood estimation to update parameters based on the complete
	  data set. It terminates when the conditional loglikelyhood, given the 
	  number of observations, stop increasing.
}
\value{
  \item{size}{
	  Estimation of the parameter size in a negative binomial distribution.
   }
  \item{mu}{
	  Estimation of the value of the parameter mu in a negative binomial
	   distribution.
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
## load library
library(preseqR)

## import data
data(FisherButterlyCaptures)

## convert a histogram into a count vector
hist.count = vector(mode = 'numeric', length = max(FisherButterlyCaptures$n))
hist.count[FisherButterlyCaptures$n] = FisherButterlyCaptures$S

## print the parameters of a fitting negative binomial distribution
preseqR.nbinom.em(hist.count)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ zerotruncated }
\keyword{ negative }
\keyword{ binomial }
\keyword{ MLE }
\keyword{ EM }
