\name{preseqR.zerotruncated.estimate}
\alias{preseqR.zerotruncated.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to predict the number of distinct items
}
\description{
preseqR.zerotruncated.estimate fits a given histogram with a negative binomial
distribution. Thus it uses the distribution to estimate the number of distinct
items if a new experiment with different size conducts.
}
\usage{
preseqR.zerotruncated.estimate(hist, n)
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
  \item{n}{
	  A positive real number represents the size of a new experiment.
}
}
\details{
	The function applies preseqR.zerotruncated.mle to fit the given histogram
	with a negative binomial. Then it estimates the total number of distinct 
	items and updates parameters size and mu. The number of distinct items in
	the new experiment is estimated by its expectation. 
}
\value{
	The expected number of distinct items in the new experiment with size n.
}
\author{
	Chao Deng
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	  See Baker, G. \& Graves-Morris, P. \emph{Pade Approximants} about the
	  form of a continued fraction
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## import data
data(SRR1041830_5M_hist)
## hist.count is a count vector of a histogram
hist.count = vector(mode="numeric", length=max(SRR1041830_5M_hist$frequencies))
hist.count[SRR1041830_5M_hist$frequencies] = SRR1041830_5M_hist$distinct.reads
preseqR.zerotruncated.estimates(hist.count, 1e8)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ estimate }
\keyword{ negative }
\keyword{ binomial }
