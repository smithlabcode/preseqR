\name{preseqR.ztnb.estimate}
\alias{preseqR.ztnb.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict the number of species through a ZTNB
}
\description{
This function calculates the expected number of new species observed
if the experiment is increased by t, relative to the number of
individuals observed in hist.
}
\usage{
preseqR.ztnb.estimate(n, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    Frequencies of the number of individuals of each species captured. 
    The data must be a two-column table.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with \eqn{j} individuals observed in the
    sample. The first column must be sorted in an ascending order.
}
  \item{t}{
    t is the relative size of the new survey to the initial survey. 
}
}
\details{
    A zero-truncated Negative Binomial distribution is fit to the number of
    individuals for each observed species. 
    We use \code{preseqR.ztnb.em} to estimate the
    parameters. 
    The function value is the gain in new species as the survey size increases.
	
	For extrapolation for the number of species through a 
	zero-truncated negative binomial distribution, 
	see supplementary of Daley, T., & Smith, A. D. (2013).
}

\value{
    a double as the function value. 
}
\author{
    Chao Deng
}

\seealso{
    \code{\link{preseqR.ztnb.em}}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
## load library
library(preseqR)

## import data
data(FisherButterflyHist)

## estimate the number of new butterfly species when the size of 
## the new sample is twice the size of the initial observed sample
preseqR.ztnb.estimate(FisherButterflyHist, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ZTNB }
