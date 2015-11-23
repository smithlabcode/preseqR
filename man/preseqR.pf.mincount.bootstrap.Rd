\name{preseqR.pf.mincount.bootstrap}
\alias{preseqR.pf.mincount.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Predict the number of species observed r or more times with confidence intervals
}
\description{
    The function estimates the number of species observed r or more times as a function of the
    number of captures. The capture counts histogram is bootstrapped to
    improve the accuracy of the estimated curve and build confidence intervals.
}
\usage{
preseqR.pf.mincount.bootstrap(n, bootstrap.times = 100,
mt = 100, ss = NULL, max.extrapolation = NULL, conf = 0.95, r=1)
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
  \item{bootstrap.times}{
    An positive integer representing the minimum required number of successful
    estimations. See detail below.
}
  \item{mt}{
    An positive integer equal to the maximum degree allowed in a continued
    fraction approximation. Default is 100.
}
  \item{ss}{
    An positive double equal to the step size between samples. Default value
    is the size of the sample in the initial survey.
}
  \item{max.extrapolation}{
    A positive double equal to the maximum possible size of a sample in a
    survey. Default value is the 100 times the size of the initial 
    survey.
}
  \item{conf}{
    A positive double in (0, 1) equal to the confidence level.
    Default value is 0.95.
  }
  \item{r}{
      A vector of positive integers equal to minimum counts.
  }
}
\details{
    This function estimates the expected number of species observed at least r 
    times through rational function approximation. The bootstrap method is used
    to construct the confidence interval. 
}
\value{
    A list contains estimates of the number of species observed at least r times 
    given samples of different sizes or  \code{NULL} representing failure of bootstrapping.
    Each element in the list is a four-column table representing estimated number of species given samples
    of different sizes or \code{NULL} representing failure of bootstrapping.
    The first column is the total size in captured individuals of a sample. 
    The second column is the estimated number of species. The third and fourth
    column are the lower and upper bounds, respectively, of the corresponding
    confidence intervals. Each  element in the list corresponds to one minimum count in the vector r. 

    NULL if bootstrapping failed.
}
\references{
}
\author{
  Chao Deng
}
\note{
    The estimator is only applied to extrapolation. For 
    estimating the expected number of species at least r times given a sample with size less than
    the size of sample in the initial survey, we use \code{\link{preseqR.interpolate.mincount}} to calculate the value.

    A global variable \code{BOOTSTRAP.factor} defines maximum resampling times 
    allowed for bootstrapping. The default value is 0.4. When resampling times
    are beyond \code{bootstrap.times} / \code{BOOTSTRAP.factor}, the function 
    will terminate.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## load library
library(preseqR)

## import data
data(ShakespeareWordHist)

## setting the random seed
set.seed(123456)
## estimate the number of distinct butterflies as a function of the number of
## captured butterflies. Minimum required successful estimation times is 10
preseqR.pf.mincount.bootstrap(ShakespeareWordHist, bootstrap.times = 10, r=c(1,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Minimum count }
