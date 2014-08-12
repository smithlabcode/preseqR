\name{presesqR.bootstrap.complexity.curve}
\alias{preseqR.bootstrap.complexity.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Predict the number of species with a confident interval
}
\description{
    The function estimates expected number of species given samples of different
    sizes. Bootstrap method is applied to improve estimated results and build 
    confident intervals.  
}
\usage{
preseqR.bootstrap.complexity.curve(hist, bootstrap.times = 100, di = 0, mt = 100, 
ss = NULL, max.extrapolation = NULL, step.adjust = TRUE, header = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
    Frequencies of count data in an inital survey. \code{hist} can be either a 
    file name of the data or a variable of a matrix containing the data. The 
    data must be a two-column table. Each row (x, y) in the table represents
    there were y species with x individuals observed in each single species. 
    Both x, y are positive integers and the first column must be sorted in an 
    ascending order.
}
  \item{header}{
    A boolean variable number representing whether or not the frequencies of 
    count data file contains a header.
}
  \item{bootstrap.times}{
    An positive integer representing the minimum required number of successful
    estimations. See detail below.
}
  \item{di}{
    An integer representing the diagonal value for a constructed continued 
    fraction. Default value is 0.
}
  \item{mt}{
    An positive integer representing the maximum degree allowed in the continued
    fraction approximation. Default is 100.
}
  \item{ss}{
    An positive double representing the step size between samples. Default value
    is the size of the sample in the initial survey.
}
  \item{max.extrapolation}{
    A positive double representing the maximum possible size of a sample in a
    survey. Default value is the \code{ss} times 100.
}
  \item{step.adjust}{
    A boolean value representing whether or not the value of \code{ss} could 
    adjust if the user-defined value of \code{ss} is relatively small comparing
    with the number of individuals observed in the inital survey. If 
    \code{step.adjust} is TRUE, \code{ss} could adjust to around one twentieth
    size of the inital sample when the user-defined \code{ss} is too small.
}
}
\details{
    First, the frequencies of count data in the initial survey is resampled. 
    Then we consider new generated data as the initial survey and apply 
    \code{preseqR.continued.fraction.estimate}, which could return either a 
    valid continued fraction approximation or NULL. When the number of valid
    estimations are beyond the \code{bootstrap.times}, 
    \code{presesqR.bootstrap.complexity.curve} terminates and returns the median
    of estimations as the expected number of species. However, if resampling 
    times have already reached the threshold, while the number of valid 
    estimations is still below the \code{bootstrap.times}, the function will
    stop and return a \code{NULL}. 

    The confident interval is estimated based on Chao, A. (1987) formula 12.
}
\value{
    A four-column table representing estimated number of species given samples
    of different sizes or \code{NULL} representing failure of boostrapping.
    The first column is the size of a sample. The second column is the estimated
    number of species. The third column is the lower bound of the 95\% confident
    interval. The fourth column is the upper bound of the 95\% confident 
    interval.
}
\references{
Efron, Bradley. "Bootstrap methods: another look at the jackknife." The annals
of Statistics (1979): 1-26.

Chao, Anne. "Estimating the population size for capture-recapture data with
unequal catchability." Biometrics (1987): 783-791.

\url{http://smithlabresearch.org/software/preseq/}
}
\author{
  Chao Deng
}
\note{
  A global variable \code{BOOTSTRAP.factor} defines maximum resampling times 
  allowed for bootstrapping. The default value is 0.1. When resampling times are
  beyond \code{bootstrap.times} / \code{BOOTSTRAP.factor}, the function will
  terminate since the quality of the sample in the inital survery is poor.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{preseqR.continued.fraction.estimate}}
}
\examples{
## load library
library(preseqR)

## import data
data(FisherButterflyCaptures)

## estimate the number of distinct butterflies if various sizes of butterflies
## are trapped
preseqR.bootstrap.complexity.curve(FisherButterflyCaptures)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ complexity }
\keyword{ bootstrap }
\keyword{ confident }
\keyword{ interval }
