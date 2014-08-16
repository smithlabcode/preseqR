\name{preseqR.bootstrap.species.richness}
\alias{preseqR.bootstrap.species.richness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Predict the number of species with a confident interval
}
\description{
    The function estimates expected number of species for samples of different
    sizes. Bootstrap method is applied to improve estimated results and build 
    confident intervals.  
}
\usage{
preseqR.bootstrap.species.richness(hist, bootstrap.times = 100, di = 0, mt = 100, 
ss = NULL, max.extrapolation = NULL, step.adjust = TRUE, header = FALSE, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
    Frequencies of count data in an initial survey. \code{hist} can be either a 
    file name of the data or a variable of a matrix containing the data. The 
    data must be a two-column table. Each row (j, n_j) in the table represents
    there were n_j species with j individuals observed in each single species. 
    Both j and n_j are positive integers and the first column must be sorted in
    an ascending order.
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
    An integer representing the diagonal value for a continued fraction 
    approximation. Default value is 0.
}
  \item{mt}{
    An positive integer representing the maximum degree allowed in a continued
    fraction approximation. Default is 100.
}
  \item{ss}{
    An positive double representing the step size between samples. Default value
    is the size of the sample in the initial survey.
}
  \item{max.extrapolation}{
    A positive double representing the maximum possible size of a sample in a
    survey. Default value is the 100 times the size of a sample in the initial 
    survey.
}
  \item{step.adjust}{
    A boolean value representing whether or not the value of \code{ss} could 
    adjust if the user-defined value of \code{ss} is relatively small comparing
    with the number of individuals observed in the initial survey. If 
    \code{step.adjust} is TRUE, \code{ss} could adjust to around one twentieth
    size of the initial sample when the user-defined \code{ss} is too small.
}
  \item{seed}{
    A single value, interpreted as an integer, or NULL for random number 
    generator. Default is NULL. Set the seed if you want to reproduce the 
    results.
}
}
\details{
In general contexts, "species" could be defined in a board way. For example, 
"species" could be reads in a sequencing experiment or words in a book. The 
function can apply to any defined "species" once the count of each single 
"species" follows a poisson distribution with rates possibly varying between
species.

According to Good, I. J., & Toulmin, G. H. (1956) and Efron & Thisted (1976), if
the count of each single species follows a compound poisson distribution, a 
power series could be derived to estimate the number of new species when
sampling continues. Coefficients of the power series are estimated through the
frequencies of count data in the initial survey. While it performs well for
small extrapolation, the power series shows large variance in general when the
size pasts twice the size of the initial survey.

Daley, T., & Smith, A. D. (2013) use rational function approximation to the
power series. The rational function approximation is locally close to the power
series of Good & Toulmin but is constructed to be globally stable. It can apply
to both small and larger extrapolation.

The confident interval is estimated through a log normal confident interval
based on Chao, A. (1987) formula 12.
}
\value{
    A four-column table representing estimated number of species given samples
    of different sizes or \code{NULL} representing failure of bootstrapping.
    The first column is the size of a sample. The second column is the estimated
    number of species. The third column is the lower bound of the 95\% confident
    interval. The fourth column is the upper bound of the 95\% confident 
    interval.
}
\references{
Good, I. J., and G. H. Toulmin. "The number of new species, and the increase 
in population coverage, when a sample is increased." Biometrika 43.1-2 (1956):
45-63.

Efron, Bradley, and Ronald Thisted. "Estimating the number of unseen species: 
How many words did Shakespeare know?." Biometrika 63.3 (1976): 435-447.

Efron, Bradley. "Bootstrap methods: another look at the jackknife." The annals
of Statistics (1979): 1-26.

Daley, Timothy, and Andrew D. Smith. "Predicting the molecular complexity of 
sequencing libraries." Nature methods 10.4 (2013): 325-327.

Chao, Anne. "Estimating the population size for capture-recapture data with
unequal catchability." Biometrics (1987): 783-791.

\url{http://smithlabresearch.org/software/preseq/}
}
\author{
  Chao Deng
}
\note{
    The rational fraction approximation is only applied to extrapolation. For 
    estimating the expected number of species given a sample with size less than
    the size of sample in the initial survey, we downsample the observed 
    individuals and count the number of species as estimation.

    A global variable \code{BOOTSTRAP.factor} defines maximum resampling times 
    allowed for bootstrapping. The default value is 0.1. When resampling times
    are beyond \code{bootstrap.times} / \code{BOOTSTRAP.factor}, the function 
    will terminate since the quality of the sample in the initial survey is
    poor.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## load library
library(preseqR)

## import data
data(FisherButterflyCaptures)

## estimate the number of distinct butterflies if various sizes of butterflies
## are captured
preseqR.bootstrap.species.richness(FisherButterflyCaptures)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rational Function Approximation }
