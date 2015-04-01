\name{preseqR.rfa.species.accum.curve}
\alias{preseqR.rfa.species.accum.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Predict the number of species with confidence intervals
}
\description{
    The function estimates the expected number of species as a function of the
    number of captures using rational function approximations to Good and
    Toulmin's power series. The capture counts histogram is bootstrapped to
    improve the accuracy of the estimated curve and build confidence intervals..
}
\usage{
preseqR.rfa.species.accum.curve(n, bootstrap.times = 100,
mt = 100, ss = NULL, max.extrapolation = NULL, conf = 0.95)
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
}
\details{
According to Good & Toulmin (1956) and Efron & Thisted (1976), under
a multinomial or independent compound Poisson model for the number of individuals
observed for each species in the population, 
a non-paramtric empirical Bayes estimator can be derived for the expected number
of new species if sampling continues in the form of an alternating power series
in t, with t equal to the relative increase in the number of individuals captured.
Coefficients of the power series are
estimated through the frequencies of the number of individuals of each observed
species in the initial survey. While it
performs well for small extrapolation, the power series shows large variance in
general when the size pasts twice the size of the initial survey.

Daley, T., & Smith, A. D. (2013) used rational function approximation to the
power series of Good, I. J., & Toulmin, G. H.. The rational function
approximation is locally close to the power series of Good & Toulmin but is
constructed to be globally stable. It can apply to both small and larger
extrapolation.

The confidence interval is estimated through a log normal confidence interval
based on Chao, A. (1987) formula 12.
}
\value{
    A four-column table representing estimated number of species given samples
    of different sizes or \code{NULL} representing failure of bootstrapping.
    The first column is the total size in captured individuals of a sample. 
    The second column is the estimated number of species. The third and fourth
    column are the lower and upper bounds, respectively, of the corresponding
    confidence intervals.

    NULL if bootstrapping failed.
}
\references{
Good, I. J., & Toulmin, G. H. (1956). The number of new species, and the
increase in population coverage, when a sample is increased.
Biometrika, 43(1-2), 45-63.

Efron, B., & Thisted, R. (1976). Estimating the number of unseen species:
How many words did Shakespeare know?. Biometrika, 63(3), 435-447.

Efron, B. (1979). Bootstrap methods: another look at the jackknife.
The annals of Statistics, 1-26.

Daley, T., & Smith, A. D. (2013). Predicting the molecular complexity of
sequencing libraries. Nature methods, 10(4), 325-327.

Chao, A. (1987). Estimating the population size for capture-recapture data with
unequal catchability. Biometrics, 783-791.

\url{http://smithlabresearch.org/software/preseq/}
}
\author{
  Chao Deng
}
\note{
    The rational fraction approximation is only applied to extrapolation. For 
    estimating the expected number of species given a sample with size less than
    the size of sample in the initial survey, we use \code{\link{preseqR.interpolate.distinct}} to calculate the value.

    A global variable \code{BOOTSTRAP.factor} defines maximum resampling times 
    allowed for bootstrapping. The default value is 0.1. When resampling times
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
preseqR.rfa.species.accum.curve(ShakespeareWordHist, bootstrap.times = 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rational Function Approximation }
