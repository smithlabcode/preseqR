\name{preseqR.rfa.species.accum.curve}
\alias{preseqR.rfa.species.accum.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Predicting the number of species in a random sample
}
\description{
  The function estimates the expected number of species represented in a random
  sample using rational function approximation to Good and Toulmin's power series.
  The initial sample is bootstrapped to improve the stability of the estimated
  curve and construct confidence intervals. 
}
\usage{
preseqR.rfa.species.accum.curve(n, bootstrap.times = 100,
mt = 100, ss = NULL, max.extrapolation = NULL, conf = 0.95, asym.linear=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column 
    is \eqn{n_j}, the number of species with each species represented by \eqn{j} 
    individuals in the initial sample. The first column must be sorted in an
    ascending order.
}
  \item{bootstrap.times}{
    An positive integer representing the minimum required times of successful
    estimation. Default is 20. See detail below.
}
  \item{mt}{
    An positive integer equal to the maximum degree allowed in a continued
    fraction approximation. Default is 100.
}
  \item{ss}{
    An positive double equal to the step size between samples. Default value
    is the size of the initial sample.
}
  \item{max.extrapolation}{
    A positive double equal to the maximum possible size of a random sample.
    Default value is the 100 times the size of the initial sample.
}
  \item{conf}{
    A positive double in (0, 1) equal to the confidence level.
    Default value is 0.95.
  }
  \item{asym.linear}{
    A logic, setting to be TRUE if the asymptotic behavior of the accumulation curve
    is linear. Otherwise setting it to be FALSE.
  }
}
\details{
According to Good & Toulmin (1956) and Efron & Thisted (1976), under
a multinomial or independent compound Poisson model for the number of individuals
represented for each species in the population, 
a non-paramtric empirical Bayes estimator can be derived for the expected number
of new species if sampling continues in the form of an alternating power series
in t, with t equal to the relative increase in the number of individuals captured.
Coefficients of the power series are
estimated through an initial sample. While the method
performs well for small extrapolation, the power series shows large variance in
general when t > 2.

Daley, T., & Smith, A. D. (2013) used rational function approximation to the
power series of Good, I. J., & Toulmin, G. H.. The rational function
approximation is locally close to the power series of Good & Toulmin but is
constructed to be globally stable. It can apply to both small and larger
extrapolation.

The confidence interval is estimated through a log normal confidence interval
based on Chao, A. (1987) formula 12.
}
\value{
  A four-column matrix for estimates of the expected number of species
  represented in a random sample. The first column is the size of the
  random sample; the second column is estimates of the expected number of species
  represented in the sample. The third and fourth column are the lower
  and upper bounds, respectively, of the corresponding confidence intervals.

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
    The rational fraction approximation can be only applied to extrapolation. For 
    estimating the expected number of species in a random sample of size less than
    the size of the initial sample, we use 
    \code{\link{preseqR.interpolate.distinct}} to calculate the value.

    A global variable \code{BOOTSTRAP.factor} defines maximum resampling times 
    allowed for bootstrapping. The default value is 0.4. 
    
    When resampling times are greater than
    \code{bootstrap.times} / \code{BOOTSTRAP.factor}, the function will terminate.

    This is a special case of \code{\link{preseqR.pf.mincount.bootstrap}} by
    setting the parameter r=1. 
}

\section{Warning}{
  The default setting for bootstrap.times (20) is not realiable in constructing
  the confidence interval. 
} 


\examples{
## load library
# library(preseqR)

## import data
# data(ShakespeareWordHist)

## estimate the number of unique words in a random sample
## minimum required successful estimation times is 100
# preseqR.rfa.species.accum.curve(ShakespeareWordHist, bootstrap.times = 100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rational Function Approximation, Bootstrap }
