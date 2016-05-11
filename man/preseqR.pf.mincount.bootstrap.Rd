\name{preseqR.pf.mincount.bootstrap}
\alias{preseqR.pf.mincount.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the number of species represented r or more times
}
\description{
    The function estimates the expected number of species represented r or more 
    times in a random sample based on the initial sample. The initial sample is
    bootstrapped to improve the stability of estimates and construct confidence
    intervals.
}
\usage{
preseqR.pf.mincount.bootstrap(n, bootstrap.times = 100, mt = 100, ss = NULL, 
                              max.extrapolation = NULL, conf = 0.95, r=1)
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
    An positive integer representing the minimum required number of successful
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
    Default value is 100 times the size of the initial sample.
}
  \item{conf}{
    The confidence level. Default value is 0.95.
  }
  \item{r}{
    A vector of positive integers. Default is 1.
  }
}
\details{
    Under a mixture of Poisson model for the number of individuals represented
    in a random sample for each species in a population, the expected number of
    species represented at least r times can be expressed as higher derivatives
    of the accumulation curve. We use rational function approximations to the 
    modified Good and Toulmin's (1956) non-parametric empirical Bayes power 
    series to construct an estimator for the acccumulation curve. we obtain the 
    estimate for the number of species represented at least r times in a random
    sample by differentiating the curve. 

    The confidence interval is estimated through a log normal confidence interval 
    based on Chao, A. (1987) formula 12. 
}
\value{
    A list of four-column matrices representing estimates of the expected number
    of species represented r or more times in a random sample. Each matrix 
    in the list corresponds to a specified value of r. The first column of a matrix is 
    the size of a random sample; the second column is the estimate of the
    number of species represented r or more times in the sample. The third 
    and fourth column are the lower and upper bounds, respectively, of the
    corresponding confidence intervals. 

    NULL if bootstrapping failed.
}

\references{
Kalinin V (1965). Functionals related to the poisson distribution and 
statistical structure of a text. Articles on Mathematical Statistics and 
the Theory of Probability pp. 202-220.

Daley, T., & Smith, A. D. (2013). Predicting the molecular complexity of
sequencing libraries. Nature methods, 10(4), 325-327.
}

\author{
  Chao Deng
}
\note{
    The estimator based on rational function approximations can be only applied
    to extrapolation. For estimating the expected number of species at least r times
    in a random sample of size less than the size of the initial sample, we use 
    \code{\link{preseqR.interpolate.mincount}} to calculate the value.

    A global variable \code{BOOTSTRAP.factor} defines maximum resampling times 
    allowed for bootstrapping. The default value is 0.4. When resampling times
    are greater than \code{bootstrap.times / BOOTSTRAP.factor}, 
    the function will terminate.
}

\section{Warning}{
  The default setting for bootstrap.times is not realiable in constructing the 
  confidence interval. 
} 

\examples{
## load library
# library(preseqR)

## import data
# data(ShakespeareWordHist)

## estimate the number of unique words appeared at least once and twice 
## as a function of the number of words
# result = preseqR.pf.mincount.bootstrap(ShakespeareWordHist, r=c(1,2, 20))
## estimates of the number of unique words appeared at least once
# result[[1]]
## estimates of the number of unique words appeared at least twice
# result[[2]]
## estimates of the number of unique words appeared at least twenty times
# result[[3]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sufficient representation, Partial fraction, Bootstrap }
