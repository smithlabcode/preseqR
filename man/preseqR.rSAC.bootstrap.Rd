\name{preseqR.rSAC.bootstrap}
\alias{preseqR.rSAC.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the number of species represented r or more times
}
\description{
The function estimates the expected number of species represented at least r
times in a random sample based on the initial sample. The initial sample is
then bootstrapped to improve the stability of estimates.
}
\usage{
preseqR.rSAC.bootstrap(n, r=1, mt=20, size=SIZE.INIT, mu=MU.INIT, times=30,
                       conf=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with each species represented \eqn{j}
    times in the initial sample. The first column must be sorted in an
    ascending order.
  }
  \item{r}{
    A positive integer. Default is 1.
  }
  \item{mt}{
    An positive integer constraining possible rational function
    approximations. Default is 20.
  }
  \item{times}{
    The number of bootstrap samples.
  }
  \item{size}{
    A positive double, the initial value of the parameter \code{size} in
    the negative binomial distribution for the EM algorithm. Default value is 1.
  }
  \item{mu}{
    A positive double, the initial value of the parameter \code{mu} in the
    negative binomial distribution for the EM algorithm. Default value is 0.5.
  }
  \item{conf}{
    The confidence level. Default is 0.95
  }
}
\details{
  This is the bootstrap version of \code{\link{preseqR.rSAC}}. The function
  provides the confidence interval for estiamtes based on
  \code{\link{preseqR.rSAC}}. The confidence interval is based on the lognormal
  model proposed by Chao (1987).
}
\value{
  \item{f}{
    The estimator for the number of species represented at least r times in 
    a random sample. The input of the estimator is a vector of sampling 
    efforts t, i.e. the relative sample sizes comparing with the initial sample.
    For example, t = 2 means a random sample that is twice the size of the 
    initial sample.
  }
  \item{v}{
    The variance for the estimator. The input is a vector of sampling 
    efforts t. The variance depends on t.
  }
  \item{lb}{
    The lower bound of the confidence interval.The input is a vector of sampling 
    efforts t.
  }
  \item{ub}{
    The upper bound of the confidence interval.The input is a vector of sampling 
    efforts t.
  }
}

\references{
Chao, A. (1987). Estimating the population size for capture-recapture data with
unequal catchability. Biometrics, 783-791.

Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.

Kalinin V (1965). Functionals related to the poisson distribution and 
statistical structure of a text. Articles on Mathematical Statistics and 
the Theory of Probability pp. 202-220.

Daley, T., & Smith, A. D. (2013). Predicting the molecular complexity of
sequencing libraries. Nature methods, 10(4), 325-327.

Deng, C., Daley, T., & Smith, A. (2015). Applications of species accumulation
curves in large-scale biological data analysis. Quantitative Biology, 3(3),
135-144.

Deng, C and Smith, AD (2016). Estimating the number of species to attain 
sufficient representation in a random sample. arXiv preprint arXiv:1607.02804
}

\author{
  Chao Deng
}

\examples{
## load library
#library(preseqR)

## import data
#data(FisherButterflyHist)

## construct the estimator for the number of species represented at least once
## in a random sample
# estimator1 <- preseqR.rSAC.bootstrap(FisherButterflyHist, r=1)
## The number of species represented at least once, when the sample size is 
## 10 or 20 times of the initial sample
# estimator1$f(c(10, 20))
## The standard error of the estiamtes
# sqrt(estimator1$v(c(10, 20)))
## The confidence interval of the estimates
# lb <- estimator1$lb(c(10, 20))
# ub <- estimator1$ub(c(10, 20))
# matrix(c(lb, ub), byrow=FALSE, ncol=2)

## construct the estimator for the number of species represented at least twice
## in a random sample
# estimator2 <- preseqR.rSAC.bootstrap(FisherButterflyHist, r=2)
## The number of species represented at least twice, when the sample size is 
## 50 or 100 times of the initial sample
# estimator2(c(50, 100))
## The standard error of the estiamtes
# sqrt(estimator2$v(c(50, 100)))
## The confidence interval of the estimates
# lb <- estimator2$lb(c(50, 100))
# ub <- estimator2$ub(c(50, 100))
# matrix(c(lb, ub), byrow=FALSE, ncol=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Estimator, r-SAC, RFA, Bootstrap, Nonparametric}
