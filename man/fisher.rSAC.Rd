\name{fisher.rSAC}
\alias{fisher.rSAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Logseries estimator
}
\description{
\code{fisher.rSAC} estimates the expected number of species represented at least 
\eqn{r} times in a random sample, based on the initial sample.
The estimator was originally proposed by Fisher et al. (1943) for estimating
the SAC. We generalize this estimator for predicting the \eqn{r}-SAC.
}
\usage{
fisher.rSAC(n, r=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{N_j}, the number of species with each species represented exactly \eqn{j}
    times in the initial sample. The first column must be sorted in an
    ascending order.
  }
  \item{r}{
    A positive integer. Default is 1.
  }
}
\value{
  The estimator for the \eqn{r}-SAC. The input of the estimator is a vector of
  sampling efforts \eqn{t}, i.e., the relative sample sizes comparing with the initial
  sample. For example, \eqn{t = 2} means a random sample that is twice the size of
  the initial sample.
}
\author{
  Chao Deng
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{
Fisher, R., Corbet, A., & Williams, C. (1943). The Relation Between the Number 
of Species and the Number of Individuals in a Random Sample of an Animal 
Population. Journal of Animal Ecology, 12(1), 42-58. doi:10.2307/1411  

Deng, C., Daley, T., Calabrese, P., Ren, J., & Smith, A.D. (2016). Estimating
the number of species to attain sufficient representation in a random sample.
arXiv preprint arXiv:1607.02804v3.
}

\examples{
## load library
library(preseqR)

## import data
data(WillButterfly)

## construct the estimator for SAC
fisher1 <- fisher.rSAC(WillButterfly, r=1)
## The number of species represented at least once in a sample, 
## when the sample size is 10 or 20 times of the initial sample
fisher1(c(10, 20))

## construct the estimator for r-SAC
fisher2 <- fisher.rSAC(WillButterfly, r=2)
## The number of species represented at least twice in a sample, 
## when the sample size is 50 or 100 times of the initial sample
fisher2(c(50, 100))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ estimator, r-SAC, parametric }
