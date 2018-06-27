\name{preseqR.rSAC.sequencing.rmdup}
\alias{preseqR.rSAC.sequencing.rmdup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Predicting \eqn{r}-SAC in WES/WGS
}
\description{
  \code{preseqR.rSAC.sequencing.rmdup} predicts the expected number of 
  nucleotides in the genome sequenced at least r times in a sequencing
  experiment, based on a shallow sequencing experiment.
}
\usage{
preseqR.rSAC.sequencing.rmdup(n_base, n_read, r=1, mt=20, times=30, conf=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n_base}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{N_j}, the number of nucleotides in the genome sequenced exactly 
    \eqn{j} times in the initial experiment. The first column must be sorted in
    an ascending order.
  }
  \item{n_read}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{N'_j}, the number of distinct reads with exactly \eqn{j} duplicates
    in the initial experiment. The first column must be sorted in an
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
    The number of bootstrap samples. Default is 30.
  }
  \item{conf}{
    The confidence level. Default is 0.95
  }
}
\details{
  \code{preseqR.rSAC.sequencing.rmdup} is designed for sequencing experiments, 
  where duplicate reads are removed. The procedure is commonly used in 
  whole-exome sequencing experiments and sometimes appeared in WGS as well.
  To use the function, one must have two histograms. The first histogram
  is the coverage histogram, which is based on distinct reads.
  The second histogram is the counts of reads with exactly \eqn{j} duplicates.
}
\value{
  \item{f}{
    The estimator for the expected number of nucleotides in the genome
    sequenced at least \eqn{r} times given the amount of sequencing.
    The input of the estimator is a vector of sequencing efforts t, i.e. 
    the relative amount of sequencing comparing with the amount in the initial
    experiment. For example, t = 2 means sequencing twice the amount of the initial
    experiment.
  }
  \item{se}{
    The standard error for the estimator. The input is a vector of sequencing
    efforts t.
  }
  \item{lb}{
    The lower bound of the confidence interval.The input is a vector of
    sequencing efforts t.
  }
  \item{ub}{
    The upper bound of the confidence interval.The input is a vector of
    sequencing efforts t.
  }
}

\references{
Deng, C., Daley, T., Calabrese, P., Ren, J., & Smith, A.D. (2016). Estimating
the number of species to attain sufficient representation in a random sample.
arXiv preprint arXiv:1607.02804v3.
}

\author{
  Chao Deng
}

\examples{
## load library
#library(preseqR)

## import data
# data(SRR1301329_1M_base)
# data(SRR1301329_1M_read)
# construct the estimator
# estimator1 <- preseqR.rSAC.sequencing.rmdup(
#                 n_base=SRR1301329_1M_base, n_read=SRR5365359_5M_read, 
#                 r=4, mt=20, times=100, conf=0.95)
## The number of nucleotides in the genome covered at least 4 times, when the
## amount of sequencing is 10 or 20 times of the intial experiment
## 10 or 20 times of the initial sample
# estimator1$f(c(10, 20))
## The standard error of the estiamtes
# estimator1$se(c(10, 20))
## The confidence interval of the estimates
# lb <- estimator1$lb(c(10, 20))
# ub <- estimator1$ub(c(10, 20))
# matrix(c(lb, ub), byrow=FALSE, ncol=2)

# construct the estimator
# estimator2 <- preseqR.rSAC.sequencing.rmdup(
#                 n_base=SRR1301329_1M_base, n_read=SRR5365359_5M_read, 
#                 r=10, mt=20, times=100, conf=0.95)
## The number of nucleotides in the genome covered at least 10 times, when the
## amount of sequencing is 10 or 20 times of the intial experiment
## 10 or 20 times of the initial sample
# estimator2$f(c(10, 20))
## The standard error of the estiamtes
# estimator2$se(c(10, 20))
## The confidence interval of the estimates
# lb <- estimator2$lb(c(10, 20))
# ub <- estimator2$ub(c(10, 20))
# matrix(c(lb, ub), byrow=FALSE, ncol=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{r-SAC, duplicates, WES, WGS, high-throughput}
