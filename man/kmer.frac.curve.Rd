\name{kmer.frac.curve}
\alias{kmer.frac.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fraction of \eqn{k}-mers observed at least \eqn{r} times 
}
\description{
\code{kmer.frac.curve} predicts the expected fraction of \eqn{k}-mers observed at
least \eqn{r} times in a high-throughput sequencing experiment given the 
amount of sequencing
}
\usage{
  kmer.frac.curve(n, k, read.len, seq, r=2, mt=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{N_j}, the number of \eqn{k}-mers observed exactly \eqn{j}
    times in the initial experiment. The first column must be sorted in an
    ascending order.
  }
  \item{k}{
    The number of nucleotides in a \eqn{k}-mer.
  }
  \item{read.len}{
    The average length of a read.
  }
  \item{seq}{
    The amount of nucleotides sequenced..
  }
  \item{r}{
    A positive integer. Default is 1.
  }
  \item{mt}{
    An positive integer constraining possible rational function
    approximations. Default is 20.
  }
}
\details{
  \code{kmer.frac.curve} is mainly designed for metagenomics to evaluate how 
  saturated a metagenomic data is.

  \code{kmer.frac.curve} is the fast version of \code{\link{kmer.frac.curve.bootstrap}}.
  The function does not provide the confidence interval. To obtain the
  confidence interval along with the estimates, one should use the function
  \code{\link{kmer.frac.curve.bootstrap}}.
}
\value{
  A two-column matrix. The first column is the amount of sequencing in an 
  experiment. The second column is the estimate of the fraction of \eqn{k}-mers observed at least
  \eqn{r} times in the experiment.
}
\references{
Deng, C and Smith, AD (2016). Estimating the number of species to attain 
sufficient representation in a random sample. arXiv preprint arXiv:1607.02804
}
\author{
  Chao Deng
}

\examples{
## load library
library(preseqR)

## import data
data(SRR061157_k31)

## the fraction of 31-mers represented at least 10 times in an experiment when
## sequencing 1M, 10M, 100M, 1G, 10G, 100G, 1T nucleotides
kmer.frac.curve(n=SRR061157_k31, k=31, read.len=100, seq=10^(6:12), r=10, mt=20) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{k-mer, sample coverage, high-throughput, metagenomics}
