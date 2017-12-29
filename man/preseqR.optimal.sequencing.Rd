\name{preseqR.optimal.sequencing}
\alias{preseqR.optimal.sequencing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicting the optimal amount of sequencing in a single-cell whole-genome
sequencing (scWGS) experiment
}
\description{
The function predicts the optimal amount of sequencing in a scWGS experiment
based on the cost-benefit ratio. The term optimal is interpreted as the maximum 
amount of sequencing with its cost-benefit ratio greater than a given threshold.
}
\usage{
preseqR.optimal.sequencing(n, efficiency=0.05, bin=1e8, r=1, mt=20, 
                           size=SIZE.INIT, mu=MU.INIT, times=30, conf=0.95)
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
  \item{efficiency}{
    The minimum cost-benefit ratio
  }
  \item{bin}{
    One unit of sequencing effort. Default is 1e8.
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
  The cost-benefit ratio is defined as the probability of a new position in the
  genome represented at least r times when one more base is sequenced. In order
  to improve the numeric stability, we use the mean of new positions 
  with coverage at least r in one unit of sequencing effort to approximate the
  ratio. The amount of sequences in one unit of sequencing effort is defined by
  the variable \code{bin}.
  The function predicts the maximum amount of
  sequencing such that the cost-benefit ratio is greater than a predefined value. 
  Note that the cost-benefit ratio is not monotonic. The ratio first increases
  and then decrease as the amount of sequencing increase. For predicting the
  optimal amount of sequencing, we consider only the part where the ratio starts
  to decrease.
}
\value{
  A vector of three dimensions. The first coordinate is the optimal amount of
  sequencing. The second and the third coordinates are the lower and upper bound
  of the confidence interval.
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
#library(preseqR)

## import data
#data(SRR611492_5M)
## the optimal amount of sequencing with the cost-benefit ratio greater than 
## 0.05 for r = 4
#preseqR.optimal.sequencing(n=SRR611492_5M, efficiency=0.05, bin=1e8, r=4)
## the optimal amount of sequencing with the cost-benefit ratio greater than 
## 0.05 for r = 10
#preseqR.optimal.sequencing(n=SRR611492_5M, efficiency=0.05, bin=1e8, r=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DNA, High-throughput sequencing}
