\name{preseqR.optimal.sequencing}
\alias{preseqR.optimal.sequencing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal amount of sequencing for scWGS
}
\description{
\code{preseqR.optimal.sequencing} predicts the optimal amount of sequencing in
a single-cell whole-genome sequencing (scWGS) experiment based on a shallow sequencing experiment.
}
\usage{
preseqR.optimal.sequencing(n, efficiency=0.05, bin=1e8, r=1, mt=20, 
                           times=30, conf=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with each species represented exactly \eqn{j}
    times in the initial sample. The first column must be sorted in an
    ascending order.
  }
  \item{efficiency}{
    The minimum benefit-cost ratio
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
  \item{conf}{
    The confidence level. Default is 0.95
  }
}
\details{
  \code{preseqR.optimal.sequencing} predicts the optimal amount of sequencing
  in a scWGS experiment. The term optimal is interpreted as the maximum 
  amount of sequencing with its benefit-cost ratio greater than a given threshold.
  The benefit-cost ratio is defined as the probability of a new nucleotide in the
  genome represented at least \eqn{r} times when one more base is sequenced. 
  In order to improve the numeric stability, we use the mean of new nucleotdies
  with coverage at least \eqn{r} in one unit of sequencing effort to approximate the
  ratio. The amount of sequences in one unit of sequencing effort is defined by
  the variable \code{bin}.

  Note that the benefit-cost ratio is not monotonic. The ratio first increases
  and then decrease as the amount of sequencing increase. To predicte the
  optimal amount of sequencing, we consider only the areas after the peak, where
  the ratio starts to decrease.
}
\value{
  A vector of three dimensions. The first coordinate is the optimal amount of
  sequencing. The second and the third coordinates are the lower and upper bound
  of the confidence interval.
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
# data(SRR611492_5M)
## the optimal amount of sequencing with the benefit-cost ratio greater than 
## 0.05 for r = 4
# preseqR.optimal.sequencing(n=SRR611492_5M, efficiency=0.05, bin=1e8, r=4)
## the optimal amount of sequencing with the benefit-cost ratio greater than 
## 0.05 for r = 10
# preseqR.optimal.sequencing(n=SRR611492_5M, efficiency=0.05, bin=1e8, r=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{WGS, high-throughput, benefit-cost ratio}
