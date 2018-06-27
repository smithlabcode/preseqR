\name{fisher.alpha}
\alias{fisher.alpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter alpha in the logseries estimator
}
\description{
\code{fisher.alpha} estimates the parameter alpha in the logseries estimator by
Fisher, R. A., et al. (1943) based on an initial sample.
}
\usage{
fisher.alpha(n)
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
}
\value{
  A double, the estimated value of the parameter alpha
}
\author{
  Chao Deng
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{
Fisher, R., Corbet, A., & Williams, C. (1943). The Relation Between the Number 
of Species and the Number of Individuals in a Random Sample of an Animal 
Population. Journal of Animal Ecology, 12(1), 42-58. doi:10.2307/1411  
}

\examples{
## load library
library(preseqR)

## import data
data(WillButterfly)

## estimating alpha
fisher.alpha <- fisher.alpha(WillButterfly)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ parametric }
