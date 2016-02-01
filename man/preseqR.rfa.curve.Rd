\name{preseqR.rfa.curve}
\alias{preseqR.rfa.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the number of species represented in a random sample
}
\description{
The function estimates the expected number of species represented in a random
sample using rational function approximation to Good and Toulmin's power series.
The only difference between this function and 
\code{\link{preseqR.rfa.species.accum.curve}} is that no bootstrapping for the 
initial sample. Therefore the function can not construct confidence intervals.
However, the function is much faster. In general, we recommand
\code{\link{preseqR.rfa.species.accum.curve}} for estimating the expected number
of species represented in a sample. 
}
\usage{
preseqR.rfa.curve(n, mt = 100, ss = NULL, max.extrapolation = NULL, asym.linear=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    A two-column matrix.
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column 
    is \eqn{n_j}, the number of species with each species represented by \eqn{j} 
    individuals in the sample. The first column must be sorted in an ascending order.
}
  \item{mt}{
    An positive integer equal to the maximum degree allowed in the continued
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
  \item{asym.linear}{
		A logic, setting to be TRUE if the asymptotic behavior of the accumulation curve
    is linear. Otherwise setting it to be FALSE.
	}  
}
\details{
    This function estimates the
    expected number of species through rational function approximation without
    bootstrapping. 
}
\value{
  \item{continued.fraction}{
    A list contains three components of a continued fraction approximation:
    ps.coeffs is the coeffiencients of the Good and Toulmin power series; 
    cf.coeffs is the coeffiencients of the associated continued
    fraction approximation; de is the degree of the continued fraction
    approximation.
  }
  \item{estimates}{
    A two-column matrix representing for estimates of the expected number of 
    species represented in a random sample. The first column is the size of the
    random sample; the second column is estimates of the expected number of species
    represented in the sample. 
  }
  \item{NULL}{
    NULL if the function failed to construct a rational functiona approximation.
}
}
\references{
Good, I. J., & Toulmin, G. H. (1956). The number of new species, and the
increase in population coverage, when a sample is increased.
Biometrika, 43(1-2), 45-63.

Efron, B., & Thisted, R. (1976). Estimating the number of unseen species:
How many words did Shakespeare know?. Biometrika, 63(3), 435-447.

Daley, T., & Smith, A. D. (2013). Predicting the molecular complexity of
sequencing libraries. Nature methods, 10(4), 325-327.

Baker, G. A., & Graves-Morris, P. (1996). Pade approximants
(Encyclopedia of Mathematics and its Applications vol 59).
}
\author{
  Chao Deng
}
\note{
    The format of a continued fraction approximation is totally determined by
    its coefficients, its diagonal value and its offset coefficients. For 
    example, when the diagonal value is zero, there are no offset coefficients.
    The format of a continued fraction should like this:
    \deqn{\frac{a_0(t-1)}{1 + \frac{a_1(t-1)}{1 + \frac{a_2(t-1)}{\ldots}}}}

    This is a special case of \code{\link{preseqR.pf.mincount}} by setting
    the parameter r=1. 
}

\seealso{
  \code{\link{preseqR.rfa.species.accum.curve}}
}

\examples{
## load library
library(preseqR)

## import data
data(ShakespeareWordHist)

## rational function approximation and estimate the number of observed word types
## as a function of total words
result = preseqR.rfa.curve(ShakespeareWordHist)

## print the constructed continued fraction approximation.
result$continued.fraction

## print the estimates
result$estimates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rational Function Approximation }
