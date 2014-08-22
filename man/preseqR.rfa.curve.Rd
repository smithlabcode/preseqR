\name{preseqR.rfa.curve}
\alias{preseqR.rfa.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the number of species observed as a function of captures
with RFA to the NPEB estimator
}
\description{
This function estimates the expected number of species observed as a function of
the number of total individuals captured. It uses rational function
approximations to Good and Toulmin's (1965) non-parametric empirical Bayes power
series.
}
\usage{
preseqR.rfa.curve(hist, di = 0, mt = 100, ss = NULL, max.extrapolation = NULL,
step.adjust = TRUE, header = FALSE, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
    Frequencies of the number of individuals of each species captured. 
    It can be either a file name of the data or a variable of a matrix 
    containing the data. The data must be a two-column table.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of species with \eqn{j} individuals observed in the
    sample. The first column must be sorted in an ascending order.
}
  \item{header}{
    A boolean variable number representing whether or not the frequencies of 
    count data file contains a header.
}
  \item{di}{
    An integer equal to the diagonal value for a constructed continued 
    fraction approximation. Default value is 0.
}
  \item{mt}{
    An positive integer equal to the maximum degree allowed in the continued
    fraction approximation. Default is 100.
}
  \item{ss}{
    An positive double equal to the step size between samples. Default value
    is the size of the sample in the initial survey.
}
  \item{max.extrapolation}{
    A positive double equal to the maximum possible size of a sample in a
    survey. Default value is 100 times the size of the sample in the initial
    survey.
}
  \item{step.adjust}{
    A boolean value representing whether or not the value of \code{ss} can be
    adjusted if the user-defined value of \code{ss} is small compared to
    the total number of individuals observed in the initial survey. If 
    \code{step.adjust} is TRUE, \code{ss} will adjust to around one twentieth
    size of the initial sample when the user-defined \code{ss} is too small.
  }
  \item{seed}{
    A single value interpreted as an integer or NULL for random number 
    generator. Default is NULL.
}
}
\details{
    This function estimates the
    expected number of species through rational function approximation without
    bootstrapping. 
}
\value{
  \item{continued.fraction}{
    A list contains five components of a continued fraction approximation:
    ps.coeffs is the coeffiencients of the Good and Toulmin power series; 
    cf.coeffs is the coeffiencients of the associated continued
    fraction approximation; offset.coeffs is the offset coeffiencients of the 
    continued fraction approximation; di is the diagonal value of the continued
    fraction approximation; de is the degree of the continued fraction
    approximation.
  }
  \item{estimates}{
    A two-column table representing estimated number of species given samples of
    different sizes.
  }
  \item{NULL}{
    NULL if the rational function approximation failed.
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
    \deqn{\frac{a_0t}{1 + \frac{a_1t}{1 + \frac{a_2t}{\ldots}}}}
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
## as a function of # words in Shakespeare's works.
result = preseqR.rfa.curve(ShakespeareWordHist)

## print the constructed continued fraction approximation.
result$continued.fraction

## print the number of observed word types as a function of # words
## in Shakespeare's works.
result$estimates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rational Function Approximation }
