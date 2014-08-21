\name{preseqR.rfa.curve}
\alias{preseqR.rfa.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict the number of species through a continued fraction
}
\description{
The function uses rational function approximation to estimate the expected
number of new species, which did not appear in an initial survey, but will be
detected if sampling is continued.
}
\usage{
preseqR.rfa.curve(hist, di = 0, mt = 100, ss = NULL, max.extrapolation = NULL,
step.adjust = TRUE, header = FALSE, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
    Frequencies of count data in an initial survey. \code{hist} can be either a 
    file name of the data or a variable of a matrix containing the data. The 
    data must be a two-column table. Each row (j, n_j) in the table represents
    there were n_j species with j individuals observed in each single species. 
    Both j and n_j are positive integers and the first column must be sorted in
    an ascending order.
}
  \item{header}{
    A boolean variable number representing whether or not the frequencies of 
    count data file contains a header.
}
  \item{di}{
    An integer representing the diagonal value for a constructed continued 
    fraction approximation. Default value is 0.
}
  \item{mt}{
    An positive integer representing the maximum degree allowed in the continued
    fraction approximation. Default is 100.
}
  \item{ss}{
    An positive double representing the step size between samples. Default value
    is the size of the sample in the initial survey.
}
  \item{max.extrapolation}{
    A positive double representing the maximum possible size of a sample in a
    survey. Default value is 100 times the size of the sample in the initial
    survey.
}
  \item{step.adjust}{
    A boolean value representing whether or not the value of \code{ss} could 
    adjust if the user-defined value of \code{ss} is relatively small comparing
    with the number of individuals observed in the initial survey. If 
    \code{step.adjust} is TRUE, \code{ss} could adjust to around one twentieth
    size of the initial sample when the user-defined \code{ss} is too small.
  }
  \item{seed}{
    A single value, interpreted as an integer, or NULL for random number 
    generator. Default is NULL. Set the seed if you want to reproduce the 
    results.
}
}
\details{
    The function is a quick mode of
    \code{preseqR.rfa.species.accumulation.curve}. It estimate the
    expected number of species through rational function approximation without
    bootstrapping. 
}
\value{
  \item{continued.fraction}{
    A list contains five components of a continued fraction approximation:
    ps.coeffs is the coeffiencies of a power series, which is the continued
    fraction approximation to; cf.coeffs is coeffiencies of the continued
    fraction approximation; offset.coeffs is offset coeffiencies of the 
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
Good, I. J., and G. H. Toulmin. "The number of new species, and the increase 
in population coverage, when a sample is increased." Biometrika 43.1-2 (1956):
45-63.

Efron, Bradley, and Ronald Thisted. "Estimating the number of unseen species: 
How many words did Shakespeare know?." Biometrika 63.3 (1976): 435-447.

Daley, Timothy, and Andrew D. Smith. "Predicting the molecular complexity of 
sequencing libraries." Nature methods 10.4 (2013): 325-327.

Baker, G. & Graves-Morris, P. Pade Approximants (Cambrige University Press,
Cambridge, UK, 1996).
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
## Shakespeare's works.
result$estimates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rational Function Approximation }
