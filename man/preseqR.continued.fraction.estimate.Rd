\name{preseqR.continued.fraction.estimate}
\alias{preseqR.continued.fraction.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to predict the number of species
}
\description{
The function constructs a continued fraction to estimate the number of new 
species, which did not appear in an intial survey, could be detected if sampling
could be continued.
}
\usage{
preseqR.continued.fraction.estimate(hist, di = 0, mt = 100, ss = NULL, 
max.extrapolation = NULL, step.adjust = TRUE, header = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
    Frequencies of frequencies data in an inital survey. \code{hist} can be 
    either a file name of the data or a variable of a matrix containing the 
    data. The data must be a two-column table. Each row (x, y) in the table
    represents there were y species with x individuals observed in each single
    species. Both x, y are positive integers and the first column must be sorted
    in an ascending order.
}
  \item{header}{
    A logic number represents whether or not the frequencies of frequencies data
    file contains a header.
}
  \item{di}{
    An integer represents the diagonal value for a constructed continued 
    fraction.
}
  \item{mt}{
    An integer represents the maximum number of parameters allowed in a 
    continued fraction.
}
  \item{ss}{
    An integer represents the step size of samples. Default value is the size of
    a sample in the initial survey.
}
  \item{max.extrapolation}{
    A positive real number represents the maximum possible size of a sample in a
    survey. Default value is the \code{ss} times 100.
}
  \item{step.adjust}{
    A logic value represents whether or not the value of \code{ss} could adjust
    if the user-defined value of \code{ss} is relatively small comparing with 
    the number of individuals observed in the inital survey. If 
    \code{step.adjust} is TRUE, \code{ss} could adjust to around one twentieth
    size of the inital sample when the user-defined \code{ss} is too small.
  }
}
\details{
In general contexts, "species" could be defined in a board way. For example, 
"species" could be reads in a sequencing experiment or words in a book. The 
function can apply to any defined "species" once the count of each single 
"species" follows a poisson distribution.

According to Good, I. J., & Toulmin, G. H. (1956), if the count of each single 
species follows a poisson distribution, a power series could be derived to 
estimate the number of new species when sampling continues. Coefficients of the
power series are estimated through the frequencies of frequencies data in the
initial survey. While it performs well for small extrapolation, the power series
shows large variance in general when the size pasts twice the size of the 
initial survey.

Daley, T., & Smith, A. D. (2013) use a continued fraction to approximate the
power series. The constructed continued fraction locally closes to the power 
series but has a larger radius of convergence. It can apply to both small and
larger extrapolation.
}
\value{
  \item{continued.fraction}{
    A list contains five components of a continued fraction: ps.coeffs is the 
    coeffiencies of a power series, which the continued fraction tries to 
    estimate; cf.coeffs is coeffiencies of the continued fraction; offset.coeffs
    is offset coeffiencies of the continued fraction; di is the diagonal value
    of the continued fraction; de is the degree of the continued fracton.
  }
  \item{yield.estimates}{
    A two-column table represents estimated number of species given samples of
    different sizes.
  }
}
\references{
Good, I. J., and G. H. Toulmin. "The number of new species, and the increase 
in population coverage, when a sample is increased." Biometrika 43.1-2 (1956):
45-63.

Daley, Timothy, and Andrew D. Smith. "Predicting the molecular complexity of 
sequencing libraries." Nature methods 10.4 (2013): 325-327.
}
\author{
  Chao Deng
}
\note{
    The format of a constructed continued fraction is totally determined by its 
    coefficients, its diagonal value and its offset coefficients. For example,
    when the diagonal value is zero, there are no offset coefficients. The 
    format of a continued fraction should like this:
    \deqn{\frac{a_0t}{1 + \frac{a_1t}{1 + \frac{a_2t}{\ldots}}}}
    See Baker, G. & Graves-Morris, P. \emph{Pade Approximants} for the detail
    of the continued fraction.
}

\examples{
## load library
library(preseqR)

## import data
data(ShakespeareWordCounts)

## construct a continued fraction and estimate the number of observed word types
## given different quantities of Shakespeare's works
result = preseqR.continued.fraction.estimate(ShakespeareWordCounts)

## print the constructed continued fraction
result$continued.fraction

## print estimated the number of word types under different quantities of 
## Shakespeare's work
result$yield.estimates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ approximation }
\keyword{ continued }
\keyword{ fraction }
