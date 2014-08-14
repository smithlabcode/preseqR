\name{preseqR-package}
\alias{preseqR-package}
\alias{preseqR}
\docType{package}
\title{
    Species richness estimation through rational function approximations
}
\description{
    This is a R package to make the functionality of Preseq available in the R
    statistical computing environment. It provides functions to construct,
    calculate and use rational function approximations to estimate the number of
    observed species for a sample with a given size.
}
\details{
\tabular{ll}{
Package: \tab preseqR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-08-12\cr
License: \tab GPL-3\cr
}

functions: preseqR.rational.function.approximation.estimate, 
           preseqR.bootstrap.species.richness,
           preseqR.calculate.continued.fraction.approximation,
           preseqR.ztnb.em,
           preseqR.ztnb.estimate,
           preseqR.ztnb.species.richness,
           preseqR.print2file;

data: FisherButterflyCaptures, ShakespeareWordCounts
}
\author{
Andrew D. Smith, Timothy Daley and Chao Deng

Maintainer: Chao Deng <chaodeng@usc.edu>
}
\references{
Fisher, R. A., Corbet, A. S., and Williams, C. B. ,1943, The Relation Between 
the Number of Species and the Number of Individuals in a Random Sample of an  
Animal Population, Journal of Animal Ecology, 12, 42-58.

Efron, Bradley, and Ronald Thisted. "Estimating the number of unseen species: 
How many words did Shakespeare know?." Biometrika 63.3 (1976): 435-447.

Good, I. J., and G. H. Toulmin. "The number of new species, and the increase 
in population coverage, when a sample is increased." Biometrika 43.1-2 (1956):
45-63.

Daley, Timothy, and Andrew D. Smith. "Predicting the molecular complexity of 
sequencing libraries." Nature methods 10.4 (2013): 325-327.

Baker, G. A., and Peter Graves-Morris. "Pade approximants (Encyclopedia of 
Mathematics and its Applications vol 59)." (1996).

Efron, Bradley. "Bootstrap methods: another look at the jackknife." The annals
of Statistics (1979): 1-26.

Chao, Anne. "Estimating the population size for capture-recapture data with
unequal catchability." Biometrics (1987): 783-791.
}

\examples{
## load library
library(preseqR)

## load data
data(FisherButterflyCaptures)

## construct a continued fraction approximation and extrapolate the number of
## butterfly species given samples of different sizes
res <- preseqR.rational.function.approximation.estimate(FisherButterflyCaptures)

## output the continued fraction approximation
res$continued.fraction

## calculate the function value when the argument equals to one
preseqR.calculate.continued.fraction.approximation(res$continued.fraction, 1)

## output the estimated number of butterflies species
res$yield.estimates

## write results down to files
preseqR.print2file(res$continued.fraction, 'continued_fraction.txt')
preseqR.print2file(res$yield.estimate, 'butterflies_estimation.txt')

## use bootstrap method to build confident intervals
res <- preseqR.bootstrap.species.richness(FisherButterflyCaptures)

## write results down to a file
preseqR.print2file(res, 'butterfliesCI.txt')

## fit the data with a zero-truncated negative binomial distribution
preseqR.ztnb.em(FisherButterflyCaptures)

## estimate the number of butterfly species when the size of sample
## reaches 10000 through ztnb
preseqR.ztnb.estimate(FisherButterflyCaptures, 10000)

## estimate the number of butterfly specie given samples of different sizes
res <- preseqR.ztnb.species.richness(FisherButterflyCaptures)

## write results down to a file
preseqR.print2file(res, 'butterflies_ZTNB_estimation.txt')
}
