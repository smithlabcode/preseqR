\name{preseqR-package}
\alias{preseqR-package}
\alias{preseqR}
\docType{package}
\title{
    Species accumulation estimation through rational function approximations
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
Date: \tab 2014-08-19\cr
License: \tab GPL-3\cr
}

functions: preseqR.rational.function.approximation.estimate, 
           preseqR.bootstrap.species.accumulation.curve,
           preseqR.calculate.continued.fraction.approximation,
           preseqR.ztnb.em,
           preseqR.ztnb.estimate,
           preseqR.ztnb.species.accumulation.curve,

data: FisherButterflyCaptures, ShakespeareWordCounts
}
\author{
Chao Deng, Timothy Dalay, and Andrew D. Smith

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
## butterfly species as a function of total captures
res <- preseqR.rational.function.approximation.estimate(FisherButterflyCaptures)

## print the continued fraction approximation
res$continued.fraction

## print the estimated number of butterfly species
res$estimates

## calculate the rational function approximation value when the argument equals
## to one
preseqR.calculate.continued.fraction.approximation(res$continued.fraction, 1)


## use bootstrap method to estimate number of butterfly species and 
## the associated confidence intervals
preseqR.bootstrap.species.accumulation.curve(FisherButterflyCaptures)

## fit a zero-truncated negative binomial distribution to the capture counts
preseqR.ztnb.em(FisherButterflyCaptures)

## estimate the number of butterfly species when the size of sample
## reaches 10000 through ztnb
preseqR.ztnb.estimate(FisherButterflyCaptures, 10000)

## estimate the number of butterfly specie as a function of captures
preseqR.ztnb.species.accumulation.curve(FisherButterflyCaptures)

}
