\name{preseqR-package}
\alias{preseqR-package}
\alias{preseqR}
\docType{package}
\title{
    An R package for estimation of species accumulation
}
\description{
    preseqR provides functions to estimate species accumulation by rational function
    approximations to Good and Toulmin's (1956) non-parametric empirical Bayes
    power series estimator and by fitting a zero-truncated Negative Binomial
    distibution to the observed species counts, both described in Daley
    and Smith (2013).
}
\details{
functions: 

preseqR.rfa.species.accum.curve, 

preseqR.ztnb.species.accum.curve,

preseqR.rfa.curve,

preseqR.rfa.estimate,

preseqR.ztnb.estimate,

preseqR.ztnb.em,

data: FisherButterflyHist, ShakespeareWordHist
}
\author{
Chao Deng, Timothy Daley, and Andrew D. Smith

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
data(FisherButterflyHist)

## estimate the accumulation of butterfly species as a function
## with rational function approximations to Good and Toulmin's
## power series.
res <- preseqR.rfa.curve(FisherButterflyHist)

## print the rational function approximation as a 
## continued fraction
res$continued.fraction

## print the butterfly species accumulation curve
res$estimates

## calculate the rational function approximation value for 
## the expected number of new butterfly species observed when 
## the experiment is doubled and 3,506 more butterflies are
## captured, following the convention of Efron and Thisted 
## (1976) that t = 0 corresponds to the initial experiment
preseqR.rfa.estimate(res$continued.fraction, t = 1)

## Estimate the butterfly species accumulation curve by
## bootstrapping the observed species count histogram
preseqR.rfa.species.accum.curve(FisherButterflyHist)

## fit a zero-truncated negative binomial distribution to 
## the observed species counts
preseqR.ztnb.em(FisherButterflyHist)

## calculate the expected number of new butterfly species
## observed when the experiment is doubled under the 
## Negative Binomial model
preseqR.ztnb.estimate(FisherButterflyHist, 1)

## estimate the number of butterfly species observed as a 
## function of total captures under the Negative Binomial model
preseqR.ztnb.species.accum.curve(FisherButterflyHist)

}
