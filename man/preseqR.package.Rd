\name{preseqR-package}
\alias{preseqR-package}
\alias{preseqR}
\docType{package}
\title{
    An R package for estimation of species accumulation
}
\description{
    preseqR provides functions to estimate species accumulation. Functions work
    through either rational function approximations to Good and Toulmin's (1956)
    non-parametric empirical Bayes power series estimator or fitting a 
    zero-truncated Negative Binomial distibution to the observed species counts.
    Both methods are described in Daley and Smith (2013).
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

Efron, B., & Thisted, R. (1976). Estimating the number of unseen species:
How many words did Shakespeare know?. Biometrika, 63(3), 435-447.

Good, I. J., & Toulmin, G. H. (1956). The number of new species, and the
increase in population coverage, when a sample is increased.
Biometrika, 43(1-2), 45-63.

Daley, T., & Smith, A. D. (2013). Predicting the molecular complexity of
sequencing libraries. Nature methods, 10(4), 325-327.

Baker, G. A., & Graves-Morris, P. (1996). Pade approximants
(Encyclopedia of Mathematics and its Applications vol 59).

Efron, B. (1979). Bootstrap methods: another look at the jackknife.
The annals of Statistics, 1-26.

Chao, A. (1987). Estimating the population size for capture-recapture data with
unequal catchability. Biometrics, 783-791.
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
