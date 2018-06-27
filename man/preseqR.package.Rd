\name{preseqR-package}
\alias{preseqR-package}
\docType{package}
\title{
    Predicting \eqn{r}-species accumulation curves
}
\description{
    The functionality of this package is to predict \eqn{r}-species
    accumulaiton curves. The method is based on a nonparametric empirical Bayes approach
    with rational function approximation. The estimator is
    excellent in accuracy for both large values of \eqn{r} and long-range
    extrapolations, which are essential to large-scale applications. Some
    examples are predicting the molecular complexity of sequencing
    libraries, estimating the minimum sufficient sequencing depths for
    whole-exome sequencing experiments and optimizing depths for single-cell
    whole-genome sequencing experiments. 
}
\details{
main functions: 

preseqR.rSAC

preseqR.rSAC.bootstrap

preseqR.optimal.sequencing

preseqR.rSAC.sequencing.rmdup

preseqR.sample.cov

preseqR.sample.cov.bootstrap
}
\author{
Chao Deng, Timothy Daley, and Andrew D. Smith

Maintainer: Chao Deng <chaodeng@usc.edu>
}
\references{
Baker, G. A., & Graves-Morris, P. (1996). Pade approximants
(Encyclopedia of Mathematics and its Applications vol 59).

Boneh, S., Boneh, A., & Caron, R. J. (1998). Estimating the prediction function
and the number of unseen species in sampling with replacement.
Journal of the American Statistical Association, 93(441), 372-379.

Chao, A., & Shen, T. J. (2004). Nonparametric prediction in species sampling.
Journal of agricultural, biological, and environmental statistics, 9(3), 253-269.

Cohen Jr, A. C. (1960). Estimating the parameters of a modified Poisson distribution.
Journal of the American Statistical Association, 55(289), 139-143.

Daley, T., & Smith, A. D. (2013). Predicting the molecular complexity of
sequencing libraries. Nature methods, 10(4), 325-327.

Deng C, Daley T & Smith AD (2015). Applications of species accumulation curves
in large-scale biological data analysis. Quantitative Biology, 3(3), 135-144. 
URL http://dx.doi.org/10.1007/s40484-015-0049-7.

Deng, C., Daley, T., Calabrese, P., Ren, J., & Smith, A.D. (2016). Estimating
the number of species to attain sufficient representation in a random sample.
arXiv preprint arXiv:1607.02804v3.

Efron, B., & Thisted, R. (1976). Estimating the number of unseen species:
How many words did Shakespeare know?. Biometrika, 63(3), 435-447.

Efron, B. (1979). Bootstrap methods: another look at the jackknife.
The annals of Statistics, 1-26.

Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.

Fisher, R. A., Corbet, A. S., and Williams, C. B. ,1943, The Relation Between 
the Number of Species and the Number of Individuals in a Random Sample of an  
Animal Population, Journal of Animal Ecology, 12, 42-58.

Good, I. J., & Toulmin, G. H. (1956). The number of new species, and the
increase in population coverage, when a sample is increased.
Biometrika, 43(1-2), 45-63.

Heck Jr, K. L., van Belle, G., & Simberloff, D. (1975). 
Explicit calculation of the rarefaction diversity measurement and the 
determination of sufficient sample size. Ecology, 1459-1461.

Kalinin V (1965). Functionals related to the poisson distribution and               
statistical structure of a text. Articles on Mathematical Statistics and            
the Theory of Probability pp. 202-220.
}
