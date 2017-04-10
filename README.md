UPDATES TO VERSION 3.0.1
========================

1. Fix a bug in Chao's estimator
2. Fix issues for a Solaris C++ compiler.

UPDATES TO VERSION 3.0.0
========================

1. We have changed the return types of many functions in the package. These
functions no longer generate estimated accumulative curves. 
Instead, they return function types, which are estimators for the number
of species represented by at least r indivdiduals in a random sample. 

2. We added several estimators for predicting the number of species represented
by at least r individuals in a random sample

UPDATES TO VERSION 2.1.1
========================

We have changed the interfaces for most of our exported functions. We add new
estimators for the number of species represented by at least r individuals in
a random sample.

preseqR
=======

Code in this repository aims to expand the functionality of Preseq available in 
the R statistical computing enviroment. There are five ways this is supposed to
work:

  1.  The basic functionality of the preseq program, initially focusing only
      on library complexity, is available. These functions contain the 
      string "rfa" as part of their names.
  2.  The mathematical routines for doing rational function approximation via
      continued fractions is implemented as a wrapper for our existing
      functionality in C++.
  3.  Fitting a zero-truncated negative binomial distribution to the sample is
      available. These functions include the string "ztnb" as part of the names.
  4.  The simulation module is used to generate samples based on mixture of Poisson.
  5.  Extra functions are provided to estimate the number of species represented
      at least r times in a random sample.

See <https://cran.r-project.org/package=preseqR> for details.
