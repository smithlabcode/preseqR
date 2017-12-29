UPDATES TO VERSION 4.0.0
========================
1. Improve the user interface for core functions
2. Add functions to optimize the depth of single-cell whole-genome sequencing
experiments and whole-exome sequencing experiments
3. Add functions to predict the sample coverage, which is the probability of
sampling an observed species from a population
4. Add functions to predict the fraction of k-mers represented at least r times
in a sequencing experiment

UPDATES TO VERSION 3.1.2
========================
1. Fix a bug for removing defects

UPDATES TO VERSION 3.1.1
========================
1. Substitute embedded c++ code with R code
2. Remove the dependencies on the software preseq 

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

INSTALLATION
============

1. We recommand everyone to install the package preseqR from CRAN. 
It can be easily done by opening an R shell and typing:

   ``>install.packages("preseqR")``

2. The following instructions are for installing the package from the source. 
Assume the source code of preseqR has been pulled from the git repo and it is 
under the current directory. Open an R shell and type:
   
   ``>install.packages("polynom")``
   
   ``>install.packages("preseqR", repos=NULL, type="source")``
   
   Note that the package polynom is required by preseqR.

