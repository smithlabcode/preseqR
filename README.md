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
