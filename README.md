preseqR
=======

Code in this repository aims to make the functionality of Preseq available
in the R statistical computing environment. There are two ways this is
supposed to work:

  1.  The basic functionality of the preseq program, initially focusing only
      on library complexity, will be available. The specific function names
      will be indicated below.
  2.  The mathematical routines for doing rational function approximation via
      continued fractions will be implemented as a wrapper for our existing
      functionality in C++.

See <http://cran.r-project.org/web/packages/preseqR/index.html> for details.

UPDATES FROM PREVIOUS RELEASE
========================================================================
We have added an option to make use of the asymptotic behavior of the target
accumulation curve.
