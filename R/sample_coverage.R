#    Copyright (C) 2016 University of Southern California and
#             Chao Deng and Andrew D. Smith and Timothy Daley
#
#    Authors: Chao Deng
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#

## predict the sample coverage of a random sample
## the parameter r means the probability of a random sampled individual 
## for which the species is represented at least r times in the sammple
preseqR.sample.cov <- function(n, r=1, mt=20)
{
  checking.hist(n)

  n[, 2] <- as.numeric(n[, 2])

  ## coefficients of average discovery rate for the first mt terms
  PS.coeffs <- discoveryrate.ps(n, mt=mt)

  if (is.null(PS.coeffs)) {
    write("the size of the initial experiment is insufficient", stderr())
    return(NULL)
  }

  ## use nonzero coefficients
  mt <- min(mt, length(PS.coeffs))
  PS.coeffs <- PS.coeffs[ 1:mt ]

  ## check whether sample size is sufficient
  if (mt < 2)
  {
    m <- paste("max count before zero is less than min required count (2)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }
  ## require odd mt
  mt <- mt - ((mt + 1) %% 2)
  valid.estimator <- FALSE
  m <- mt
  while (valid.estimator == FALSE && m >= 1) {

    rfa <- rfa.sample.cov(n, mt=m)
    rfa <- rfa.simplify(rfa)

    if (is.null(rfa)) {
      m <- m - 2
      next
    }
    ## check stability
    if (any(Re(rfa$poles) >= 0)) {
      m <- m - 2
      next
    }

    coefs <- rfa$coefs
    poles <- rfa$poles
    ## check whether the estimator is non-decreased
    ## NOTE: it only checks for t >= 1 !!!
    deriv.f <- function(t) {
      Re(sapply(t, function(x) {-(coefs*poles) %*% ( 1 / ((x-poles)^2))}))}
    if (any( deriv.f(seq(1, 100, by=0.05)) < 0 )) {
      m <- m - 2
      next
    } else {
      f <- function(x) {
        term1 <- coefs * (x / (x - poles))^r
        term2 <- r / (x - poles) - 1 / poles
        -Re( term1 %*% term2 / sum(coefs / poles) )
      }
      f.sample.cov <- function(t) {sapply(t, function(x) return(f(x)))}
      valid.estimator <- TRUE
    }
  }

  if (valid.estimator == TRUE) {
    return(f.sample.cov)
  } else {
    return(NULL)
  }
}


## predict the fraction of k-mers represented at least r times in the sample
kmer.frac <- function(n, r=2, mt=20) {
  return(preseqR.sample.cov(n=n, r=r-1, mt=mt))
}


## the fraction of k-mers represented at least r times as a function of 
## sample sizes
kmer.frac.curve <- function(n, k, read.len, seq.gb, r=2, mt=20) {
  f <- kmer.frac(n, r=r, mt=mt)
  if (is.null(f))
    return(NULL)
  n[, 2] <- as.numeric(n[, 2])
  N <- n[, 1] %*% n[, 2]
  ## average number of k-mers per read
  m <- read.len - k + 1
  unit.gb <- N / m * read.len / 1e9
  seq.effort <- seq.gb / unit.gb
  result <- matrix(c(seq.gb, f(seq.effort)), ncol=2, byrow=FALSE)
  colnames(result) <- c("bases(GB)", paste("frac(X>=", r, ")", sep=""))
  return(result)
}
