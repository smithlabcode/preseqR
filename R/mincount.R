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
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

## the pair, a pole and a root, is a defect if the distance is less than 
## the value defined by the variable PRECISION
PRECISION <- 1e-3

### interpolating for r-SAC
### ss, step size
preseqR.interpolate.rSAC <- function(n, ss, r=1)
{
  checking.hist(n)

  n[, 2] <- as.numeric(n[, 2])
  ## sample size
  N <- n[, 1] %*% n[, 2]

  ## the number of species in the initial sample
  initial.distinct <- sum(n[, 2])
  step.size <- as.double(ss)

  ## l is the number of sampled points for interpolation
  l <- N / step.size

  ## the step size is too large or too small
  if (l < 1 || ss < 1 || r < 1)
    return(NULL)
  ## the step size is the sample size
  ## count the number of species observed r or more times
  else if (l == 1) {
    index <- which(n[, 1] >= r)
    result <- matrix(c(step.size, sum(n[index, 2])), ncol = 2, byrow = FALSE)
    colnames(result) <- c('sample.size', 'interpolation')
    return(result)
  }

  ## explicitly calculating the expected species observed at least r times
  ## based on sampling without replacement
  ## see K.L Heck 1975
  ## N, the number of individuals in a sample
  ## S, the number of species in a sample
  ## size, the size of the subsample
  expect.distinct <- function(n, N, size, S, r) {
    denom <- lchoose(N, size)
    p <- sapply(n[, 1], function(x) {
       sum(exp(lchoose(N - x, size - 0:(r-1)) + lchoose(x, 0:(r-1)) - denom))})
    return(S - p %*% n[, 2])
  }

  ## subsample sizes
  x <- step.size * ( 1:l )

  ## the number of species represented at least r times in a subsample
  yields <- sapply(x, function(x) {
      expect.distinct(n, N, x, initial.distinct, r)})

  result <- matrix(c(x, yields), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'interpolation')

  return(result)
}


## coefficients for the power series of E(S_1(t)) / t
discoveryrate.ps <- function(n, mt) {
  ## transform a histogram into a vector of frequencies
  hist.count <- vector(length=max(n[, 1]), mode="numeric")
  hist.count[n[, 1]] <- n[, 2]

  PS.coeffs <- sum(hist.count)
  change.sign <- 0

  ## preserve extra precision mt+1
  for (i in 1:(min(mt+1, length(hist.count)))) {
    PS.coeffs <- c(
            PS.coeffs, 
            (-1)^change.sign * hist.count[i] - PS.coeffs[length(PS.coeffs)])
    change.sign <- change.sign + 1
  }

  ## truncate at coefficients where it is zero
  zero.index <- which(PS.coeffs == 0)
  if (length(zero.index) > 0) {
    PS.coeffs[1:(min(zero.index) - 1)]
  } else {
    PS.coeffs
  }
}


ds.mincount <- function(n, r=1, mt=20)
{
  checking.hist(n)

  n[, 2] <- as.numeric(n[, 2])

  ## constructing the power series
  PS.coeffs <- discoveryrate.ps(n, mt=mt)

  if (is.null(PS.coeffs)) {
    write("the size of the initial experiment is insufficient", stderr())
    return(NULL)
  }

  ## use only the first mt terms in the power series
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

  ## construct the continued fraction approximation to the power seies
  cf <- ps2cfa(coef=PS.coeffs, mt=mt)
  rf <- cfa2rf(CF=cf)
  ## the length of cf could be less than mt
  ## even if ps do not have zero terms, coefficients of cf may have
  mt <- length(cf)
  ## select rational function approximants [m-1/m] m=mt, mt-2, ..., 2
  ## asymptotically ~ C / t
  mt <- mt - (mt %% 2)
  valid.estimator <- FALSE
  m <- mt
  while (valid.estimator == FALSE && m >= 2) {

    rfa <- rf2rfa(RF=rf, m=m)
    ## solving roots
    numer.roots <- solve(rfa[[1]])
    denom.roots <- solve(rfa[[2]])

    ## finite
    if (any(!is.finite(c(numer.roots, denom.roots)))) {
      m = m - 2
      next;
    }

    ## record roots in the numerator that are significantly similar to
    ## roots in the denominator
    tmp.roots <- c()

    ## simplify the rational function approximation
    ## two roots are same if the difference is less than the 
    ## predefined PRECISION
    if (length(denom.roots) > 0) {
      for (i in 1:length(denom.roots)) {
        if (length(numer.roots) > 0) {
          d <- Mod(denom.roots[i] - numer.roots)
          ind <- which.min(d)
          if (d[ind] < PRECISION) {
            numer.roots <- numer.roots[-ind]
            tmp.roots <- c(tmp.roots, denom.roots[i])
          }
        }
      }
    }

    ## roots in simplified RFA
    denom.roots <- denom.roots[!denom.roots %in% tmp.roots]

    ## convert roots from t - 1 to t
    roots <- denom.roots + 1

    ## pacman rule checking
    if (any(Re(roots) >= 0)) {
      m <- m - 2
      next
    } else {
      if (length(numer.roots) == 0) {
        poly.numer <- as.function(polynomial(1))
      } else {
        poly.numer <- as.function(poly.from.roots(numer.roots))
      }
      l <- length(denom.roots)
      ## treat polynomials in the rational function to be monic
      ## the difference to the original RFA is a multiplier C

      ## c_i in the estimator
      coefs <- sapply(1:l, function(x) {
        poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
      ## calculate the constant C
      C <- coef(rfa[[1]])[length(coef(rfa[[1]]))] / 
           coef(rfa[[2]])[length(coef(rfa[[2]]))]
      coefs <- coefs * C

      ## check whether the estimator is non-decreased
      ## NOTE: it only checks for t >= 1 !!!

      deriv.f <- function(t) {
        Re(sapply(t, function(x) {-(coefs*roots) %*% ( 1 / ((x-roots)^2))}))}
      if (any( deriv.f(seq(1, 100, by=0.05)) < 0 )) {
        m <- m - 2
        next
      } else {
        f.mincount <- function(t) {
          sapply(r, function(x) {
            Re(coefs %*% (t / (t - roots))^x)})}
        f.mincount(1)
        valid.estimator <- TRUE
      }
    }
  }
  ## remove M, M.adjust in the future
  if (valid.estimator == TRUE) {
    return(list(FUN=f.mincount, M=m / 2, M.adjust=length(roots), 
                FUN.elements=list(coefs=coefs, roots=roots)))
  } else {
    ## the case m = 0
    f.mincount <- function(t) {
      sapply(r, function(x) sum(n[, 2]))}
    return(list(FUN=f.mincount, M=1, M.adjust=1, 
                FUN.elements=list(coefs=sum(n[, 2]), roots=0)))
  }
}
