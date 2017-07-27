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

## Two roots are the same if the difference is less than the PRECISION
PRECISION <- 1e-3

### interpolating for species accumulation curve with minimum count
### ss step size
### n two-column histogram
### r minimum count
preseqR.interpolate.mincount <- function(ss, n, r=1)
{
  checking.hist(n)

  n[, 2] <- as.numeric(n[, 2])
  ## total individuals captured
  total.sample <- n[, 1] %*% n[, 2]
  N <- total.sample

  ## total species
  initial.distinct <- sum(n[, 2])
  step.size <- as.double(ss)

  ## l is the number of sampled points for interpolation
  l <- N / step.size

  ## if the sample size is larger than the size of experiment or 
  ## the step size is too small, return NULL
  if (l < 1 || ss < 1 || r < 1)
    return()
  ## if the sample size is the size of the experiment
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
  ## N total individuals
  ## S the number of species
  ## size the size of the subsample
  expect.distinct <- function(n, N, size, S, r) {
    denom <- lchoose(N, size)
    p <- sapply(n[, 1], function(x) {
       sum(exp(lchoose(N - x, size - 0:(r-1)) + lchoose(x, 0:(r-1)) - denom))})
    return(S - p %*% n[, 2])
  }

  ## sample sizes
  x <- step.size * ( 1:l )

  ## calculate the number of distinct reads based on each sample size
  yield.estimates <- sapply(x, function(x) {
      expect.distinct(n, N, x, initial.distinct, r)})

  ## put size and yield together into a matrix
  result <- matrix(c(x, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'interpolation')

  return(result)
}


### power series based on count frequencies starting from frequency j
### when j = 1, it is the power series expansion of E(S_1(t)) / t at t = 1
### the maximum number of terms
generating.ps <- function(n, mt, j=1) {
  if (j >= max(n[, 1])) return(NULL)
  ## transform a histogram into a vector of frequencies
  hist.count <- vector(length=max(n[, 1]), mode="numeric")
  hist.count[n[, 1]] <- n[, 2]

  ## shift to required count frequencies
  hist.count <- hist.count[j: length(hist.count)]

  PS.coeffs <- sum(hist.count)
  change.sign <- 0

  ## preserve extra precision mt+1
  for (i in 1:(min(mt+1, length(hist.count)))) {
    PS.coeffs <- c(PS.coeffs, 
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


## species accum curves based on parital fraction expansion
## the function is used whenever count frequency 1 is unavaible or the sample
## size is saturated
## using count frequencies starting from a given count frequency instead of 1
## when start.freq = 1, it is identical to the function preseqR.pf.mincount
## CHAO: save for a rainy day
general.ds.mincount <- function(n, r=1, mt=100, start.freq=1)
{
  # check the input format of the histogram
  checking.hist(n)

  n[, 2] <- as.numeric(n[, 2])
  ## constructing the power series
  PS.coeffs <- generating.ps(n, j=start.freq, mt=mt)

  if (is.null(PS.coeffs)) {
    write("the size of the initial experiment is insufficient", stderr())
    return(NULL)
  }

  ## constrain the continued fraction approximation with even degree 
  ## asymptotically ~ C / t
  mt <- min(mt, length(PS.coeffs))
  PS.coeffs <- PS.coeffs[ 1:mt ]

  ## check whether sample size is sufficient
  if (mt < 2)
  {
    m <- paste("max count before zero is les than min required count (2)",
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
  ## select rational function approximants [M-1/M] m=2M
  ## asymptotically ~ C / t
  mt <- mt - (mt %% 2)
  valid.estimator <- FALSE
  m <- mt
  while (valid.estimator == FALSE) {

    rfa <- rf2rfa(RF=rf, m=m)
    ## solving roots
    numer.roots <- solve(rfa[[1]])
    denom.roots <- solve(rfa[[2]])
    ## seperating roots by their real parts
    numer.roots.neg <- numer.roots[which(Re(numer.roots) < 0)]
    numer.roots.pos <- numer.roots[which(Re(numer.roots) >= 0)]
    denom.roots.neg <- denom.roots[which(Re(denom.roots) < 0)]
    denom.roots.pos <- denom.roots[which(Re(denom.roots) >= 0)]

    ## record roots in the numerator that are significantly similar to
    ## roots in the denominator
    tmp.roots <- c()

    ## simplify the rational function approximation
    ## two roots are same if the difference is less than the 
    ## predefined PRECISION
    if (length(numer.roots.pos) > 0) {
      for (i in 1:length(numer.roots.pos)) {
        if (length(denom.roots.pos) > 0) {
          d <- Mod(denom.roots.pos - numer.roots.pos[i])
          for (j in 1:length(d)) {
            if (d[j] < PRECISION) {
              denom.roots.pos <- denom.roots.pos[-j]
              tmp.roots <- c(tmp.roots, numer.roots.pos[i])
              break
            }
          }
        }
      }
    }

    ## roots in simplified RFA
    numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
    denom.roots <- c(denom.roots.neg, denom.roots.pos)

    ## convert roots from t - 1 to t
    roots <- denom.roots + 1
    ## pacman rule checking
    if (length(which(roots == 0)) || length(which(Re(roots) > 0))) {
      m <- m - 2
      next
    } else {
      poly.numer <- as.function(poly.from.roots(numer.roots))
      l <- length(denom.roots)
      ## treat polynomials in the rational function to be monic
      ## the difference to the original RFA is a multiplier C

      ## c_i in the estimator
      coef <- sapply(1:l, function(x) {
        poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
      ## calculate the constant C
      C <- coef(rfa[[1]])[length(coef(rfa[[1]]))] / 
           coef(rfa[[2]])[length(coef(rfa[[2]]))]
      ## species accum curves with minimum count r
      ## using parital fraction expansion
      denom.roots <- denom.roots + 1
      coef <- coef * C
      ## modify the coefficients
      coef <- coef * (1 - denom.roots)^(start.freq - 1)
      ## check whether the estimator is non-decreased                             
      deriv.f <- function(t) {
        Re(sapply(t, function(x) {-(coef*denom.roots) %*% ( 1 / ((x-denom.roots)^2))}))} 
      if (length(which( deriv.f(seq(0.05, 100, by=0.05)) < 0 ) != 0)) {
        m <- m - 2
        next
      }
      f.mincount <- function(t) {
        sapply(r, function(x) {
          Re(coef %*% (t / (t - denom.roots))^x)})}
      f.mincount(1)
      valid.estimator <- TRUE
    }
  }
  ## remove M, M.adjust in the future
  list(FUN=f.mincount, M=m / 2, M.adjust=length(denom.roots), FUN.elements=list(coef=coef, roots=denom.roots))
}

## nonparametric approach Deng & Smith 2016
ds.mincount.bootstrap <- function(n, r=1, mt=100, times=100)
{
  n[, 2] <- as.numeric(n[, 2])
  ## total individuals
  total <- n[, 1] %*% n[, 2]

  ## returned function
  f.mincount <- vector(length=times, mode="list")

  ds.estimator <- function(n, r, mt, t.scale) {
    f <- ds.mincount(n, r=r, mt=mt)
    if (f$M == 1) {
      f <- ztnb.mincount(n, r=r)
      function(t) {f(t * t.scale)}
    } else {
      function(t) {f$FUN(t * t.scale)}
    }
  }

  while (times > 0) {
    n.bootstrap <- matrix(c(n[, 1], rmultinom(1, sum(n[, 2]), n[, 2])), ncol=2)
    total.bootstrap <- n.bootstrap[, 1] %*% n.bootstrap[, 2]
    t.scale <- total / total.bootstrap
    f <-  ds.estimator(n.bootstrap, r=r, mt=mt, t.scale=t.scale) 

    f.mincount[[times]] <- f
    ## prevent later binding!!!
    f.mincount[[times]](1)
    times <- times - 1
  }
  f.estimator <- ds.mincount(n=n, r=r, mt=mt)
  if (length(r) == 1) {
    median.estimators <- function(t) {median( sapply(f.mincount, function(x) x(t)) )}
    var.estimator <- function(t) {var( sapply(f.mincount, function(x) x(t)) )}
  } else {
    median.estimators <- function(t) {apply(sapply(f.mincount, function(x) x(t)), FUN=median, MARGIN=1)}
    var.estimator <- function(t) {apply(sapply(f.mincount, function(x) x(t)), FUN=var, MARGIN=1)}
  }
  ## prevent later binding!!!
  f.estimator$FUN(1); median.estimators(1); var.estimator(1)
  return(list(FUN.nobootstrap=f.estimator, FUN.bootstrap=median.estimators, var=var.estimator))
}

ds.mincount <- function(n, r=1, mt=100)
{
  checking.hist(n)

  n[, 2] <- as.numeric(n[, 2])

  ## constructing the power series
  PS.coeffs <- generating.ps(n, mt=mt, j=1)

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
  ## select rational function approximants [M-1/M] m=2M
  ## asymptotically ~ C / t
  mt <- mt - (mt %% 2)
  valid.estimator <- FALSE
  m <- mt
  while (valid.estimator == FALSE) {

    rfa <- rf2rfa(RF=rf, m=m)
    ## solving roots
    numer.roots <- solve(rfa[[1]])
    denom.roots <- solve(rfa[[2]])

    ## finite
    if (any(!is.finite(c(numer.roots, denom.roots)))) {
      m = m - 2
      next;
    }

    ## seperating roots by their real parts
    numer.roots.neg <- numer.roots[which(Re(numer.roots) < 0)]
    numer.roots.pos <- numer.roots[which(Re(numer.roots) >= 0)]
    denom.roots.neg <- denom.roots[which(Re(denom.roots) < 0)]
    denom.roots.pos <- denom.roots[which(Re(denom.roots) >= 0)]

    ## record roots in the numerator that are significantly similar to
    ## roots in the denominator
    tmp.roots <- c()

    ## simplify the rational function approximation
    ## two roots are same if the difference is less than the 
    ## predefined PRECISION
    if (length(numer.roots.pos) > 0) {
      for (i in 1:length(numer.roots.pos)) {
        if (length(denom.roots.pos) > 0) {
          d <- Mod(denom.roots.pos - numer.roots.pos[i])
          for (j in 1:length(d)) {
            if (d[j] < PRECISION) {
              denom.roots.pos <- denom.roots.pos[-j]
              tmp.roots <- c(tmp.roots, numer.roots.pos[i])
              break
            }
          }
        }
      }
    }

    ## roots in simplified RFA
    numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
    denom.roots <- c(denom.roots.neg, denom.roots.pos)

    ## convert roots from t - 1 to t
    roots <- denom.roots + 1
    ## pacman rule checking
    if (length(which(roots == 0)) || length(which(Re(roots) > 0))) {
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
  ## remove M, M.adjust in the future)
  if (valid.estimator == TRUE) {
    return(list(FUN=f.mincount, M=m / 2, M.adjust=length(roots), FUN.elements=list(coefs=coefs, roots=roots)))
  } else {
    ## the case m = 0
    f.mincount <- function(t) {
      sapply(r, function(x) sum(n[, 2]))}
    return(list(FUN=f.mincount, M=1, M.adjust=1, FUN.elements=list(coefs=sum(n[, 2]), roots=0)))
  }
}
