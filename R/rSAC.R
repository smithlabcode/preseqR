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

## interpolating for r-SAC
## ss, step size
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
  if (l < 1 || ss < 1 || r < 1) {
    return(NULL)
  ## the step size is the sample size
  ## count the number of species observed r or more times
  }  else if (l == 1) {
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
      if (x <= r - 1) {
        return(1)
      } else {
        logp = lchoose(N - x, size - 0:(r-1)) + lchoose(x, 0:(r-1)) - denom
        return(sum(exp(logp)))
      }})
    return(S - p %*% n[, 2])
  }

  ## subsample sizes
  X <- step.size * ( 1:l )

  ## the number of species represented at least r times in a subsample
  yields <- sapply(X, function(x) {
      expect.distinct(n, N, x, initial.distinct, r)})

  result <- matrix(c(X, yields), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'interpolation')

  return(result)
}

## interpolate rSAC based on sampling without replacement
#preseqR.interpolate.rSAC <- function(n, ss, r=1) {
#  checking.hist(n)
#  n[, 2] <- as.double(n[, 2])
#  N <- n[, 1] %*% n[, 2]
#  p <- ss / N
#  i <- dim(n)[1]
#  frequencies <- lapply(1:i, function(x) rbinom(n[x, 2], n[x, 1], p))
#  f <- unlist(frequencies)
#  h <- hist(f, breaks=-1:max(f), plot=FALSE)$counts[-1]
#  matrix(c(which(h != 0), h[which(h != 0)]), byrow = FALSE, ncol=2) 
#}


## coefficients for the power series of E(S_1(t)) / t
## return the first mt terms
discoveryrate.ps <- function(n, mt)
{
  ## transform a histogram into a vector of frequencies
  hist.count <- vector(length=max(n[, 1]), mode="numeric")
  hist.count[n[, 1]] <- n[, 2]

  PS.coeffs <- sum(hist.count)
  if (mt == 1) {
    return(PS.coeffs)
  }

  change.sign <- 0
  for (i in 1:(min(mt-1, length(hist.count)))) {
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


ds.rSAC <- function(n, r=1, mt=20)
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

  ## construct the continued fraction approximation to the power seies
  cf <- ps2cfa(coefs=PS.coeffs, mt=mt)
  rf <- cfa2rf(CF=cf)
  ## the length of cf could be less than mt
  ## even if ps do not have zero terms, coefficients of cf may have
  mt <- length(cf)
  mt <- mt - (mt %% 2)
  valid.estimator <- FALSE
  m <- mt
  while (valid.estimator == FALSE && m >= 2) {

    ## rational function approximants [m / 2 - 1,  m / 2]
    rfa <- rf2rfa(RF=rf, m=m)
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
      f.rSAC <- function(t) {
        sapply(t, function(x) {
          Re(coefs %*% (x / (x - poles))^r)})}
      valid.estimator <- TRUE
    }
  }

  if (valid.estimator == TRUE) {
    return(f.rSAC)
  } else {
    ## the case S1 = S2 where the numbe of species represented exactly once
    ## is 0
    return(function(t) {sapply(t, function(x) return(sum(n[, 2])))})
  }
}


## the bootstrap version of ds.rSAC
## with confidence interval
ds.rSAC.bootstrap <- function(n, r=1, mt=20, times=30, conf=0.95)
{
  n[, 2] <- as.numeric(n[, 2])
  ## individuals in the sample
  N <- n[, 1] %*% n[, 2]

  ## returned function
  f.rSACs <- vector(length=times, mode="list")

  f.bootstrap <- function(n, r, mt) {
    n.bootstrap <- matrix(c(n[, 1], rmultinom(1, sum(n[, 2]), n[, 2])), ncol=2)
    N.bootstrap <- n.bootstrap[, 1] %*% n.bootstrap[, 2]
    N <- n[, 1] %*% n[, 2]
    t.scale <- N / N.bootstrap
    f <- ds.rSAC(n.bootstrap, r=r, mt=mt)
    return(function(t) {f(t * t.scale)})
  }

  while (times > 0) {
    f.rSACs[[times]] <- f.bootstrap(n=n, r=r, mt=mt)
    ## prevent later binding!!!
    f.rSACs[[times]](1)
    times <- times - 1
  }

  estimator <- function(t) {
    result <- sapply(f.rSACs, function(f) f(t))
    if (length(t) == 1) {
      return(median(result))
    } else {
      return(apply(result, FUN=median, MARGIN=1))
    }
  }

  variance <- function(t) {
    result <- sapply(f.rSACs, function(f) f(t))
    if (length(t) == 1) {
      return(var(result))
    } else {
      return(apply(result, FUN=var, MARGIN=1))
    }
  }
  
  se <- function(x) sqrt(variance(x))

  ## prevent later binding!!!
  estimator(1); estimator(1:2)
  variance(1); variance(1:2)
  ## confidence interval using lognormal
  q <- (1 + conf) / 2
  lb <- function(t) {
    C <- exp(qnorm(q) * sqrt(log( 1 + variance(t) / (estimator(t)^2) )))
    return(estimator(t) / C)
  }
  ub <- function(t) {
    C <- exp(qnorm(q) * sqrt(log( 1 + variance(t) / (estimator(t)^2) )))
    return(estimator(t) * C)
  }
  return(list(f=estimator, se=se, lb=lb, ub=ub))
}

## Best practice
preseqR.rSAC <- function(n, r=1, mt=20, size=SIZE.INIT, mu=MU.INIT)
{
  para <- preseqR.ztnb.em(n)
  shape <- para$size
  mu <- para$mu
  ## the population is heterogeneous
  ## because the coefficient of variation is large $1 / sqrt(shape)$
  if (shape <= 1) {
    f.rSAC <- ds.rSAC(n=n, r=r, mt=mt)
  } else {
    ## the population is close to be homogeneous
    ## the ZTNB approach is applied

    ## the probability of a species observed in the initial sample
    p <- 1 - dnbinom(0, size = size, mu = mu)
    ## L is the estimated number of species in total
    L <- sum(as.numeric(n[, 2])) / p
    ## ZTNB estimator
    f.rSAC <- function(t) {
      L * pnbinom(r - 1, size=size, mu=mu*t, lower.tail=FALSE)
    }
  }
  return(f.rSAC)
}


## the bootstrap version of preseqR.rSAC
## with confidence interval
preseqR.rSAC.bootstrap <- function(n, r=1, mt=20, size=SIZE.INIT, 
                                   mu=MU.INIT, times=30, conf=0.95)
{
  checking.hist(n)
  n[, 2] <- as.numeric(n[, 2])

  para <- preseqR.ztnb.em(n)
  shape <- para$size
  mu <- para$mu
  if (shape <= 1) {
    f.bootstrap <- function(n, r, mt, size, mu) {
      n.bootstrap <- matrix(c(n[, 1], rmultinom(1, sum(n[, 2]), n[, 2])), ncol=2)
      N.bootstrap <- n.bootstrap[, 1] %*% n.bootstrap[, 2]
      N <- n[, 1] %*% n[, 2]
      t.scale <- N / N.bootstrap
      f <- ds.rSAC(n.bootstrap, r=r, mt=mt)
      return(function(t) {f(t * t.scale)})
    }
  } else {
    f.bootstrap <- function(n, r, mt, size, mu) {
      n.bootstrap <- matrix(c(n[, 1], rmultinom(1, sum(n[, 2]), n[, 2])), ncol=2)
      N.bootstrap <- n.bootstrap[, 1] %*% n.bootstrap[, 2]
      N <- n[, 1] %*% n[, 2]
      t.scale <- N / N.bootstrap
      f <- ztnb.rSAC(n.bootstrap, r=r, size=size, mu=mu)
      return(function(t) {f(t * t.scale)})
    }
  }

  ## returned function
  f.rSACs <- vector(length=times, mode="list")

  while (times > 0) {
    f.rSACs[[times]] <- f.bootstrap(n=n, r=r, mt=mt, size=size, mu=mu)
    ## prevent later binding!!!
    f.rSACs[[times]](1)
    times <- times - 1
  }

  estimator <- function(t) {
    result <- sapply(f.rSACs, function(f) f(t))
    if (length(t) == 1) {
      return(median(result))
    } else {
      return(apply(result, FUN=median, MARGIN=1))
    }
  }

  variance <- function(t) {
    result <- sapply(f.rSACs, function(f) f(t))
    if (length(t) == 1) {
      return(var(result))
    } else {
      return(apply(result, FUN=var, MARGIN=1))
    }
  }

  se <- function(x) sqrt(variance(x))

  ## prevent later binding!!!
  estimator(1); estimator(1:2)
  variance(1); variance(1:2)
  ## confidence interval using lognormal
  q <- (1 + conf) / 2
  lb <- function(t) {
    C <- exp(qnorm(q) * sqrt(log( 1 + variance(t) / (estimator(t)^2) )))
    ## if var and estimates are 0
    C[which(!is.finite(C))] = 1
    return(estimator(t) / C)
  }
  ub <- function(t) {
    C <- exp(qnorm(q) * sqrt(log( 1 + variance(t) / (estimator(t)^2) )))
    ## if var and estimates are 0
    C[which(!is.finite(C))] = 1
    return(estimator(t) * C)
  }
  return(list(f=estimator, se=se, lb=lb, ub=ub))
}
