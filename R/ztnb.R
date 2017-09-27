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

### initial settings of two parameters size and mu in a negative binomial
### distribution for the numeric optimal searching function optim in R
SIZE.INIT <- 1
MU.INIT <- 0.5

### termination conditions for EM algorithm
TOLERANCE <- 1e-10
ITER.TOLERANCE <- 1e5


### density function of a zero-truncated negative binomial distribution
### size and mu are two parameters for the negative binomial
dztnb <- function(x, size, mu, log = FALSE)
{
  ## the density of x in negative binomial
  p <- dnbinom(x, size = size, mu = mu, log = log)

  ## set zeros in x with zero probability
  if (log == FALSE) {
    p[ which(x == 0) ] <- 0
  } else {
    p[ which(x == 0) ] <- -Inf
  }

  ## the density of non-zero in negative binomial
  q <- 1 - dnbinom(0, size = size, mu = mu)

  ## normalize all non-zero values in negrative binomial
  if (log == FALSE) {
    return( p/q )
  } else {
    return( p - log(q) )
  }
}


### zero-truncated negative loglikelihood
ztnb.minus.loglikelihood <- function(n, size, mu)
{
  prob <- dztnb(n[, 1], size, mu, log = TRUE)

  ## negative loglikelihood
  prob <- -prob
  return( prob %*% n[, 2] )
}


### calculate the negative binomial loglikelihood
### zero.count is the number of unobserved species
nb.loglikelihood <- function(n, zero.count, size, mu)
{
  ## loglikelihood for nonzero counts
  log.prob <- dnbinom(n[, 1], size = size, mu = mu, log = TRUE)
  loglikelihood <- log.prob %*% n[, 2]

  ## add loglikelihood for zero count
  log.zero.prob <- dnbinom(0, size = size, mu = mu, log = TRUE)
  loglikelihood <- loglikelihood + zero.count * log.zero.prob

  return(loglikelihood)
}


### EM algorithm to fit the histogram with a negative binomial distribution
### n is the histogram for observed species
### the number of unobserved species is the missing data
preseqR.ztnb.em <- function(n, size=SIZE.INIT, mu=MU.INIT)
{
  checking.hist(n)

  n[, 1] <- as.numeric(n[, 1])
  n[, 2] <- as.numeric(n[, 2])
  zero.prob <- exp(dnbinom(0, size = size, mu = mu, log = TRUE))

  S <- sum(n[, 2])
  ## estimate the total number of species
  L <- S / ( 1 - zero.prob )

  ## expected the number of zero counts
  zero.counts <- L*zero.prob

  ## estimated mean and variance
  m <- (n[, 1] %*% n[, 2]) / L
  v <- ( (n[, 1] - m)^2 %*% n[, 2] + m^2 * zero.counts )/(L - 1)

  ## target function f
  f <- function(x) {
        return( -nb.loglikelihood(n, zero.counts, size = x, mu = m)/L )
  }

  ## derivative of f
  ## zero.counts is an external variable that are updated by the EM algorithm
  ## CHECK IT!
  gr <- function(x)
  {
    first.term <- ( digamma(x) * zero.counts +
                    digamma(n[, 1] + x) %*% n[, 2] )/L
    second.term <- digamma(x)
    third.term <- log(x) - log(x + m)
    result <- first.term - second.term + third.term
    # f is negative loglikelihood
    return(-result)
  }

  ## estimate size and mu based on first and second moments
  if (v > m) {
    res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B",
           lower = 0.0001, upper = 10000)
  } else {
    res <- optim(size, f, gr, method = "L-BFGS-B",
           lower = 0.0001, upper = 10000)
  }

  ## count the times of iteration
  iter <- as.double(1)

  ## initialize the negative loglikelihood
  loglikelihood.pre <- Inf

  ## zero-truncated loglikelihood
  loglikelihood <- ztnb.minus.loglikelihood(n, res$par, m)

  ## EM algorithm
  while (( loglikelihood.pre - loglikelihood ) / S > TOLERANCE &&
           iter < ITER.TOLERANCE)
  {
    ## update negative loglikelihood
    loglikelihood.pre <- loglikelihood

    ## update parameters
    size <- res$par
    mu <- m

### E-step: estimate the number of unobserved species

    zero.prob <- exp(dnbinom(0, size = size, mu = mu, log = TRUE))
    L <- S / ( 1 - zero.prob )
    zero.counts <- L * zero.prob
    m <- (n[, 1] %*% n[, 2])/L
    v <- ( (n[, 1] - m)^2 %*% n[, 2] + m^2 * zero.counts )/(L - 1)

### M step: estimate the parameters size and mu

    if (v > m) {
      res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B",
             lower = 0.0001, upper = 10000)
    } else {
      res <- optim(size, f, gr, method = "L-BFGS-B",
             lower = 0.0001, upper = 10000)
    }
    iter <- iter + 1
    loglikelihood <- ztnb.minus.loglikelihood(n, res$par, m)
  }
  return(list(size = size, mu = mu, loglik = -loglikelihood.pre))
}


## fitting the negative binoimal distribution to the data by EM algorithm
ztnb.rSAC <- function(n, r=1, size=SIZE.INIT, mu=MU.INIT)
{
  checking.hist(n)

  n[, 2] <- as.numeric(n[, 2])
  S <- sum(n[, 2])

  ## estimate parameters
  opt <- preseqR.ztnb.em(n, size, mu)
  size <- opt$size
  mu <- opt$mu

  ## the probability of a species observed in the initial sample
  p <- 1 - dnbinom(0, size = size, mu = mu)

  ## L is the estimated number of species in total
  L <- S / p

  f.rSAC <- function(t) {
    L * pnbinom(r - 1, size=size, mu=mu*t, lower.tail=FALSE)
  }
  f.rSAC
}
