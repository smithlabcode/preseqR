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

## Zero Truncated Poisson (ZTP) estimator
## Cohen, A. Clifford. "Estimating the parameter in a conditional Poisson 
## distribution." Biometrics 16, no. 2 (1960): 203-211.
ztp.rSAC <- function(n, r=1) 
{
  ## mle based on zero-truncated Poisson
  n[, 2] <- as.numeric(n[, 2])
  C <- n[, 1] %*% n[, 2] / sum(n[, 2])
  f <- function(x) {x / (1 - exp(-x))}
  result <- uniroot(function(x) f(x) - C, c(0.001, 1e9), tol = 0.0001,
                    extendInt="upX")
  lambda <- result$root
  L <- sum(n[, 2]) / (1 - ppois(0, lambda))

  ## estimator 
  function(x) {
    L * ppois(q=r - 1, lambda=lambda * x, lower.tail=FALSE)}
}


## BBC estimator
## Boneh, Shahar, Arnon Boneh, and Richard J. Caron. "Estimating the prediction
## function and the number of unseen species in sampling with replacement." 
## Journal of the American Statistical Association 93, no. 441 (1998): 372-379.
bbc.rSAC <- function(n, r=1) {
  n[, 2] <- as.numeric(n[, 2])
  S <- sum(n[, 2])
  ## BBC estimator without bias correction
  tmp <- function(t) {sapply(t, function(x) {
            n[, 2] %*% (exp(-n[, 1]) - ppois(r-1, n[, 1] * x)) + S})}
  ## bias correction
  index.f1 <- which(n[, 1] == 1)
  f1 <- n[index.f1, 2]
  U0 <- n[, 2] %*% exp(-(n[, 1]))
  ## if satisfy the condition, adjust the bias
  ## otherwise return the estimator without correction
  if (length(index.f1) == 1 && f1 > U0) {
    result <- uniroot(function(x) x*(1 - exp(-f1 / x)) - U0, c(0.001, 1e9),
                      tol=0.0001, extendInt="upX")
    U <- result$root
    f.rSAC <- function(t) {tmp(t) + U*(exp(-(f1 / U)) - ppois(r-1, f1 * t / U))}
  } else {
    f.rSAC <- tmp
  }
  f.rSAC
}


## CS estimator
## Chao, Anne, and Tsung-Jen Shen. "Nonparametric prediction in species 
## sampling." Journal of agricultural, biological, and environmental 
## statistics 9, no. 3 (2004): 253-269.
cs.rSAC <- function(n, r=1, k=10) {
  n[, 2] <- as.numeric(n[, 2])
  S <- sum(n[, 2])
  index.f1 <- which(n[, 1] == 1)
  ## no species observed once
  if (length(index.f1) != 1)
    return(NULL)
  f1 <- n[index.f1, 2]
  ## rare species if frequency no more than k
  index.rare <- which(n[, 1] <= k)
  S.rare <- sum(n[index.rare, 2])
  C.rare <- 1 - f1 / (n[index.rare, 1] %*% n[index.rare, 2])
  ## estimate parameters
  gamma.rare <- max(S.rare / C.rare * 
    ((n[index.rare, 1] * (n[index.rare, 1] - 1)) %*% n[index.rare, 2]) /
    (n[index.rare, 1] %*% n[index.rare, 2])^2 - 1, 0)
  f0 = S.rare / C.rare + f1 / C.rare * gamma.rare - S.rare
  ## estimator
  function(t) {f0 + S  - f0 * ppois(r-1, f1 * t / f0) * exp(f1/f0) }
}


## The following estimator is based on the logseries estimator by Fisher
## Fisher, R. A., A. Steven Corbet, and C. B. Williams. "The Relation Between
## the Number of Species and the Number of Individuals in a Random Sample of an
## Animal Population." Journal of Animal Ecology 12, no. 1 (1943): 42-58.
## doi:10.2307/1411.

## parametric for the logseries
fisher.alpha <- function(n) {
  n[, 2] <- as.numeric(n[, 2])
  N <- n[, 1] %*% n[, 2]
  S <- sum(n[, 2])
  result <- uniroot(function(x) (exp(x) - 1) / x - N / S, c(0.001, 1e9), 
                    tol=0.0001, extendInt="upX")
  alpha <- S / result$root
  return(alpha)
}

## logseries estimator for the number of species represented at least r times
fisher.rSAC <- function(n, r=1) {
  n[, 2] <- as.numeric(n[, 2])
  alpha <- fisher.alpha(n)
  N <- n[, 1] %*% n[, 2]
  f.rSAC <- function(t) {sapply(t, function(x) alpha *
    integrate(function(z) (z^(r-1) / (1 - z)), lower=0, 
                upper=N*x / (N*x + alpha))$value)}
  return(f.rSAC)
}


## fit a Poisson distribution
pois.rSAC <- function(n, L, r=1) {
    lambda <- n[, 1] %*% n[, 2] / L
    f.rSAC <- function(t) {
      L * ppois(q=r - 1, lambda=lambda * t, lower.tail=FALSE)
    }
    return(f.rSAC)
}

## estimate the parameter for negative binomial distribuiton
nb.fitting <- function(n, L, size=SIZE.INIT)
{
  n[, 2] <- as.numeric(n[, 2])

  ## the number of unobservations
  zero.counts <- L - sum(n[, 2])

  ## estimated mean and variance
  m <- (n[, 1] %*% n[, 2]) / L
  v <- ( (n[, 1] - m)^2 %*% n[, 2] + m^2 * zero.counts )/(L - 1)

  ## target function f
  f <- function(x) {
        return( -nb.loglikelihood(n, zero.counts, size = x, mu = m)/L )
  }

  ## derivative of f
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

  ## optimization
  if (v > m) {
    res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B",
           lower = 0.00001, upper = 100000)
  } else {
    res <- optim(size, f, gr, method = "L-BFGS-B",
           lower = 0.00001, upper = 100000)
  }

  loglikelihood <- nb.loglikelihood(n, zero.counts, size=res$par, mu=m)
  ## update parameters
  size <- res$par
  mu <- m

  return(list(size = size, mu = mu, loglik = -loglikelihood))
}


## fitting the negative binoimal distribution to the data
## L is the total number of species
nb.rSAC <- function(n, L, r=1, size=SIZE.INIT)
{
  n[, 2] <- as.numeric(n[, 2])

  ## estimate parameters
  opt <- nb.fitting(n, L, size=size)
  size <- opt$size
  mu <- opt$mu

  f.rSAC <- function(t) {
    L * pnbinom(r-1, size=size, mu=mu*t, lower.tail=FALSE)
  }
  return(f.rSAC)
}
