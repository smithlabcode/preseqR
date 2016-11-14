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


pois.mincount <- function(n, L, r=1) {
    lambda <- n[, 1] %*% n[, 2] / L
    f.mincount <- function(t) {
      L * ppois(q=r - 1, lambda=lambda * t, lower.tail=FALSE)
    }
    f.mincount(1); f.mincount
}


nb.fitting <- function(n, L, size=SIZE.INIT)
{
  n[, 2] <- as.numeric(n[, 2])

  ## the number of unobservations
  zero.items <- L - sum(n[, 2])

  ## estimated mean and variance
  m <- (n[, 1] %*% n[, 2]) / L
  v <- ( (n[, 1] - m)^2 %*% n[, 2] + m^2 * zero.items )/(L - 1)

  ## target function f
  f <- function(x) {
        return( -nb.loglikelihood(n, zero.items, size = x, mu = m)/L )
  }

  ## derivative of f
  gr <- function(x)
  {
    first.term <- ( digamma(x) * zero.items +
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
           lower = 0.00001, upper = 100000)
  } else {
    res <- optim(size, f, gr, method = "L-BFGS-B",
           lower = 0.00001, upper = 100000)
  }

  loglikelihood <- nb.loglikelihood(n, zero.items, size=res$par, mu=m)
  ## update parameters
  size <- res$par
  mu <- m

  return(list(size = size, mu = mu, loglik = -loglikelihood))
}


## fitting the negative binoimal distribution to the data
## ss is the step.size
## max.extrapoltion is the maximum value for extrapolation
## r is a vector of frequencies
## L is the total number of species
nb.mincount <- function(n, L, r=1, size=SIZE.INIT)
{
  n[, 2] <- as.numeric(n[, 2])

  ## estimate parameters
  opt <- nb.fitting(n, L, size=size)
  size <- opt$size
  mu <- opt$mu

  f.mincount <- function(t) {
    L * pnbinom(r-1, size=size, mu=mu*t, lower.tail=FALSE)
  }
  f.mincount(1); f.mincount
}


# write out the information about the experiment and the number of reads needs
# to be sequenced
preseqR.depthseq <- function(n.reads, FUN=NULL, LIB.FUN=NULL, L=NULL, rho=0.85, r=8, uniq=TRUE)
{
  checking.hist(n.reads)

  if (is.null(FUN)) {
    cat(paste("No estimator for the number of nucleotides with at least ", r, 
        " aligned reads as a function of sequencing effort\n", sep=""))
    return(NULL)
  }
  if (is.null(LIB.FUN)) {
    cat("No estimator for the library complexity curve\n") 
    return(NULL)
  }
  if (is.null(L)) {
    cat("The length of the targeted regions is not specified\n")
    return(NULL)
  }
  N <- n.reads[, 1] %*% n.reads[, 2]
  cat(paste("The number of reads is ", N, " in the initial experiment\n", sep=""))
  if (FUN(10000) <= L * rho) {
    cat(paste("More than 10000 times the size of the initial experiment is",
              " needed to achieve the required standard\n", sep=""))
    return(NULL)
  }
  t0 <- uniroot(function(x) {FUN(x) - L * rho}, interval=c(0, 10000), tol=0.00001)$root
  ## whether remove the duplicates when counting coverage depth
  ## Counting with duplicates
  if (uniq==FALSE) {
    reads.total <- ceiling(N * t0 / 1000000.0)
    cat(paste("In order to attain ", rho*100, "% of the targeted regions of length ", 
        L, " with ", r, "X or greater coverage depth\n", sep=""))
    cat(paste("A total of ", reads.total, " million reads are needed in the full experiment\n", sep=""))
    return(reads.total)
  }
  ## duplicates are removed
  N.uniq <- sum(n.reads[, 2])
  if (LIB.FUN(10000) <= N.uniq * t0) {
    cat(paste("More than 10000 times the size of the initial experiment is",
              " needed to achieve the required standard\n", sep=""))
    return(NULL)
  }
  t1 <- uniroot(function(x) {LIB.FUN(x) - N.uniq * t0}, interval=c(0, 10000), tol=0.00001)$root
  reads.total <- ceiling(N * t1 / 1000000.0)
  cat(paste("In order to attain ", rho*100, "% of the targeted regions of length ", 
      L, " with ", r, "X or greater coverage depth\n", sep=""))
  cat(paste("A total of ", reads.total, " million reads are needed in the full experiment\n", sep=""))
  return(reads.total)
}
