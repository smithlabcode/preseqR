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

## get the expectation of a distribution through Law of Large Number
get.expectation <- function(FUN) {
  ## magic number
  N <- 1e6
  return(sum(FUN(N)) / N)
}


## generate histograms based on the mixture of Poisson distributions model
## L, the total number of species in a population
## FUN, the RNG used to generate the relative abundance for each species 
preseqR.simu.hist <- function(L=1e8, N, FUN) {
  if (L <= 0) {
    write("L has to be a positive number", stderr())
    return(NULL)
  }
  L <- as.integer(L)
  ## save the poisson parameters for each individual in the population
  lambda <- FUN(L)
  ## S saves samples
  S <- rmultinom(n=1, size=N, prob=lambda)
  t <- table(S)
  ## histogram
  h <- matrix(c(as.integer(names(t)), as.numeric(t)), ncol=2)
  ## exclude 0 counts
  if (h[1, 1] == 0) {
    h = h[2:dim(h)[1], ]
  }
  return(h)
}


## simulate a sample in a population
## L, the total number of species in a population
## FUN, the RNG used to generate the relative abundance for each species 
## t, a vector representing relative sample sizes 
## output 
## 1. the histogram of the sample
## 2. the probabilty of a species represented at least r times in a random
##    sample that is t times the size of the given sample
preseqR.simu.sc <- function(L=1e8, N, t, r, FUN) {
  if (L > 1e8) {
    L <- 1e8
  } else if (L <= 0) {
    write("L has to be a positive number", stderr())
    return(NULL)
  }
  L <- as.integer(L)
  E <- get.expectation(FUN)
  t0 <- N / (L * E)
  ## save the Poisson rate for each individual in the population
  lambda <- FUN(L)
  ## the number of individuals for each species
  counts <- rpois(L, lambda * t0)
  hist.count <- vector(length = max(counts), mode = "numeric")
  for (freq in counts) {
    if (freq > 0) {
      hist.count[freq] <- hist.count[freq] + 1
    }
  }
  nonzero.index <- which(hist.count != 0)
  nonzero <- hist.count[nonzero.index]
  n <- matrix(c(nonzero.index, nonzero), ncol = 2)
  coverage <- sapply(t, function(x) {
                counts <- rpois(L, lambda * t0 * x)
                sum(lambda[which(counts >= r)]) / sum(lambda)})
  return(list(histogram=n, sc=coverage))
}


## simulating an interpolation curve
## FUN should always generate positive number
## OBSOLATE
preseqR.simu.interpolate <- function(L=1e7, ss, max.size, r, FUN) {
  ## too much memory if L is too large
  if (L > 1e8) {
    L <- 1e8
  } else if (L <= 0) {
    write("L has to be a positive number", stderr())
    return(NULL)
  }
  L <- as.integer(L)
  E <- get.expectation(FUN)
  t <- ss / (L * E)
  ## save the poisson parameters for each individual in the population
  ## assume all species in the library have positive probability to be sampled
  lambda <- FUN(L)
  
  ## the interpolation curve
  N <- floor(max.size / ss)
  S <- vector(length = L, mode = "numeric")
  points <- c(0, 0)
  hists <- list()
  for (i in 1:N) {
    res <- rpois(L, lambda * t)
    S <- S + res

    hist.count <- vector(length = max(S), mode = "numeric")
    for (freq in S) {
      if (freq > 0) {
        hist.count[freq] <- hist.count[freq] + 1
      }
    }
    ## add the new point into the interpolation curve
    point <- c(1:length(hist.count) %*% hist.count, 
               sum(hist.count[r:length(hist.count)]))
    points <- rbind(points, point)
    ## add the new generated histogram
    nonzero.index <- which(hist.count != 0)
    nonzero <- hist.count[nonzero.index]
    n <- matrix(c(nonzero.index, nonzero), ncol = 2)
    hists <- c(hists, list(n))
  }
  list(histograms = hists, interpolation.curve = unname(points))
}
