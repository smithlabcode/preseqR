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

### Generate histograms through simulation
### The basic assumption is that the number of each individual sampled follows
### a poisson process.

### get the expectation of a distribution through Law of Large Number
get.expectation <- function(FUN) {
  ## magic number
  N <- 1e6
  return(sum(FUN(N)) / N)
}

### generate the histogram based on the simulation
### L is the total number of species in a population
### t is the relative sample size
### FUN is an RNG. It must take one argument as the number of random numbers
### generated and return the number of positive random numbers
### FUN can be defined by users
preseqR.simu.hist <- function(L=1e8, t, FUN) {
  if (L > 1e8) {
    L <- 1e8
  } else if (L <= 0) {
    write("L has to be a positive number", stderr())
    return(NULL)
  }
  L <- as.integer(L)
  ## save the poisson parameters for each individual in the population
  lambda <- FUN(L)
  ## S saves samples
  S <- rpois(L, lambda * t)
  hist.count <- vector(length = max(S), mode = "numeric")
  for (freq in S) {
    if (freq > 0) {
      hist.count[freq] <- hist.count[freq] + 1
    }
  }
  nonzero.index <- which(hist.count != 0)
  nonzero <- hist.count[nonzero.index]
  matrix(c(nonzero.index, nonzero), ncol = 2)
}

### simulating an interpolation curve
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
  ## assume all items in the library have positive probability to be sampled
  ## FUN should always generate positive number
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
