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
### sample.size is the size of the sample
### FUN is an RNG. It must take one argument as the number of random numbers
### generated and return the number of positive random numbers
### FUN can be defined by users
preseqR.simu.hist <- function(L=1e8, size, FUN) {
  if (L > 1e8) {
    L <- 1e8
  } else if (L <= 0) {
    write("L has to be a positive number", stderr())
    return(NULL)
  }
  L <- as.integer(L)
  E <- get.expectation(FUN)
  t <- size / (L * E)
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
  hist <- matrix(c(nonzero.index, nonzero), ncol = 2)
  colnames(hist) <- c("n", "n_j")
  return(hist) 
}

### simulating an interpolation curve
preseqR.simu.interpolate <- function(L=1e7, ss, max.size, k, FUN) {
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
    point <- c(1:length(hist.count) %*% hist.count, sum(hist.count[k:length(hist.count)]))
    points <- rbind(points, point)
    ## add the new generated histogram
    nonzero.index <- which(hist.count != 0)
    nonzero <- hist.count[nonzero.index]
    hist <- matrix(c(nonzero.index, nonzero), ncol = 2)
    colnames(hist) <- c("n", "n_j")
    hists <- c(hists, list(hist))
  }
  res <- list(histograms = hists, interpolation.curve = unname(points))
  return(res)
}

### Construct the histogram given the size of the initial sample,
### the underlining Gamma-Poisson distribution and the total number of distinct
### species. Truncated the frequency when it is smaller than one and all frequencies
### are the integeral parts of calculated expectations.
simu.hist <- function(L, alpha, beta, n) {
  ## coefficients of the constructed power series
  PS <- c()
  ## the ratio between the size of the sample to the expected size of an initial sample
  t <- n / (L * alpha * beta)
  ## set a upperbound for the number of items in the power series
  mt <- 100

  i = 1
  while (i < mt) {
    co.eff = L * dnbinom(x=i, size=alpha, prob=1/(1 + t*beta))
    if (floor(co.eff) < 1) {
      break
    } else {
      PS <- c(PS, floor(co.eff))
    }
    i <- i + 1
  }
  l <- length(PS)
  if (l > 0) {
    hist <- matrix(c(1:l, PS), byrow=FALSE, ncol=2)
    return(hist)
  } else {
   return(NULL)
  }
}
