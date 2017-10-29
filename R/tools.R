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


### check the input histogram in an appropriate format
checking.hist <- function(n) 
{
  if (ncol(n)!=2 || is.numeric(n[,1])==FALSE || is.numeric(n[,2])==FALSE) {
    stop("Input must be a two-column matrix")
  }
  ## the first column is the frequency i
  ## the second column is the number of species represented i times 
  ## in the sample
  freq <- n[, 1]
  num <- n[, 2]

  ## check whether frequencies are at least one and the histogram is sorted
  ## based on frequencies
  for (i in 1:length(freq))
    if (freq[i] <= 0 || freq[i] != floor(freq[i])) {
      stop("The first column must be positive integers!")
    } else if (num[i] < 0) {
      stop("The second column must be non negative")
    }
    else {
        if (i > 1 && freq[i - 1] >= freq[i])
          stop("The first column is not sorted in the ascending order")
    }

  return(n)
}


### check determinants of matrix M_{m-1,m-1},M_{m-1,m},M_{m,m-1},M_{m,m}
## OBSOLETE
checking.matrix.det <- function(n, m) 
{
  ps <- discoveryrate.ps(n, mt=2*m + 1)
  ps <- c(0, ps)
  matrix.dets <- vector(length=4, mode="numeric")
  count <- 1
  for (i in (m-1):m)
    for (j in (m-1):m) {
      pade.matrix <- sapply(1:(i+1), function(x) {
      start <- j + 1 + x
      end <- j - i + 1 + x
      indexes <- seq(start, end, -1)
      ps[indexes]})

      matrix.dets[count] <- det(pade.matrix / max(abs(pade.matrix)))
      count <- count + 1
    }
  matrix.dets
}


## sampling from a histogram without replacement
nonreplace.sampling <- function(n, size)
{
  ## make sure the number of species are integers
  n[, 2] <- floor(n[, 2])
  ## the number of species in total
  S <- sum(n[, 2])

  ## identifier for each species
  ind <- 1:S

  ## construct a sample space X composed by all individuals
  N <- rep(n[, 1], n[, 2])
  X <- rep(ind, N)

  return(sample(X, size, replace = FALSE))
}


## sampling from a histogram without replacement
## both input and output are histograms
preseqR.nonreplace.sampling <- function(n, size)
{
  ## check the input histogram file
  checking.hist(n)
  ## subsampling
  X <- nonreplace.sampling(n, size)
  ## record the freq of each sampled species
  freq <- hist(X, breaks=0:max(X), plot=FALSE)$count
  ## the frequency of each frequency
  T <- hist(freq, breaks=-1:max(freq), plot=FALSE)$counts[-1]
  ## histogram
  matrix(c(which(T != 0), T[which(T != 0)]), byrow = FALSE, ncol=2)
}

## OBSOLETE
## log of N choose k
#lchoose <- function(N, k) 
#{
#  result <- vector(length=max(length(N), length(k)), mode="numeric")
#  index <- which(N - k + 1 > 0)
#  if (length(index) == 0) {
#    result[] <- -Inf }
#  else {
#    result[index] <- (lgamma(N + 1) - lgamma(k + 1))[index] - 
#                     lgamma((N - k + 1)[index])
#    result[-index] <- -Inf
#  }
#  result
#}
