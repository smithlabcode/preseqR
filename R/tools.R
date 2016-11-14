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


### checking the input histogram in an appropariat format
checking.hist <- function(n)
{
  if (ncol(n)!=2 || is.numeric(n[,1])==FALSE || is.numeric(n[,2])==FALSE) {
    stop("Input must be a two-column matrix")
  }
  ## the first column is the frequencies of observed items
  freq <- n[, 1]

  ## the second column is the number of observed distinct items for each
  ## frequency
  number.items <- n[, 2]

  ## check whether frequencies are at least one and the histogram is sorted
  for (i in 1:length(freq))
    if (freq[i] <= 0 || freq[i] != floor(freq[i])) {
      stop("The first column must be positive integers!")
    } else if (number.items[i] < 0) {
      stop("The second column must be non negative")
    }
    else {
        if (i > 1 && freq[i - 1] >= freq[i])
          stop("The first column is not sorted in the ascending order")
    }

  return(n)
}


## sampling without replacement
## n frequencies counts
nonreplace.sampling <- function(size, n)
{
  ## make sure frequencies are integers
  n[, 2] <- floor(n[, 2])
  ## the number of distinct items
  distinct <- sum(n[, 2])

  ## identifier for each distinct item
  ind <- 1:distinct

  ## the size of each read in the library
  N <- rep(n[, 1], n[, 2])

  ## construct a sample space X 
  ## the whole library represents by its indexes. If a read presents t
  ## times in the library, its indexes presents t times in X
  X <- rep(ind, N)

  return(sample(X, size, replace = FALSE))
}


## sampling without replacement
## input frequencies counts; output subsample as a frequencies counts
preseqR.nonreplace.sampling <- function(size, n)
{
  ## check the input histogram file
  checking.hist(n)
  ## sub sampling
  X <- nonreplace.sampling(size, n)
  ## record the freq of each sampled species
  freq.counts <- hist(X, breaks=0:max(X), plot=FALSE)$count
  ## frequencies counts; frequency 0 excluded
  T <- hist(freq.counts, breaks=-1:max(freq.counts), plot=FALSE)$counts[-1]
  matrix(c(which(T != 0), T[which(T != 0)]), byrow = FALSE, ncol=2)
}


lchoose <- function(N, k) {
  result <- vector(length=max(length(N), length(k)), mode="numeric")
  index <- which(N - k + 1 > 0)
  if (length(index) == 0) {
    result[] <- -Inf }
  else {
    result[index] <- (lgamma(N + 1) - lgamma(k + 1))[index] - lgamma((N - k + 1)[index])
    result[-index] <- -Inf
  }
  result
}
