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
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#


## predict the fraction of k-mers represented at least r times in the sample
kmer.frac <- function(n, r=2, mt=20) {
  return(preseqR.sample.cov(n=n, r=r-1, mt=mt))
}


## the fraction of k-mers represented at least r times as a function of 
## sample sizes
kmer.frac.curve <- function(n, k, read.len, seq, r=2, mt=20) {
  f <- kmer.frac(n, r=r, mt=mt)
  if (is.null(f))
    return(NULL)
  n[, 2] <- as.numeric(n[, 2])
  N <- n[, 1] %*% n[, 2]
  ## average number of k-mers per read
  m <- read.len - k + 1
  unit <- N / m * read.len
  seq.effort <- seq / unit
  result <- matrix(c(seq, f(seq.effort)), ncol=2, byrow=FALSE)
  colnames(result) <- c("bases", paste("frac(X>=", r, ")", sep=""))
  return(result)
}


## predict the fraction of k-mers represented at least r times in the sample
kmer.frac.bootstrap <- function(n, r=2, mt=20, times=30, conf=0.95) {
  return(preseqR.sample.cov.bootstrap(n=n, r=r-1, mt=mt, times=times, conf=conf))
}


## the fraction of k-mers represented at least r times as a function of 
## sample sizes
kmer.frac.curve.bootstrap <- function(n, k, read.len, seq, r=2, mt=20,
                                      times=30, conf=0.95)
{
  f <- kmer.frac.bootstrap(n, r=r, mt=mt, times=times, conf=conf)
  if (is.null(f))
    return(NULL)
  n[, 2] <- as.numeric(n[, 2])
  N <- n[, 1] %*% n[, 2]
  ## average number of k-mers per read
  m <- read.len - k + 1
  unit <- N / m * read.len
  seq.effort <- seq / unit
  result <- matrix(c(seq, f$f(seq.effort), f$lb(seq.effort), 
                     f$ub(seq.effort)), ncol=4, byrow=FALSE)
  colnames(result) <- c("bases", paste("frac(X>=", r, ")", sep=""), 
                        "lb", "ub")
  return(result)
}
