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


## predict the optimal number of sequenced bases using cost-benefit ratio
preseqR.optimal.sequencing <- function(
  n, efficiency=0.05, bin=1e8, r=1, mt=20, times=30, conf=0.95)
{
  find.start <- function(f, N, bin, efficiency) {
    y = sapply(1:100, function(x) (f(x + bin / N) - f(x)) / bin - efficiency)
    if (length(which(y > 0)) == 0) {
      return(-1)
    }
    index = min(which(y > 0))
    return(index)
  }

  n[, 2] <- as.numeric(n[, 2])
  N <- n[, 1] %*% n[, 2]

  ## r-species accumulation curve as a function of relative sample size
  f.rSAC <- ds.rSAC.bootstrap(
    n=n, r=r, mt=mt, times=times, conf=conf)

  ## hint: using r-SAC as a function of the number of sequenced bases
  f <- f.rSAC$f
  lb <- f.rSAC$lb
  ub <- f.rSAC$ub
  start = find.start(f=f, N=N, bin=bin, efficiency=efficiency)
  if (start == -1) {
    t0 = 0
  } else {
    t0 <- uniroot(function(x) {(f(x+bin/N) - f(x)) / bin - efficiency}, 
                  interval=c(start, 10000), tol=0.01)$root 
  }
  start = find.start(f=lb, N=N, bin=bin, efficiency=efficiency)
  if (start == -1) {
    t1 = 0
  } else {
    t1 <- uniroot(function(x) {(lb(x+bin/N) - lb(x)) / bin - efficiency},
                     interval=c(start, 10000), tol=0.01)$root
  }
  start = find.start(f=ub, N=N, bin=bin, efficiency=efficiency)
  if (start == -1) {
    t2 = 0
  } else {
    t2 <- uniroot(function(x) {(ub(x+bin/N) - ub(x)) / bin - efficiency},
                     interval=c(start, 10000), tol=0.01)$root
  }
  t0_lb = min(t1, t2)
  t0_ub = max(t1, t2)
  return(c(t0, t0_lb, t0_ub) * N)
}


## the function is designed for EXOME sequencing, where aligned reads that
## map to the same location are removed to avoid potential duplicate
preseqR.rSAC.sequencing.rmdup <- function(
  n_base, n_read, r=1, mt=20, times=30, conf=0.95)
{
  checking.hist(n_read)
  checking.hist(n_base)
  n_read[, 2] <- as.numeric(n_read[, 2])
  N <- n_read[, 1] %*% n_read[, 2]

  ## the bootstrappedestiamtor
  f.bootstrap <- function(n, r, mt) {
    n.bootstrap <- matrix(c(n[, 1], rmultinom(1, sum(n[, 2]), n[, 2])), ncol=2)
    N.bootstrap <- n.bootstrap[, 1] %*% n.bootstrap[, 2]
    N <- n[, 1] %*% n[, 2]
    t.scale <- N / N.bootstrap
    f <- ds.rSAC(n.bootstrap, r=r, mt=mt)
    return(function(t) {f(t * t.scale)})
  }

  ## returned function
  f.rSACs <- vector(length=times, mode="list")
  ## the number of bases with coverage at least r as a function of sequenced
  ## bases from uniquely aligned reads
  f.bases <- vector(length=times, mode="list")
  ## the number of unique reads as a function of the number of reads
  f.reads <- vector(length=times, mode="list")

  ## the number of nucleotides with coverage at least r based on uniquely 
  ## aligned reads as a function of sequencing effort
  f <- function(n_read, f.base, f.read) {
    N <- sum(as.double(n_read[, 2]))
    return(function(t) {t0 = f.read(t) / N; f.base(t0)})
  }

  while (times > 0) {
    f.bases[[times]] <- f.bootstrap(n=n_base, r=r, mt=mt)
    f.reads[[times]] <- f.bootstrap(n=n_read, r=1, mt=mt)
    f.rSACs[[times]] <- f(n_read=n_read, f.base=f.bases[[times]], f.read=f.reads[[times]])
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
