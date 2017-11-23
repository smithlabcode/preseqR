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
  n, efficiency=0.05, bin=1e8, r=1, mt=20, size=SIZE.INIT,
  mu=MU.INIT, times=30, conf=0.95)
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
  f.rSAC <- preseqR.rSAC.bootstrap(
    n=n, r=r, mt=mt, size=size, mu=mu,times=times, conf=conf)

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
