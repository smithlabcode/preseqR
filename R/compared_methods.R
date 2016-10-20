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

## zero truncated Poisson
## method ref Cohen, A. Clifford. (1960): 203-211.

ztpois.mincount <- function(n, r=1) {
    total.sample <- floor(n[, 1] %*% n[, 2])
    distinct <- sum(n[, 2])
    
    C <- n[, 1] %*% n[, 2] / sum(n[, 2])
    f <- function(x) {x / (1 - exp(-x))}
    result = uniroot(function(x) f(x) - C, c(0.001, 1e9), tol = 0.0001, extendInt="upX")
    lambda = result$root
    L <- sum(n[, 2]) / (1 - ppois(0, lambda))
    f.mincount <- function(t) {
      L * ppois(q=r - 1, lambda=lambda * t, lower.tail=FALSE)
    }
    f.mincount(1); f.mincount
}

## Boneh (1998)

boneh.mincount <- function(n, r=1) {
  total.sample <- floor(n[, 1] %*% n[, 2])
  distinct <- sum(n[, 2])

  tmp <- function(t) { sapply(r-1, function(x) {n[, 2] %*% (exp(-n[, 1]) - ppois(x, n[, 1] * t)) + distinct}) }
  
  index.f1 <- which(n[, 1] == 1)
  f1 <- n[index.f1, 2]
  U0 <- n[, 2] %*% exp(-(n[, 1]))
  if (length(index.f1) == 1 && f1 > U0) {
    result <- uniroot(function(x) x*(1 - exp(-f1 / x)) - U0, c(0.001, 1e9), tol=0.0001, extendInt="upX")
    U <- result$root
    f.mincount <- function(t) {tmp(t) + sapply(r-1, function(x) {U * (exp(-(f1 / U)) - ppois(x, f1 * t / U))})}
  } else {
    f.mincount <- tmp
  }
  f.mincount(1); f.mincount
}

## Chao and Shen (2004)

chao.mincount <- function(n, r=1, k=10) {
  total.sample <- floor(n[, 1] %*% n[, 2])
  distinct <- sum(n[, 2])
  index.f1 <- which(n[, 1] == 1)
  ## something wrong with the histogram
  if (length(index.f1) != 1)
    return(NULL)
  f1 <- n[index.f1, 2]
  index.rare <- which(n[, 1] <= k)
  S.rare <- sum(n[index.rare, 2])
  C.rare <- 1 - f1 / (n[index.rare, 1] %*% n[index.rare, 2])
  gamma.rare <- max(S.rare / C.rare * 
              ((n[index.rare, 1] * (n[index.rare, 1] - 1)) %*% n[index.rare, 2]) /
              (n[index.rare, 1] %*% n[index.rare, 2])^2 - 1, 0)
  f0 = S.rare / C.rare + f1 / C.rare * gamma.rare - S.rare
  f.mincount <- function(t) { sapply(r-1, function(x) { f0 + distinct  - f0 * ppois(x, f1 * (t-1) / f0) })}
  f.mincount(1); f.mincount
}
