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

## the method refers to 
## Fisher, R. A., A. Steven Corbet, and C. B. Williams. "The Relation Between
## the Number of Species and the Number of Individuals in a Random Sample of an
## Animal Population." Journal of Animal Ecology 12, no. 1 (1943): 42-58.
## doi:10.2307/1411.

fisher.alpha <- function(n) {
  N <- n[, 1] %*% n[, 2]
  S <- sum(n[, 2])
  result <- uniroot(function(x) (exp(x) - 1) / x - N / S, c(0.001, 1e9), tol=0.0001, extendInt="upX")
  alpha <- S / result$root
  return(alpha)
}

fisher.mincount <- function(n, r=1) {
  alpha <- fisher.alpha(n)
  N <- n[, 1] %*% n[, 2]
  f.mincount <- function(t) {sapply(r, function(x) { alpha *
    integrate(function(z) (z^(x-1) / (1 - z)), lower=0, upper=N*t / (N*t + alpha))$value})}
  f.mincount(1)
  return(f.mincount)
}
