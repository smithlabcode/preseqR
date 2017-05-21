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

## continued fraction approximant to a power series
## QD algorithm
## input: coefficients of the power series; begin with the constant
## mt: the maximum number of terms used in the power series
ps2cfa <- function(coef, mt) {
  index <- which(coef == 0)
  if (length(index) == 0) {
    mt <- min(mt, length(coef))
  } else {
    mt <- min(mt, index[1] - 1)
  }
  if (mt == 1) {
    return(coef[1])
  }
  qd.table <- matrix(data=0, nrow=mt, ncol=mt)
  ## initialize the table
  ## the first column is 0
  qd.table[1:(mt-1), 2] <- coef[2:mt] / coef[1:(mt-1)]
  if (mt == 2) {
    return(c(coef[1], -qd.table[1, 2]))
  }
  ## two types of columns e or q
  for (i in 3:mt) {
    n <- mt + 1 - i
    if (i %% 2 == 1) {
      ## number of entries in the column
      qd.table[1:n, i] <- qd.table[2:(n+1), i-1] - qd.table[1:n, i-1] + qd.table[2:(n+1), i-2]
      if (!is.finite(qd.table[1, i]) || qd.table[1, i] == 0) 
        return(c(coef[1], -qd.table[1, 2:(i-1)]))
    } else {
      qd.table[1:n, i] <- qd.table[2:(n+1), i-1] / qd.table[1:n, i-1] * qd.table[2:(n+1), i-2]
      if (!is.finite(qd.table[1, i]) || qd.table[1, i] == 0)
        return(c(coef[1], -qd.table[1, 2:(i-1)]))
    }
  }
  return(c(coef[1], -qd.table[1, 2:mt]))
}


## convert truncated continued fraction to a series of rational functions
## output two sets: A for numerators and B for denumerators
cfa2rf <- function(CF) {
  ## A, B are sets of polynomials based on recursive formula
  A <- list()
  B <- list()
  if (length(CF) < 2) {
    return(polynomial(CF))
  }
  A[[1]] <- polynomial(CF[[1]])
  A[[2]] <- polynomial(CF[[1]])
  B[[1]] <- polynomial(1)
  B[[2]] <- polynomial(c(1, CF[[2]]))
  if (length(CF) == 2) {
    return(list(A=A, B=B))
  }
  for (i in 3:length(CF)) {
    A[[i]] <- A[[i-1]] + polynomial(c(0, CF[[i]])) * A[[i-2]]
    B[[i]] <- B[[i-1]] + polynomial(c(0, CF[[i]])) * B[[i-2]]
  }
  return(list(A=A, B=B))
}

## Pad\'{e} approximant by picking out the numerator and the denominator
## input: two sets of polynomials for numerators and denominators
##        the degree m
## output: rational function approximant or Pad\'{e} approximant
rf2rfa <- function(RF, m) {
  return(polylist(RF$A[[m]], RF$B[[m]]))
}
