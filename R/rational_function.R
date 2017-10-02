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

## continued fraction approximant to a power series based on
## QD algorithm
## coefs, coefficients of the power series; 
## mt, the number of terms in the power series used for constructing
## the continued fraction approximation
## ref pp. 131, 147 and 148 in the book Pad\'{e} Approximants 2ed
ps2cfa <- function(coefs, mt) {
  ## use nonzero terms which are required by QD algorithm
  index <- which(coefs == 0)
  if (length(index) == 0) {
    mt <- min(mt, length(coefs))
  } else {
    mt <- min(mt, index[1] - 1)
  }
  if (mt == 1) {
    return(coefs[1])
  }
  ## QD algorithm
  qd.table <- matrix(data=0, nrow=mt, ncol=mt)
  ## initialize the table
  ## the first column is 0
  qd.table[1:(mt-1), 2] <- coefs[2:mt] / coefs[1:(mt-1)]
  if (mt == 2) {
    return(c(coefs[1], -qd.table[1, 2]))
  }
  ## two types of columns e or q
  for (i in 3:mt) {
    n <- mt + 1 - i
    if (i %% 2 == 1) {
      ## number of entries in the column
      qd.table[1:n, i] <- qd.table[2:(n+1), i-1] - qd.table[1:n, i-1] + 
                          qd.table[2:(n+1), i-2]
      if (!is.finite(qd.table[1, i]) || qd.table[1, i] == 0) 
        return(c(coefs[1], -qd.table[1, 2:(i-1)]))
    } else {
      qd.table[1:n, i] <- qd.table[2:(n+1), i-1] / qd.table[1:n, i-1] * 
                          qd.table[2:(n+1), i-2]
      if (!is.finite(qd.table[1, i]) || qd.table[1, i] == 0)
        return(c(coefs[1], -qd.table[1, 2:(i-1)]))
    }
  }
  return(c(coefs[1], -qd.table[1, 2:mt]))
}


## convert truncated continued fraction to a series of rational functions
## numerators are stored in set A and denumerators are stored in set B
## equation (2.14a), (2.14b), (2.15) in the book Pad\'{e} Approximants 2ed
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
## output: Pad\'{e} approximant
rf2rfa <- function(RF, m) {
  return(polylist(RF$A[[m]], RF$B[[m]]))
}


## simplify the rational function, eliminate defects and partial-fraction
## decompoistion
rfa.simplify <- function(rfa) {
  ## solving roots
  numer.roots <- solve(rfa[[1]])
  denom.roots <- solve(rfa[[2]])

  ## finite
  if (any(!is.finite(c(numer.roots, denom.roots))))
    return(NULL)

  ## identify defects
  ## the root and the pole is a defect if the difference is less than 
  ## the predefined precision, which is defined by the variable PRECISION
  tmp.roots <- c()
  for (i in 1:length(denom.roots)) {
    if (length(numer.roots) > 0) {
      d <- Mod(denom.roots[i] - numer.roots)
      ind <- which.min(d)
      if (d[ind] < PRECISION) {
        numer.roots <- numer.roots[-ind]
        tmp.roots <- c(tmp.roots, denom.roots[i])
      }
    }
  }

  ## eliminate defects
  denom.roots <- denom.roots[!denom.roots %in% tmp.roots]
  ## convert roots from t - 1 to t
  poles <- denom.roots + 1
    
  ## treat both numerator and denuminator in the rational function as
  ## monic polynomials
  ## the difference from the original rational function is up to a factor
  if (length(numer.roots) == 0) {
    poly.numer <- as.function(polynomial(1))
  } else {
    poly.numer <- as.function(poly.from.roots(numer.roots))
  }
  l <- length(denom.roots)

  ## coefficients in the partial fraction
  coefs <- sapply(1:l, function(x) {
    poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
  ## calculate the factor
  C <- coef(rfa[[1]])[length(coef(rfa[[1]]))] / 
       coef(rfa[[2]])[length(coef(rfa[[2]]))]
  coefs <- coefs * C
  return(list(coefs=coefs, poles=poles))
}


## modified Pad\'{e} approximant
## close to the average discovery rate and satisfies
## the sum of estimates of E(S_r(t)) for r >= 1 is equal to Nt
## require mt to be odd
rfa.sample.cov <- function(n, mt) {
  mt <- mt - (mt + 1) %% 2
  PS.coeffs <- discoveryrate.ps(n, mt)
  m <- (mt + 1) / 2
  N <- n[, 1] %*% n[, 2]
  if (mt == 1) {
    a = sum(n[, 2])
    b = c(1, (N - a) / N)
  } else {
    ## system equation ax = b
    A <- t(sapply(1:(m - 1), function(x) PS.coeffs[x:(x + m - 1)]))
    ## the last equation is adjust to make sure the sum
    last.eqn.coefs <- N
    for (i in 1:(m-1))
      last.eqn.coefs <- c(last.eqn.coefs, PS.coeffs[i] - last.eqn.coefs[length(last.eqn.coefs)])
    A <- rbind(A, last.eqn.coefs)
    b0 <- -c(PS.coeffs[(m+1):mt], PS.coeffs[m] - last.eqn.coefs[length(last.eqn.coefs)])
    b <- solve(A, b0)
    b <- c(1, rev(b))
    a <- sapply(1:m, function(x) b[1:x] %*% rev(PS.coeffs[1:x]))
  }
  return(polylist(polynomial(a), polynomial(b)))
}


## discriminant of the quadratic polynomial, which is
## the denominator of the discovery rate at m = 2
## OBSOLATE 
discriminant <- function(n) {
  if (max(n[, 1]) < 3) {
    return(NULL)
  }
  n[, 2] <- as.numeric(n[, 2])
  S1 <- sum(n[, 2])
  if (length(which(n[, 1] == 1))) {
    S2 <- S1 - n[which(n[, 1] == 1), 2]
  } else {
  	S2 <- S1
  }
  if (length(which(n[, 1] == 2))) {
    S3 <- S2 - n[which(n[, 1] == 2), 2]
  } else {
  	S3 <- S2
  }
  if (length(which(n[, 1] == 3))) {
    S4 <- S3 - n[which(n[, 1] == 3), 2]
  } else {
  	S4 <- S3
  }
  a <- S2*S4 - S3^2
  b <- S1*S4 - S2*S3
  c <- S1*S3 - S2^2
  return((b / a)^2 - 4 * (c / a))
}
