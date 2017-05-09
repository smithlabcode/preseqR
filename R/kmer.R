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

kmer.frac <- function(n, r=2, mt=100)
{
  checking.hist(n)

  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation

  n[, 2] <- as.numeric(n[, 2])
  N <- n[, 1] %*% n[, 2]
  
  ## constructing the power series
  PS.coeffs <- generating.ps(n, mt=mt, j=1)

  if (is.null(PS.coeffs)) {
    write("the size of the initial experiment is insufficient", stderr())
    return(NULL)
  }

  ## use only the first mt terms in the power series
  mt <- min(mt, length(PS.coeffs))
  PS.coeffs <- PS.coeffs[ 1:mt ]

  ## check whether sample size is sufficient
  if (mt < 2)
  {
    m <- paste("max count before zero is les than min required count (2)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }
  
  ## construct the continued fraction approximation to the power seies
  cf <- ps2cfa(coef=PS.coeffs, mt=mt)
  rf <- cfa2rf(CF=cf)
  ## the length of cf could be less than mt
  ## even if ps do not have zero terms, coefficients of cf may have
  mt <- length(cf)
  ## select rational function approximants [M-1/M] m=2M
  ## asymptotically ~ C / t
  mt <- mt - (mt %% 2)
  valid.estimator <- FALSE
  m <- mt
  while (valid.estimator == FALSE) {

    rfa <- rf2rfa(RF=rf, m=m)
    ## solving roots
    numer.roots <- solve(rfa[[1]])
    denom.roots <- solve(rfa[[2]])
    ## seperating roots by their real parts
    numer.roots.neg <- numer.roots[which(Re(numer.roots) < 0)]
    numer.roots.pos <- numer.roots[which(Re(numer.roots) >= 0)]
    denom.roots.neg <- denom.roots[which(Re(denom.roots) < 0)]
    denom.roots.pos <- denom.roots[which(Re(denom.roots) >= 0)]

    ## record roots in the numerator that are significantly similar to
    ## roots in the denominator
    tmp.roots <- c()

    ## simplify the rational function approximation
    ## two roots are same if the difference is less than the 
    ## predefined PRECISION
    if (length(numer.roots.pos) > 0) {
      for (i in 1:length(numer.roots.pos)) {
        if (length(denom.roots.pos) > 0) {
          d <- Mod(denom.roots.pos - numer.roots.pos[i])
          for (j in 1:length(d)) {
            if (d[j] < PRECISION) {
              denom.roots.pos <- denom.roots.pos[-j]
              tmp.roots <- c(tmp.roots, numer.roots.pos[i])
              break
            }
          }
        }
      }
    }

    ## roots in simplified RFA
    numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
    denom.roots <- c(denom.roots.neg, denom.roots.pos)

    ## convert roots from t - 1 to t
    roots <- denom.roots + 1
    ## pacman rule checking
    if (length(which(roots == 0)) || length(which(Re(roots) > 0))) {
      m <- m - 2
      next
    } else {
      poly.numer <- as.function(poly.from.roots(numer.roots))
      l <- length(denom.roots)
      ## treat polynomials in the rational function to be monic
      ## the difference to the original RFA is a multiplier C

      ## c_i in the estimator
      coef <- sapply(1:l, function(x) {
        poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
      ## check whether the estimator is non-decreased
      ## NOTE: it only checks for t >= 1 !!!
      deriv.f <- function(t) {
        Re(sapply(t, function(x) {-(coef*roots) %*% ( 1 / ((x-denom.roots)^2))}))} 
      if (length(which( deriv.f(seq(0.05, 100, by=0.05)) < 0 ) != 0)) {
        m <- m - 2
        next
      }
      ## calculate the constant C
      C <- coef(rfa[[1]])[length(coef(rfa[[1]]))] / 
           coef(rfa[[2]])[length(coef(rfa[[2]]))]
      ## species accum curves with minimum count r
      ## using parital fraction expansion
      denom.roots <- denom.roots + 1
      coef <- coef * C
      frac.bias <- -Re(sum(coef / denom.roots / N))

      f.frac <- function(t) {
        sapply(r, function(x) {
          Re((x * coef) %*% ((t / (t - denom.roots))^(x-1) / (t - denom.roots)) / N - 
             (coef / denom.roots) %*% (t / (t - denom.roots))^x / N)})}
      ## adjust the bias
      ## make sure the value of the function is 1 when r = 1
      f.frac.adjust <- function(t) {
        f.frac(t) / frac.bias
      }
      valid.estimator <- TRUE
    }
  }
  f.frac.adjust(1)
  return(f.frac.adjust)
}

## N, number of nucleotides sequenced in the initial experiment
## seq.size.GB, the number of nucleotides plan to sequence
## the parameter can take an array of numbers
kmer.frac.curve <- function(n, N, seq.size.GB, r=2, mt=100)
{
  f <- kmer.frac(n, r=r, mt=mt)
  if (is.null(f)) return(NULL)
  seq.effort <- seq.size.GB * 1e9 / N
  frac <- t(sapply(seq.effort, function(x) f(x)))
  curves <- cbind(seq.size.GB, frac)
  colnames(curves) <- c("bases(GB)", paste("frac(X>=", r, ")", sep=""))
  return(curves)
}

kmer.frac.bootstrap <- function(n, r=2, mt=100, times=100)
{
  n[, 2] <- as.numeric(n[, 2])
  ## total individuals
  total <- n[, 1] %*% n[, 2]

  ## returned function
  frac <- vector(length=times, mode="list")

  ds.estimator <- function(n, r, mt, t.scale) {
    f <- kmer.frac(n, r=r, mt=mt)
    function(t) {f(t * t.scale)}
  }

  while (times > 0) {
    n.bootstrap <- matrix(c(n[, 1], rmultinom(1, sum(n[, 2]), n[, 2])), ncol=2)
    total.bootstrap <- n.bootstrap[, 1] %*% n.bootstrap[, 2]
    t.scale <- total / total.bootstrap
    f <-  ds.estimator(n.bootstrap, r=r, mt=mt, t.scale=t.scale) 

    frac[[times]] <- f
    ## prevent later binding!!!
    frac[[times]](1)
    times <- times - 1
  }
  f.estimator <- kmer.frac(n=n, r=r, mt=mt)
  if (length(r) == 1) {
    median.estimators <- function(t) {median( sapply(frac, function(x) x(t)) )}
    var.estimator <- function(t) {var( sapply(frac, function(x) x(t)) )}
  } else {
    median.estimators <- function(t) {apply(sapply(frac, function(x) x(t)), FUN=median, MARGIN=1)}
    var.estimator <- function(t) {apply(sapply(frac, function(x) x(t)), FUN=var, MARGIN=1)}
  }
  ## prevent later binding!!!
  f.estimator(1); median.estimators(1); var.estimator(1)
  return(list(FUN.nobootstrap=f.estimator, FUN.bootstrap=median.estimators, var=var.estimator))
}
