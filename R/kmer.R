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

preseqR.kmer.frac <- function(n, r=2, mt=100)
{
  ## setting the diagonal value
  di <- 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  n[, 2] <- as.numeric(n[, 2])
  N <- n[, 1] %*% n[, 2]
  
  ## constructing the power series
  PS.coeffs <- generating.ps(n, 1, mt=mt)

  if (is.null(PS.coeffs)) {
    write("the size of the initial experiment is insufficient", stderr())
    return(NULL)
  }

  ## constrain the continued fraction approximation with even terms
  ## asymptotically ~ C / t
  mt <- min(mt, length(PS.coeffs))
  mt <- mt - (mt %% 2)
  PS.coeffs <- PS.coeffs[ 1:mt ]

  ## check whether sample size is sufficient
  if (mt < 2)
  {
    m <- paste("max count before zero is les than min required count (4)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }

  ## construct a continued fraction approximation including as many as possible
  ## terms
  valid <- FALSE
  if (mt >= MIN_REQUIRED_TERMS) {
    DE <- seq(mt, MIN_REQUIRED_TERMS, by=-2)
    for (de in DE) {
      ## continued fraction approximation to a power series
      out <- .C('c_PS2CF', as.integer(di), 
              as.integer(de), as.double(PS.coeffs[1:de]), 
              as.integer(length(PS.coeffs[1:de])),
              ps.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
              ps.coeffs.l = as.integer(0),
              cf.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
              cf.coeffs.l = as.integer(0),
              offset.coeffs =as.double(vector(mode='numeric',length=MAXLENGTH)),
              diagonal.idx = as.integer(0),
              degree = as.integer(0), is.valid = as.integer(0));
      if (out$is.valid) {break}
    }
    if (out$is.valid) {
      length(out$ps.coeffs) <- out$ps.coeffs.l
      length(out$cf.coeffs) <- out$cf.coeffs.l
      length(out$offset.coeffs) <- as.integer(abs(out$diagonal.idx))
      CF.space <- list(out$ps.coeffs, out$cf.coeffs, out$offset.coeffs, 
                     out$diagonal.idx, out$degree)
      names(CF.space) <- c('ps.coeffs', 'cf.coeffs', 'offset.coeffs', 'diagonal.idx',
                         'degree')
      DE = seq(CF.space$degree, MIN_REQUIRED_TERMS, by=-2)

      for (de in DE) {
        CF <- list(CF.space$ps.coeffs[1:de], CF.space$cf.coeffs[1:de],
                 CF.space$offset.coeffs, CF.space$diagonal.idx, de)
        names(CF) <- c('ps.coeffs', 'cf.coeffs', 'offset.coeffs', 'diagonal.idx',
                     'degree')
        class(CF) <- 'CF'
        ## convert the continued fraction to the RFA 
        RF <- CF2RFA(CF)
        RF[[1]] <- RF[[1]] / polynomial(c(0, 1))

        ## solving roots
        numer.roots <- solve(RF[[1]])
        denom.roots <- solve(RF[[2]])
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
            next
          }
          ## the estimator passes the requirement
          valid <- TRUE
          ## calculate the constant C
          C <- coef(RF[[1]])[length(coef(RF[[1]]))] / 
               coef(RF[[2]])[length(coef(RF[[2]]))]
          ## species accum curves with minimum count r
          ## using parital fraction expansion
          denom.roots <- denom.roots + 1
          coef <- coef * C

          frac.bias <- -Re(sum(coef / denom.roots / N))

          f.frac <- function(t) {
            sapply(r, function(x) {
              Re((x * coef) %*% ((t / (t - denom.roots))^(x-1) / (t - denom.roots)) / N - 
                 (coef / denom.roots) %*% (t / (t - denom.roots))^x / N)})}
          break
        }
      }
    }
  }
  if (valid==FALSE) {
    s1 <- sum(n[, 2])
    s2 <- sum(n[n[, 1] > 1, 2])
    coef <- s1^2 / s2
    denom.roots <- (s1 - s2) / s2
    frac.bias <- -Re(sum(coef / denom.roots / N))
    f.frac <- function(t) {
      sapply(r, function(x) {
        Re((x * coef) %*% ((t / (t - denom.roots))^(x-1) / (t - denom.roots)) / N - 
          (coef / denom.roots) %*% (t / (t - denom.roots))^x / N)})}
  }

  ## add a correction
  d <- 1 - frac.bias
  if (d >= 0) {
    f.frac.adjust <- function(t) {
      f.frac(t) + ifelse(t <= 1, d * t, d)
    }
  } else {
    return(NULL)
  }
}


preseqR.kmer.frac.curve <- function(n, N, seq.size.GB, r=2, mt=100)
{
  f <- preseqR.kmer.frac(n, r=r, mt=mt)
  if (is.null(f)) return(NULL)
  seq.effort <- seq.size.GB * 1e9 / N
  kmer.frac <- t(sapply(seq.effort, function(x) f(x)))
  curves <- cbind(seq.size.GB, kmer.frac)
  colnames(curves) <- c("bases(GB)", paste("frac(X>=", r, ")", sep=""))
  return(curves)
}


preseqR.kmer.frac.bootstrap <- function(n, r=1, mt=100, times=100)
{
  BOOTSTRAP.factor <- 0.4
  n[, 2] <- as.numeric(n[, 2])
  ## total individuals
  total <- n[, 1] %*% n[, 2]

  ## the number of resampling times                                                 
  counter <- 0
  ## returned function
  f.mincount <- vector(length=times, mode="list")

  ## upperbound of times of iterations for bootstrapping
  upper.limit <- times / BOOTSTRAP.factor

  ds.estimator <- function(n, r, mt, t.scale) {
    f <- preseqR.kmer.frac(n, r=r, mt=mt)
    if (is.null(f)) {
      return(NULL)
    } else {
      function(t) {f(t * t.scale)}
    }
  }

  while (times > 0) {
    n.bootstrap <- matrix(c(n[, 1], rmultinom(1, sum(n[, 2]), n[, 2])), ncol=2)
    total.bootstrap <- n.bootstrap[, 1] %*% n.bootstrap[, 2]
    t.scale <- total / total.bootstrap
    f <-  ds.estimator(n.bootstrap, r=r, mt=mt, t.scale=t.scale) 
    counter <- counter + 1
    if (!is.null(f)) {
      f.mincount[[times]] <- f
      ## avoid late binding
      f.mincount[[times]](1)
      times <- times - 1
    }
    if (counter > upper.limit)
      break
  }
  if (times > 0) {
     write("fail to bootstrap!", stderr())
    return(NULL)
  } else {
    if (length(r) == 1)
      return(function(t) {median( sapply(f.mincount, function(x) x(t)) )})
    return( function(t) {apply(sapply(f.mincount, function(x) x(t)), FUN=median, MARGIN=1)} )
  }
}
