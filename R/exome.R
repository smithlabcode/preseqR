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


pois.mincount <- function(n, L, r=1) {
    lambda <- n[, 1] %*% n[, 2] / L
    f.mincount <- function(t) {
      L * ppois(q=r - 1, lambda=lambda * t, lower.tail=FALSE)
    }
    f.mincount(1); f.mincount
}


### calculate the negative binomial loglikelihood
### zero.items is number of items unobserved
### size and mu are parameters in a negative binomial distribution
nb.loglikelihood <- function(n, zero.items, size, mu)
{
  ## likelihood of nonzero terms
  log.prob <- dnbinom(n[, 1], size = size, mu = mu, log = TRUE)
  loglikelihood <- log.prob %*% n[, 2]

  ## add items with zero count
  log.zero.prob <- dnbinom(0, size = size, mu = mu, log = TRUE)
  loglikelihood <- loglikelihood + zero.items * log.zero.prob

  return(loglikelihood)
}

nb.fitting <- function(n, L)
{
  n[, 2] <- as.numeric(n[, 2])

  ## the number of unobservations
  zero.items <- L - sum(n[, 2])

  ## estimated mean and variance
  m <- (n[, 1] %*% n[, 2]) / L
  v <- ( (n[, 1] - m)^2 %*% n[, 2] + m^2 * zero.items )/(L - 1)

  ## target function f
  f <- function(x) {
        return( -nb.loglikelihood(n, zero.items, size = x, mu = m)/L )
  }

  ## derivative of f
  gr <- function(x)
  {
    first.term <- ( digamma(x) * zero.items +
                    digamma(n[, 1] + x) %*% n[, 2] )/L
    second.term <- digamma(x)
    third.term <- log(x) - log(x + m)
    result <- first.term - second.term + third.term
    # f is negative loglikelihood
    return(-result)
  }

  ## estimate size and mu based on first and second moments
  if (v > m) {
    res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B",
           lower = 0.00001, upper = 100000)
  } else {
    res <- optim(size, f, gr, method = "L-BFGS-B",
           lower = 0.00001, upper = 100000)
  }

  loglikelihood <- nb.loglikelihood(n, zero.items, size=res$par, mu=m)
  ## update parameters
  size <- res$par
  mu <- m

  return(list(size = size, mu = mu, loglik = -loglikelihood))
}


## fitting the negative binoimal distribution to the data
## ss is the step.size
## max.extrapoltion is the maximum value for extrapolation
## r is a vector of frequencies
## L is the total number of species
nb.mincount <- function(n, L, r=1)
{
  n[, 2] <- as.numeric(n[, 2])

  ## estimate parameters
  opt <- nb.fitting(n, L)
  size <- opt$size
  mu <- opt$mu

  f.mincount <- function(t) {
    L * pnbinom(r-1, size=size, mu=mu*t, lower.tail=FALSE)
  }
  f.mincount(1); f.mincount
}


ds.mincount <- function(n, r=1, mt=100)
{
  ## setting the diagonal value
  di <- 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  n[, 2] <- as.numeric(n[, 2])

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
  if (mt < MIN_REQUIRED_TERMS)
  {
    m <- paste("max count before zero is les than min required count (4)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }

  ## construct a continued fraction approximation including as many as possible
  ## terms
  valid <- FALSE
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
        f.mincount <- function(t) {
          sapply(r, function(x) {
              Re(coef %*% (t / (t - denom.roots))^x)})}
        break
      }
    }
  }
  if (valid==FALSE)
    return(NULL)
  f.mincount(1)
  list(FUN=f.mincount, m=de)
}


## nonparametric approach Deng & Smith 2016
ds.mincount.bootstrap <- function(n, r=1, mt=100, times=100)
{
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
    f <- ds.mincount(n, r=r, mt=mt)
    if (is.null(f)) {
      return(NULL)
    } else {
      function(t) {f$FUN(t * t.scale)}
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
      ## prevent later binding!!!
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
    f.estimator <- ds.mincount(n=n, r=r, mt=mt)
    if (length(r) == 1) {
      median.estimators <- function(t) {median( sapply(f.mincount, function(x) x(t)) )}
      var.estimator <- function(t) {var( sapply(f.mincount, function(x) x(t)) )}
    } else {
      median.estimators <- function(t) {apply(sapply(f.mincount, function(x) x(t)), FUN=median, MARGIN=1)}
      var.estimator <- function(t) {apply(sapply(f.mincount, function(x) x(t)), FUN=var, MARGIN=1)}
    }
    if (!is.null(f.estimator)) f.estimator$FUN(1); median.estimators(1); var.estimator(1)
    return(list(f=f.estimator, median=median.estimators, var=var.estimator))
  }
}

# write out the information about the experiment and the number of reads needs
# to be sequenced
preseqR.depthseq <- function(n.reads, FUN=NULL, LIB.FUN=NULL, L=NULL, rho=0.85, r=8, uniq=TRUE)
{
  checking.hist(n.reads)

  if (is.null(FUN)) {
    cat(paste("No estimator for the number of nucleotides with at least ", r, 
        " aligned reads as a function of sequencing effort\n", sep=""))
    return(NULL)
  }
  if (is.null(LIB.FUN)) {
    cat("No estimator for the library complexity curve\n") 
    return(NULL)
  }
  if (is.null(L)) {
    cat("The length of the targeted regions is not specified\n")
    return(NULL)
  }
  N <- n.reads[, 1] %*% n.reads[, 2]
  cat(paste("The number of reads is ", N, " in the initial experiment\n", sep=""))
  if (FUN(10000) <= L * rho) {
    cat(paste("More than 10000 times the size of the initial experiment is",
              " needed to achieve the required standard\n", sep=""))
    return(NULL)
  }
  t0 <- uniroot(function(x) {FUN(x) - L * rho}, interval=c(0, 10000), tol=0.00001)$root
  ## whether remove the duplicates when counting coverage depth
  ## Counting with duplicates
  if (uniq==FALSE) {
    reads.total <- ceiling(N * t0 / 1000000.0)
    cat(paste("In order to attain ", rho*100, "% of the targeted regions of length ", 
        L, " with ", r, "X or greater coverage depth\n", sep=""))
    cat(paste("A total of ", reads.total, " million reads are needed in the full experiment\n", sep=""))
    return(reads.total)
  }
  ## duplicates are removed
  N.uniq <- sum(n.reads[, 2])
  if (LIB.FUN(10000) <= N.uniq * t0) {
    cat(paste("More than 10000 times the size of the initial experiment is",
              " needed to achieve the required standard\n", sep=""))
    return(NULL)
  }
  t1 <- uniroot(function(x) {LIB.FUN(x) - N.uniq * t0}, interval=c(0, 10000), tol=0.00001)$root
  reads.total <- ceiling(N * t1 / 1000000.0)
  cat(paste("In order to attain ", rho*100, "% of the targeted regions of length ", 
      L, " with ", r, "X or greater coverage depth\n", sep=""))
  cat(paste("A total of ", reads.total, " million reads are needed in the full experiment\n", sep=""))
  return(reads.total)
}
