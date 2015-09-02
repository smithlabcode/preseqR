## Two roots are the same if the difference is less than the PRECISION
PRECISION = 1e-3

### calculate the coefficients of the power series with freq >= r given the
### vector of count frequencies
mincount.ps <- function(hist.count, r) {
  max.freq <- length(hist.count)
  res <- vector(length = max.freq, mode = 'numeric')
  for (i in 1:max.freq) {
    co.eff <- 0
    for ( l in 0:(r-1) ) {
      for (j in 0:i) {
        index <- j + l
        if (i <= index && index <= max.freq) {
          co.eff <- co.eff + (-1)^(j + 1) * hist.count[index] *
                   exp( lgamma(j + l + 1) - lgamma(j + 1) - lgamma(l + 1) +
                        lgamma(l + 1) - lgamma(i-j+1) - lgamma(l - i + j + 1) )
        }
      }
    }
    res[i] <- co.eff
  }
  res
}

### minimum count interpolation
### ss step size
### n two-column histogram
preseqR.interpolate.mincount <- function(ss, n, r=1)
{
  checking.hist(n)

  ## total individuals captured
  total.sample <- n[, 1] %*% n[, 2]
  N <- floor(total.sample)

  ## total species
  initial.distinct <- sum(n[, 2])
  step.size <- as.double(ss)

  ## l is the number of points for interpolation
  l <- as.integer(N / step.size)

  r <- floor(r)
  ## if the sample size is larger than the size of experiment or 
  ## the step size is too small, return NULL
  if (l == 0 || ss < 1 || r < 1)
    return()

  ## explicitly calculating the expectation species observed at least r times
  ## based on sampling without replacement
  ## see K.L Heck 1975
  ## N total individuals
  ## S the number of species
  ## size the size of the subsample
  expect.distinct <- function(n, N, size, S, r) {
    denom <- lchoose(N, size)
    p <- sapply(n[, 1], function(x) {
           sum(exp(lchoose(N - x, size - 0:(r-1)) + lchoose(x, 0:(r-1)) - denom))})
    return(S - p %*% n[, 2])
  }

  ## sample size vector
  x <- step.size * ( 1:l )

  ## calculate the number of distinct reads based on each sample size
  yield.estimates <- sapply(x, function(x) expect.distinct(n, N, x, 
                            initial.distinct, r))

  ## put size and yield together into a matrix
  result <- matrix(c(x, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'interpolation')

  return(result)
}


### convert a continued fraction to a rational function
### The form of the continued fraction refers to the Supplementary of 
### Daley, T., & Smith, A. D.(2013)
### the form of continued fractions e.g. a_0(t - 1) / (1 + a_1(t - 1))
CF2RFA <- function(CF)
{
  poly.numer <- polynomial(1)
  poly.denom <- polynomial(1)
  cf <- CF$cf.coeffs
  for (i in rev(cf)) {
    tmp <- poly.numer
    poly.numer <- poly.numer + polynomial(c(0, i)) * poly.denom
    poly.denom <- tmp
  }
  ## according to the representation of the continued fraction
  ## the first item is a_0x/... not 1 + a_0x/...
  ## thus we need to substract one from the result
  poly.numer <- poly.numer - poly.denom
  polylist(poly.numer, poly.denom)
}


### the kth derivative of a continued fraction
deriv.CF <- function(CF, k = 0)
{
  ## check CF is a continued fraction with CF attribute
  if (class(CF) != "CF")
    return(NULL)
  rational.f <- CF2RFA(CF)

  while (k > 0) {
    ## quotient rule
    rational.f[[1]] <- deriv(rational.f[[1]]) * rational.f[[2]] - 
                       rational.f[[1]] * deriv(rational.f[[2]])
    rational.f[[2]] <- rational.f[[2]]^2
    k <- k - 1
  }
  rational.f
}


### species accumu curve with minimum count r
### by RFA to a generalized Good-Toulmin power series
### RFA is based on pacman rule
### n two-column histogram
### ss step size
### mt maximum terms used
### r minimum count
### CHAO: old method
preseqR.rfa.mincount <- function(n, mt = 50, ss = NULL,
                                 max.extrapolation = NULL, r=1)
{
  checking.hist(n)
  ## setting the diagonal value
  di = 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  ## calculate total number of sample
  total.sample <- n[, 1] %*% n[, 2]
  total.sample <- floor(total.sample)

  ## set step.size as the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- floor(total.sample)
    step.size <- ss
  } else if (ss < 1) {
    write("The step size is too small", stderr())
    return()
  } else {
    step.size <- floor(ss)
  }

  ## no interpolation if step.size is larger than the size of experiment
  ## set the starting sample size as the step.size
  if (step.size > total.sample) {
    yield.estimates <- vector(mode = 'numeric', length = 0)

    ## starting sample size for extrapolation
    starting.size <- step.size
  } else {
      ## interpolation when sample size is no more than total sample size
      ## interpolate and set the size of the sample for an initial extrapolation
      out <- preseqR.interpolate.mincount(step.size, n, r)
      yield.estimates <- out[, 2]

      ## starting sample size for extrapolation
      starting.size <- ( floor(total.sample/step.size) + 1 )*step.size
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  ## transform a histogram into a vector of frequencies
  hist.count <- vector(length=max(n[, 1]), mode="numeric")
  hist.count[n[, 1]] <- n[, 2]

  ## coefficients of the generalized power series
  PS.coeffs <- mincount.ps(hist.count, r)

  ## only use power series with non-zero coefficients
  ## effective coefficients from begining up to the first zero
  counts.before.first.zero = 1
  while (as.integer(counts.before.first.zero) <= length(PS.coeffs) &&
         PS.coeffs[counts.before.first.zero] != 0)
    counts.before.first.zero <- counts.before.first.zero + 1

  ## for r > 1, the jth coefficient in the power series requires
  ## n_{i} i = j, j + 1, ..., j + r - 1 itmes 
  ## constrain the continued fraction approximation with even degree 
  ## conservatively estimates
  mt <- min(mt, counts.before.first.zero - r)
  mt <- mt - (mt %% 2)

  ## pre-check whether sample is sufficient
  if (mt < MIN_REQUIRED_TERMS)
  {
    m <- paste("max count before zero is les than min required count (4)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }

  PS.coeffs <- PS.coeffs[ 1:mt ]

  ## construct a continued fraction approximation from
  ## the maximum available degree
  valid = FALSE
  DE = seq(mt, MIN_REQUIRED_TERMS, by=-2)
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

    if (out$is.valid) {
      ## pass results into R variables
      length(out$ps.coeffs) <- out$ps.coeffs.l
      length(out$cf.coeffs) <- out$cf.coeffs.l
      length(out$offset.coeffs) <- as.integer(abs(out$diagonal.idx))
      CF <- list(out$ps.coeffs, out$cf.coeffs, out$offset.coeffs, 
                 out$diagonal.idx, out$degree)
      names(CF) <- c('ps.coeffs', 'cf.coeffs', 'offset.coeffs', 'diagonal.idx',
                     'degree')
      class(CF) <- 'CF'

      ## convert the continued fraction into a RFA
      RF <- CF2RFA(CF)
      ## roots of polynomials of numerator and denominator
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

      ## remove roots in both numerator and denominator if the real part >= 0 
      ## and the difference is less than the predefined PRECISION
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
      ## eliminate roots in numerator that are similar to the denominator
      numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
      denom.roots <- c(denom.roots.neg, denom.roots.pos)

      ## assume 1 <= t <= 100
      ## convert roots from t - 1 to t
      roots <- denom.roots + 1
      ## restrict RFA to zero as r goes to infinite
      if (length(which(Re(roots) > 0 & 
          Mod(roots) / Re(roots) / 2 <= max.extrapolation / total.sample))) {
        next
      } else {
        ## found a valid RFA
        valid = TRUE
        ## treat numerator and denomimator to be monic
        poly.numer <- as.function(poly.from.roots(numer.roots))
        poly.denom <- as.function(poly.from.roots(denom.roots))
        ## calculate the constant C
        C <- coef(RF[[1]])[length(coef(RF[[1]]))] / 
             coef(RF[[2]])[length(coef(RF[[2]]))]
        ## species accum curves with minimum count r
        ## using parital fraction expansion
        ## f is the function in terms of (t - 1)
        f <- function(t) { 
          sapply(t, function(x) {poly.numer(x) / poly.denom(x)}) * C}
        mincount.accum.curve.f <- function(x) {f(x - 1)}
        break
      }
    }
  }

  if (valid == FALSE) {
    return(NULL)
  }

  ## if the sample size is larger than max.extrapolation
  ## stop extrapolation
  ## MINOR.correction prevents machinary precision from biasing comparison
  ## result
  if (starting.size > (max.extrapolation + MINOR.correction))
  {
    index <- as.double(step.size) * (1:length(yield.estimates))
    yield.estimates <- list(sample.size = index, yields = yield.estimates)
    result <- list(continued.fraction = CF, yield.estimates = yield.estimates)
    return(result)
  }

  ## extrapolation for experiment with larger sample sizes
  start <- starting.size / total.sample
  end <- (max.extrapolation + MINOR.correction) / total.sample
  step <- step.size / total.sample

  ## extrapolating for the general accumulation curves
  extrap <- mincount.accum.curve.f(seq(start, end, by=step)) + 
            sum(n[r:length(n[, 1]), 2])

  ## combine results from interpolation/extrapolation
  yield.estimates <- c(yield.estimates, extrap)
  index <- as.double(step.size) * (1: length(yield.estimates))

  ## put index and estimated yields together into a two-colunm matrix
  yield.estimates <- matrix(c(index, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(yield.estimates) <- c('sample.size', 'yield.estimate')

  return(yield.estimates)
}


## species accum curves basis on parital fraction expansion estimation
## Rational function approximation to E(S_1(t))
## the function contains a constant c_0
### CHAO: old method
preseqR.pf.mincount.c0 <- function(n, mt = 100, ss = NULL, 
                                   max.extrapolation = NULL, r=1)
{
  checking.hist(n)
  ## setting the diagonal value
  di = 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  ## calculate total number of sample
  total.sample <- n[, 1] %*% n[, 2]
  total.sample <- floor(total.sample)

  ## set step.size as the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- floor(total.sample)
    step.size <- ss
  } else if (ss < 1) {
    write("The step size is too small", stderr())
    return()
  } else {
    step.size <- floor(ss)
  }

  ## no interpolation if step.size is larger than the size of experiment
  ## set the starting sample size as the step.size
  if (step.size > total.sample) {
    yield.estimates <- vector(mode = 'numeric', length = 0)

    ## starting sample size for extrapolation
    starting.size <- step.size
  } else {
      ## interpolation when sample size is no more than total sample size
      ## interpolate and set the size of the sample for an initial extrapolation
      out <- preseqR.interpolate.mincount(step.size, n, r)
      yield.estimates <- out[, 2]

      ## starting sample size for extrapolation
      starting.size <- ( floor(total.sample/step.size) + 1 )*step.size
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  ## transform a histogram into a vector of frequencies
  hist.count <- vector(length=max(n[, 1]), mode="numeric")
  hist.count[n[, 1]] <- n[, 2]

  ## coefficients of Good-Toulmin's power series
  PS.coeffs <- mincount.ps(hist.count, r=1)

  ## only use power series with non-zero coefficients
  ## effective coefficients from begining up to the first zero
  counts.before.first.zero = 1
  while (as.integer(counts.before.first.zero) <= length(PS.coeffs) &&
         PS.coeffs[counts.before.first.zero] != 0)
    counts.before.first.zero <- counts.before.first.zero + 1

  ## constrain the continued fraction approximation with even degree 
  ## conservatively estimates
  mt <- min(mt, counts.before.first.zero - 1)
  mt <- mt - (mt %% 2)

  ## pre-check whether sample size is sufficient
  if (mt < MIN_REQUIRED_TERMS)
  {
    m <- paste("max count before zero is les than min required count (4)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }

  PS.coeffs <- PS.coeffs[ 1:mt ]

  ## construct a continued fraction approximation with 
  ## maximum available degree
  valid = FALSE
  DE = seq(mt, MIN_REQUIRED_TERMS, by=-2)
  for (de in DE) {
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

    if (out$is.valid) {
      ## pass results into R variables
      length(out$ps.coeffs) <- out$ps.coeffs.l
      length(out$cf.coeffs) <- out$cf.coeffs.l
      length(out$offset.coeffs) <- as.integer(abs(out$diagonal.idx))
      CF <- list(out$ps.coeffs, out$cf.coeffs, out$offset.coeffs, 
                 out$diagonal.idx, out$degree)
      names(CF) <- c('ps.coeffs', 'cf.coeffs', 'offset.coeffs', 'diagonal.idx',
                     'degree')
      class(CF) <- 'CF'

      RF <- CF2RFA(CF)
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

      ## remove roots in both numerator and denominator if the real part >= 0 
      ## and the difference is less than the predefined PRECISION
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

      ## cancel similar roots in both numerator and denominator
      numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
      denom.roots <- c(denom.roots.neg, denom.roots.pos)

      ## convert roots from t - 1 to t
      roots <- denom.roots + 1
      if (length(which(roots == 0)) || length(which(Re(roots) > 0 & 
          Mod(roots) / Re(roots) / 2 <= max.extrapolation / total.sample))) {
        next
      } else {
        valid = TRUE
        ## -1 root of (t - 1) + 1
        denom.roots <- c(-1, denom.roots)
        poly.numer <- as.function(poly.from.roots(numer.roots))
        l <- length(denom.roots)
        ## treat polynomials in the rational function as monic
        ## the difference to the original rational function is a constant C
        coef <- sapply(1:l, function(x) {
          poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
        ## calculate the constant C
        C <- coef(RF[[1]])[length(coef(RF[[1]]))] / 
             coef(RF[[2]])[length(coef(RF[[2]]))]
        ## species accum curves with minimum count r
        ## using parital fraction expansion
        ## f is actually a function of "t - 1"
        f <- function(t) { 
          Re(sapply(t, function(x) {coef %*% ( (x+1)/(x-denom.roots) )^r})) * C}
        mincount.accum.curve.f <- function(x) {f(x - 1)}
        break
      }
    }
  }

  if (valid == FALSE) {
    return(NULL)
  }

  ## if the sample size is larger than max.extrapolation
  ## stop extrapolation
  ## MINOR.correction prevents machinary precision from biasing comparison
  ## result
  if (starting.size > (max.extrapolation + MINOR.correction))
  {
    index <- as.double(step.size) * (1:length(yield.estimates))
    yield.estimates <- list(sample.size = index, yields = yield.estimates)
    result <- list(continued.fraction = CF, yield.estimates = yield.estimates)
    return(result)
  }

  ## extrapolation for experiment with large sample size
  start <- starting.size / total.sample
  end <- (max.extrapolation + MINOR.correction) / total.sample
  step <- step.size / total.sample

  ## extrapolating for the general accumulation curves
  extrap <- mincount.accum.curve.f(seq(start, end, by=step)) + sum(n[, 2])

  ## combine results from interpolation/extrapolation
  yield.estimates <- c(yield.estimates, extrap)
  index <- as.double(step.size) * (1: length(yield.estimates))

  ## put index and estimated yields together into a two-colunm matrix
  yield.estimates <- matrix(c(index, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(yield.estimates) <- c('sample.size', 'yield.estimate')

  result <- list(continued.fraction = CF, estimates = yield.estimates)
  return(result)
}


## funcion of species accum curves with minimum count k based on a rational
## function approximation to E(S_1(t))
## CF is the RFA with continued fraction representation
### CHAO: old method
mincount.f <- function(CF, k) {
  RF <- CF2RFA(CF)
  numer.roots <- solve(RF[[1]])
  denom.roots <- solve(RF[[2]])
  denom.roots <- c(-1, denom.roots)

  ## record roots in the numerator that are significant different from 
  ## roots in the denominator
  tmp.roots <- c()

  ## remove roots in both numerator and denominator if the difference is
  ## less than the predefined PRECISION
  for (i in 1:length(numer.roots)) {
    d <- Mod(denom.roots - numer.roots[i])
    REMOVE.ROOTS <- FALSE
    for (j in 1:length(d)) {
      if (d[j] < PRECISION) {
        denom.roots <- denom.roots[-j]
        REMOVE.ROOTS <- TRUE
        break
      }
    }
    if (!REMOVE.ROOTS) {
      tmp.roots <- c(tmp.roots, numer.roots[i])
    }
  }

  numer.roots <- tmp.roots
  poly.numer <- as.function(poly.from.roots(numer.roots))
  l <- length(denom.roots)

  ## treat polynomials in the rational function as monic
  ## the difference is a constant C
  coef <- sapply(1:l, function(x) { poly.numer(denom.roots[x]) / 
                                    prod(denom.roots[x] - denom.roots[-x]) } )
  C <- coef(RF[[1]])[length(coef(RF[[1]]))] / coef(RF[[2]])[length(coef(RF[[2]]))]

  
  species.accum.curve.f <- function(t) {
    Re(sapply(t, function(x) {coef %*% ( (x+1)/(x-denom.roots) )^k})) * C
  }
  return(function(x) { species.accum.curve.f(x - 1) })
}


### species accumulation curve with minimum count k basis on partial fraction
### decomposition; 
### restrict RFA to be f' > 0 and f'' < 0
### CHAO: old method
preseqR.rfa.curve.derivSelect <- function(n, mt=100, ss=NULL, 
                                          max.extrapolation=NULL, k=1)
{
  checking.hist(n)
  ## setting the diagonal value
  di = 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  ## calculate total number of sample
  total.sample <- n[, 1] %*% n[, 2]
  total.sample <- floor(total.sample)

  ## set step.size as the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- floor(total.sample)
    step.size <- ss
  } else if (ss < 1) {
    write("step size is should be at least one", stderr())
    return(NULL)
  } else {
    step.size <- floor(ss)
  }

  ## no interpolation if step.size is larger than the size of experiment
  ## set the starting sample size as the step.size
  if (step.size > total.sample) {
    yield.estimates <- vector(mode = 'numeric', length = 0)

    ## starting sample size for extrapolation
    starting.size <- step.size
  } else {
      ## interpolation when sample size is no more than total sample size
      ## interpolate and set the size of sample for initial extrapolation
      out <- preseqR.interpolate.mincount(step.size, n, k)
      yield.estimates <- out[, 2]

      ## starting sample size for extrapolation
      starting.size <- ( as.integer(total.sample/step.size) + 1 )*step.size
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  ## transform a histogram into a vector of frequencies
  hist.count <- vector(length=max(n[, 1]), mode="numeric")
  hist.count[n[, 1]] <- n[, 2]
  ## only use non zeros items in histogram from begining up to the first zero
  counts.before.first.zero = 1
  while (as.integer(counts.before.first.zero) <= length(hist.count) &&
         hist.count[counts.before.first.zero] != 0)
    counts.before.first.zero <- counts.before.first.zero + 1

  ## constrain the continued fraction approximation with even degree 
  ## conservatively estimates
  mt <- min(mt, counts.before.first.zero - 1)
  mt <- mt - (mt %% 2)

  ## pre-check to make sure the sample is good for prediction
  if (mt < MIN_REQUIRED_TERMS)
  {
    m <- paste("max count before zero is les than min required count (4)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }

  ## adjust the format of count vector of the histogram in order to
  ## call the c-encoded function
  hist.count <- c(0, hist.count)

  ## allocate spaces to store constructed continued fraction approximation
  ## construct a continued fraction approximation with minimum degree
  out <- .C('c_continued_fraction_estimate', as.double(hist.count), 
            as.integer(length(hist.count)), as.integer(di), as.integer(mt),
            ps.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
            ps.coeffs.l = as.integer(0),
            cf.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
            cf.coeffs.l = as.integer(0),
            offset.coeffs =as.double(vector(mode='numeric',length=MAXLENGTH)),
            diagonal.idx = as.integer(0),
            degree = as.integer(0),
            is.valid = as.integer(0));

  if (!out$is.valid)
  {
    return(NULL)
  }

  ## restore the hist.count into the R coded format
  hist.count <- hist.count[-1]

  ## pass results into R variables
  length(out$ps.coeffs) <- out$ps.coeffs.l
  length(out$cf.coeffs) <- out$cf.coeffs.l
  length(out$offset.coeffs) <- as.integer(abs(out$diagonal.idx))
  CF <- list(out$ps.coeffs, out$cf.coeffs, out$degree)
  names(CF) <- c('ps.coeffs', 'cf.coeffs', 'degree')
  class(CF) <- 'CF'

  ## if the sample size is larger than max.extrapolation
  ## stop extrapolation
  ## MINOR.correction prevents machinary precision from biasing comparison
  ## result
  if (starting.size > (max.extrapolation + MINOR.correction))
  {
    index <- as.double(step.size) * (1:length(yield.estimates))
    yield.estimates <- list(sample.size = index, yields = yield.estimates)
    result <- list(continued.fraction = CF, yield.estimates = yield.estimates)
    return(result)
  }

  ## extrapolation for experiment with large sample size
  start <- starting.size / total.sample
  end <- (max.extrapolation + MINOR.correction) / total.sample
  step <- step.size / total.sample

  ## construct a partial fraction for extrapolating
  f <- mincount.f(CF, k)

  ## extrapolating for the general accumulation curves
  extrap <- f(seq(start, end, by=step)) + sum(n[, 2])

  ## combine results from interpolation/extrapolation
  yield.estimates <- c(yield.estimates, extrap)
  index <- as.double(step.size) * (1: length(yield.estimates))

  ## put index and estimated yields together into a two-colunm matrix
  yield.estimates <- matrix(c(index, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(yield.estimates) <- c('sample.size', 'yield.estimate')

  result <- list(continued.fraction = CF, estimates = yield.estimates)
  return(result)
}


### the power series based on count frequencies starting from frequency j
### when j = 1, it is the power series expansion of E(S_1(t)) / t at t = 1
generating.ps <- function(n, j) {
  ## transform a histogram into a vector of frequencies
  hist.count <- vector(length=max(n[, 1]), mode="numeric")
  hist.count[n[, 1]] <- n[, 2]
  
  if (j >= length(hist.count)) {
    return(NULL)
  }

  ## shift to concerned count frequencies
  hist.count <- hist.count[j: length(hist.count)]

  PS.coeffs <- sum(n[ j:dim(n)[1], 2])
  change.sign = 0

  for (i in hist.count) {
    PS.coeffs <- c(PS.coeffs, (-1)^change.sign * i - PS.coeffs[length(PS.coeffs)])
    change.sign = change.sign + 1
  }

  ## truncate at coefficients where it is zero
  zero.index = which(PS.coeffs == 0)
  if (length(zero.index) > 0) {
    PS.coeffs[1:(zero.index - 1)]
  } else {
    PS.coeffs
  }
}


### species accum curves based on parital fraction expansion estimation
### rational function approximation to E(S_1(t)) / t instead of E(S_1(t))
### CHAO: the main function
preseqR.pf.mincount <- function(n, mt = 100, ss = NULL, 
                                max.extrapolation = NULL, r=1)
{
  checking.hist(n)
  ## setting the diagonal value
  di = 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  ## total individuals
  total.sample <- n[, 1] %*% n[, 2]
  total.sample <- floor(total.sample)

  ## set step.size as the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- floor(total.sample)
    step.size <- ss
  } else if (ss < 1) {
    write("The step size is too small", stderr())
    return()
  } else {
    step.size <- floor(ss)
  }

  ## no interpolation if step.size is larger than the size of experiment
  ## set the starting sample size as the step.size
  if (step.size > total.sample) {
    yield.estimates <- vector(mode = 'numeric', length = 0)

    ## starting sample size for extrapolation
    starting.size <- step.size
  } else {
      ## interpolation when sample size is no more than total sample size
      ## interpolate and set the size of the sample for an initial extrapolation
      out <- preseqR.interpolate.mincount(step.size, n, r)
      yield.estimates <- out[, 2]

      ## starting sample size for extrapolation
      starting.size <- ( floor(total.sample/step.size) + 1 )*step.size
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  PS.coeffs <- generating.ps(n, 1)

  if (is.null(PS.coeffs)) {
    write("the size of the initial experiment is insufficient", stderr())
    return(NULL)
  }

  ## constrain the continued fraction approximation with even degree 
  ## asymptotically ~ 1 / t
  mt <- min(mt, length(PS.coeffs))
  mt <- mt - (mt %% 2)

  ## pre-check whether sample size is sufficient
  if (mt < MIN_REQUIRED_TERMS)
  {
    m <- paste("max count before zero is les than min required count (4)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }

  PS.coeffs <- PS.coeffs[ 1:mt ]

  ## construct a continued fraction approximation starting from
  ## the maximum available degree
  valid = FALSE
  DE = seq(mt, MIN_REQUIRED_TERMS, by=-2)
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

    if (out$is.valid) {
      ## pass results into R variables
      length(out$ps.coeffs) <- out$ps.coeffs.l
      length(out$cf.coeffs) <- out$cf.coeffs.l
      length(out$offset.coeffs) <- as.integer(abs(out$diagonal.idx))
      CF <- list(out$ps.coeffs, out$cf.coeffs, out$offset.coeffs, 
                 out$diagonal.idx, out$degree)
      names(CF) <- c('ps.coeffs', 'cf.coeffs', 'offset.coeffs', 'diagonal.idx',
                     'degree')
      class(CF) <- 'CF'

      ## convert the continued fraction to the RFA
      RF <- CF2RFA(CF)
      ## compatible with the constructed continued fraction in C extension
      RF[[1]] <- RF[[1]] / polynomial(c(0, 1))

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

      ## cancel roots in both numerator and denominator if the real part >= 0 
      ## and the difference is less than the predefined PRECISION
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

      ## unique roots
      numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
      denom.roots <- c(denom.roots.neg, denom.roots.pos)

      ## convert roots from t - 1 to t
      roots <- denom.roots + 1
      ## assume 1 <= t <= max range
      ## max.range = max.extrapolation / total.sample
      if (length(which(roots == 0)) || length(which(Re(roots) > 0 & 
          Mod(roots) / Re(roots) / 2 <= as.double(max.extrapolation / total.sample)))) {
        next
      } else {
        poly.numer <- as.function(poly.from.roots(numer.roots))
        l <- length(denom.roots)
        ## treat polynomials in the rational function to be monic
        ## the difference to the original RFA is a constant C
        coef <- sapply(1:l, function(x) {
          poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
        ## check whether the estimator is non-decreased
        deriv.f <- function(t) {
          Re(sapply(t, function(x) {-(coef*roots) %*% ( 1 / ((x-denom.roots)^2))}))} 
        if (length(which( deriv.f(seq(0.05, as.double(max.extrapolation / total.sample), by=0.05)) < 0 ) != 0)) {
          next
        }
        valid = TRUE
        ## calculate the constant C
        C <- coef(RF[[1]])[length(coef(RF[[1]]))] / 
             coef(RF[[2]])[length(coef(RF[[2]]))]
        ## species accum curves with minimum count r
        ## using parital fraction expansion
        ## f is actually a function of "t - 1"
        f <- function(t) { 
          Re(sapply(t, function(x) {coef %*% ( (x+1)/(x-denom.roots) )^r})) * C}
        mincount.accum.curve.f <- function(x) {f(x - 1)}
        break
      }
    }
  }

  if (valid == FALSE) {
    return(NULL)
  }

  ## if the sample size is larger than max.extrapolation
  ## stop extrapolation
  ## MINOR.correction prevents machinary precision from biasing comparison
  ## result
  if (starting.size > (max.extrapolation + MINOR.correction))
  {
    index <- as.double(step.size) * (1:length(yield.estimates))
    yield.estimates <- list(sample.size = index, yields = yield.estimates)
    result <- list(continued.fraction = CF, yield.estimates = yield.estimates)
    return(result)
  }

  ## extrapolation for experiment with large sample size
  start <- starting.size / total.sample
  end <- (max.extrapolation + MINOR.correction) / total.sample
  step <- step.size / total.sample

  ## extrapolating for the general accumulation curves
  extrap <- mincount.accum.curve.f(seq(start, end, by=step))

  ## combine results from interpolation and extrapolation
  yield.estimates <- c(yield.estimates, extrap)
  index <- as.double(step.size) * (1: length(yield.estimates))

  ## put index and estimated yields together into a two-colunm matrix
  yield.estimates <- matrix(c(index, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(yield.estimates) <- c('sample.size', 'yield.estimate')

  result <- list(continued.fraction = CF, estimates = yield.estimates)
  return(result)
}


## species accum curves based on parital fraction expansion
## using count frequencies starting from a given count frequency instead of 1
## CHAO: currently maybe not interested 
general.preseqR.pf.mincount <- function(n, mt = 100, ss = NULL, 
                                        max.extrapolation = NULL, 
                                        r=1, start.freq=1)
{
  checking.hist(n)
  ## setting the diagonal value
  di = 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  ## total individuals
  total.sample <- n[, 1] %*% n[, 2]
  total.sample <- floor(total.sample)

  ## set step.size as the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- floor(total.sample)
    step.size <- ss
  } else if (ss < 1) {
    write("The step size is too small", stderr())
    return()
  } else {
    step.size <- floor(ss)
  }

  ## no interpolation if step.size is larger than the size of experiment
  ## set the starting sample size as the step.size
  if (step.size > total.sample) {
    yield.estimates <- vector(mode = 'numeric', length = 0)

    ## starting sample size for extrapolation
    starting.size <- step.size
  } else {
      ## interpolation when sample size is no more than total sample size
      ## interpolate and set the size of the sample for an initial extrapolation
      out <- preseqR.interpolate.mincount(step.size, n, r)
      yield.estimates <- out[, 2]

      ## starting sample size for extrapolation
      starting.size <- ( floor(total.sample/step.size) + 1 )*step.size
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  PS.coeffs <- generating.ps(n, start.freq)

  if (is.null(PS.coeffs)) {
    write("the size of the initial experiment is insufficient", stderr())
    return(NULL)
  }

  ## constrain the continued fraction approximation with even degree 
  ## asymptotically ~ 1 / t
  mt <- min(mt, length(PS.coeffs))
  mt <- mt - (mt %% 2)

  ## pre-check whether the size of the initial experiment is sufficient
  if (mt < MIN_REQUIRED_TERMS)
  {
    m <- paste("max count before zero is les than min required count (4)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }

  PS.coeffs <- PS.coeffs[ 1:mt ]

  ## construct a continued fraction approximation with 
  ## maximum available degree
  valid = FALSE
  DE = seq(mt, MIN_REQUIRED_TERMS, by=-2)
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

    if (out$is.valid) {
      ## pass results into R variables
      length(out$ps.coeffs) <- out$ps.coeffs.l
      length(out$cf.coeffs) <- out$cf.coeffs.l
      length(out$offset.coeffs) <- as.integer(abs(out$diagonal.idx))
      CF <- list(out$ps.coeffs, out$cf.coeffs, out$offset.coeffs, 
                 out$diagonal.idx, out$degree)
      names(CF) <- c('ps.coeffs', 'cf.coeffs', 'offset.coeffs', 'diagonal.idx',
                     'degree')
      class(CF) <- 'CF'

      ## convert the continued fraction to the RFA
      RF <- CF2RFA(CF)
      ## compatible with constructing continued fraction in C extension
      RF[[1]] <- RF[[1]] / polynomial(c(0, 1))

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

      ## cancel roots in both numerator and denominator if the real part >= 0 
      ## and the difference is less than the predefined PRECISION
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

      ## simplified roots
      numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
      denom.roots <- c(denom.roots.neg, denom.roots.pos)

      ## convert roots from t - 1 to t for checking
      roots <- denom.roots + 1
      ## assume 1 <= t <= max range
      ## max.range = max.extrapolation / total.sample
      if (length(which(roots == 0)) || length(which(Re(roots) > 0 & 
          Mod(roots) / Re(roots) / 2 <= as.double(max.extrapolation / total.sample)))) {
        next
      } else {
        poly.numer <- as.function(poly.from.roots(numer.roots))
        l <- length(denom.roots)
        ## treat polynomials in the rational function to be monic
        ## the difference to the original RFA is a constant C
        coef <- sapply(1:l, function(x) {
          poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})

        ## modify the coefficients when using count frequencies not starting
        ## count frequency 1
        coef <- coef * (-denom.roots)^(start.freq - 1)
        ## check whether the estimator is non-decreased
        deriv.f <- function(t) {
          Re(sapply(t, function(x) {-(coef*roots) %*% ( 1 / ((x-denom.roots)^2))}))} 
        if (length(which( deriv.f(seq(0.05, as.double(max.extrapolation / total.sample), by=0.05)) < 0 ) != 0)) {
          next
        }
        valid = TRUE
        ## calculate the constant C
        C <- coef(RF[[1]])[length(coef(RF[[1]]))] / 
             coef(RF[[2]])[length(coef(RF[[2]]))]
        ## species accum curves with minimum count r
        ## using parital fraction expansion
        ## f is actually a function of "t -1" 
        f <- function(t) { 
          Re(sapply(t, function(x) {coef %*% ( (x+1)/(x-denom.roots) )^r})) * C}
        ## a function of "t"
        mincount.accum.curve.f <- function(x) {f(x - 1)}
        break
      }
    }
  }

  if (valid == FALSE) {
    return(NULL)
  }

  ## if the sample size is larger than max.extrapolation
  ## stop extrapolation
  ## MINOR.correction prevents machinary precision from biasing comparison
  ## result
  if (starting.size > (max.extrapolation + MINOR.correction))
  {
    index <- as.double(step.size) * (1:length(yield.estimates))
    yield.estimates <- list(sample.size = index, yields = yield.estimates)
    result <- list(continued.fraction = CF, yield.estimates = yield.estimates)
    return(result)
  }

  ## extrapolation for experiment with large sample size
  start <- starting.size / total.sample
  end <- (max.extrapolation + MINOR.correction) / total.sample
  step <- step.size / total.sample

  ## extrapolating for the general accumulation curves
  extrap <- mincount.accum.curve.f(seq(start, end, by=step))

  ## combine results from interpolation/extrapolation
  yield.estimates <- c(yield.estimates, extrap)
  index <- as.double(step.size) * (1: length(yield.estimates))

  ## put index and estimated yields together into a two-colunm matrix
  yield.estimates <- matrix(c(index, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(yield.estimates) <- c('sample.size', 'yield.estimate')

  result <- list(continued.fraction = CF, estimates = yield.estimates)
  return(result)
}
