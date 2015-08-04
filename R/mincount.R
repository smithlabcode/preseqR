## Two roots are the same if the difference is less than the PRECISION
PRECISION = 1e-3

### calculate the coefficients of the power series with freq >= r given the
### two-column matrix of a histogram
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
preseqR.interpolate.mincount <- function(ss, n, r=1)
{
  hist <- n

  checking.hist(hist)

  ## calculate total number of sample
  total.sample <- hist[, 1] %*% hist[, 2]
  N <- floor(total.sample)

  initial.distinct <- sum(hist[, 2])
  ## the total individuals captured
  step.size <- as.double(ss)

  ## l is the number of interpolation points
  l <- as.integer(N / step.size)

  r <- floor(r)
  ## if the sample size is larger than the size of experiment or 
  ## the step size is too small, return NULL
  if (l == 0 || ss < 1 || r < 1)
    return()

  ## explicit calculating the expectation for sampling without replacement
  ## see K.L Heck 1975
  ## N is the size of population; n is the size of the sample;
  ## S is the number of unique species
  ## n is the size of the sub sample
  expect.distinct <- function(hist, N, n, S, r) {
    denom <- lchoose(N, n)
    p <- sapply(hist[, 1], function(x) {
           sum(exp(lchoose(N - x, n - 0:(r-1)) + lchoose(x, 0:(r-1)) - denom))})
    return(S - p %*% hist[, 2])
  }

  ## sample size vector
  x <- step.size * ( 1:l )

  ## calculate the number of distinct reads based on each sample size
  yield.estimates <- sapply(x, function(x) expect.distinct(hist, N, x, 
                            initial.distinct, r))

  ## put size and yield together into a matrix
  result <- matrix(c(x, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'interpolation')

  return(result)
}


## species accumu curve with minimum count r
## by RFA to Good-Toulmin power series
preseqR.rfa.mincount <- function(n, mt = 50, ss = NULL,
                                 max.extrapolation = NULL, r=1)
{
  hist <- n

  checking.hist(hist)
  ## setting the diagonal value
  di = 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  ## calculate total number of sample
  total.sample <- hist[, 1] %*% hist[, 2]
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
      out <- preseqR.interpolate.mincount(step.size, hist, r)
      yield.estimates <- out[, 2]

      ## starting sample size for extrapolation
      starting.size <- ( floor(total.sample/step.size) + 1 )*step.size
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  ## transform a histogram into a vector of frequencies
  hist.count <- vector(length=max(hist[, 1]), mode="numeric")
  hist.count[hist[, 1]] <- hist[, 2]
  ## only use non zeros items in histogram from begining up to the first zero
  counts.before.first.zero = 1
  while (as.integer(counts.before.first.zero) <= length(hist.count) &&
         hist.count[counts.before.first.zero] != 0)
    counts.before.first.zero <- counts.before.first.zero + 1

  ## for r > 1, the jth coefficient in the power series requires
  ## n_{i} i = j, j + 1, ..., j + r - 1 itmes 
  PS.coeffs <- mincount.ps(hist.count, r)

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

  ## pre-check to make sure the sample is good for prediction
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
  DE = seq(mt, 4, by=-2)
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
      class(CF) <- 'RFA'

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
      numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
      denom.roots <- c(denom.roots.neg, denom.roots.pos)

      ## assume 1 <= t <= 100
      ## convert roots from t - 1 to t
      roots <- denom.roots + 1
      ## restrict RFA to zero as r goes to infinite
      if (length(which(Re(roots) > 0 & Mod(roots) / Re(roots) / 2 <= 100))) {
        next
      } else {
        valid = TRUE
        poly.numer <- as.function(poly.from.roots(numer.roots))
        poly.denom <- as.function(poly.from.roots(denom.roots))
        ## calculate the constant C
        C <- coef(RF[[1]])[length(coef(RF[[1]]))] / 
             coef(RF[[2]])[length(coef(RF[[2]]))]
        ## species accum curves with minimum count r
        ## using parital fraction expansion
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

  ## extrapolation for experiment with large sample size
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
preseqR.pf.mincount <- function(n, mt = 100, ss = NULL, 
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
  ## only use non zeros items in histogram from begining up to the first zero
  counts.before.first.zero = 1
  while (as.integer(counts.before.first.zero) <= length(hist.count) &&
         hist.count[counts.before.first.zero] != 0)
    counts.before.first.zero <- counts.before.first.zero + 1

  ## for r > 1, the jth coefficient in the power series requires
  ## n_{i} i = j, j + 1, ..., j + r - 1 itmes 
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

  ## pre-check to make sure the sample is good for prediction
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
  DE = seq(mt, 4, by=-2)
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
      class(CF) <- 'RFA'

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
      numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
      denom.roots <- c(denom.roots.neg, denom.roots.pos)

      ## assume 1 <= t <= 100
      ## convert roots from t - 1 to t
      roots <- denom.roots + 1
      if (length(which(roots == 0)) || length(which(Re(roots) > 0 & 
          Mod(roots) / Re(roots) / 2 <= 100))) {
        next
      } else {
        valid = TRUE
        denom.roots <- c(-1, denom.roots)
        poly.numer <- as.function(poly.from.roots(numer.roots))
        l <- length(denom.roots)
        ## treat polynomials in the rational function as monic
        ## the difference is a constant C
        coef <- sapply(1:l, function(x) {
          poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
        ## calculate the constant C
        C <- coef(RF[[1]])[length(coef(RF[[1]]))] / 
             coef(RF[[2]])[length(coef(RF[[2]]))]
        ## species accum curves with minimum count r
        ## using parital fraction expansion
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
