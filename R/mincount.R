### calculate expected number of species appearing at least k times
### (k >= 2) through non-parameteric emprical estimator

### calculate the coefficients of the power series with freq >= k given the freq
### count of a histogram
mincount.ps <- function(hist.count, k) {
  max.freq <- length(hist.count)
  res <- vector(length = max.freq, mode = 'numeric')
  for (i in 1:max.freq) {
    co.eff <- 0
    for ( l in 0:(k-1) ) {
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


### count the unique items with freq >= k
mincount.distinct <- function(sample, k)
{
  max.value <- max(sample)
  sample.table <- vector(mode = "numeric", length = max.value)
  for (i in sample) {
    sample.table[i] <- sample.table[i] + 1
  }
  return(length(which(sample.table >= k)))
}

### interpolate when the sample size is no more than the size of
##  the initial experiment
preseqR.mincount.interpolate <- function(hist.count, ss, k)
{
  ## calculate total number of sample
  freq <- 1:length(hist.count)
  total.sample <- freq %*% hist.count

  inital.distint <- sum(hist.count)
  upper.limit <- as.integer(total.sample)
  step.size <- ss

  ## l is the number of interpolation points
  l <- as.integer(upper.limit / step.size)

  ## if the sample size is larger than the size of experiment, return NULL
  if (l == 0)
    return()

  ## sample size vector
  x <- step.size * ( 1:l )

  ## dimesion must be defined in order to use R apply
  dim(x) <- length(x)

  ## do sampling without replacement 
  s <- lapply(x, function(x) nonreplace.sampling(x, hist.count))

  ## calculate the number of distinct reads with freq k>=k based on each sample sizes
  dim(s) <- length(s)
  yield.estimates <- sapply(s, function(x) mincount.distinct(x, k))

  ## put sample.size and yield.estimates together into a matrix
  result <- matrix(c(x, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'interpolation')

  return(result)
}

### extrapolate the mincount continued fraction function
### given a histogram and a continued fraction 
preseqR.mincount.extrapolate <- function(hist.count, CF, k, start.size = NULL,
                                         step.size = NULL, max.size = NULL)
{
  ## check CF is a continued fraction with CF attribute
  if (class(CF) != "RFA")
    return(NULL)

  ## parameters for calling the c-encode function c_extrapolate_distinct
  cf.coeffs <- as.double(CF$cf.coeffs)
  cf.coeffs.l <- as.integer(length(CF$cf.coeffs))
  offset.coeffs <- as.double(CF$offset.coeffs)
  di <- as.integer(CF$diagonal.idx)
  de <- as.integer(CF$degree)
  hist.count <- as.double(hist.count)

  ## the styles of the histogram count vector are different between R code
  ## and c++ code; The first line of the histogram is always [0  0] in c++
  ## but the line is removed in R-encoded function
  hist.count <- c(0, hist.count)
  hist.count.l <- as.integer(length(hist.count))

  ## record the size of the sample based on the histogram count
  total.reads = 0.0
  for (i in 1:length(hist.count))
    total.reads <- total.reads + i*as.integer(hist.count[i])

  ## set start.size, step.size, max.size if they are not defined by user
  if (is.null(start.size))
    start.size <- total.reads
  if (start.size > max.size)
  {
    write("start position has already beyond the maximum prediction", stderr())
    return(NULL)
  }
  if (is.null(step.size))
    step.size <- total.reads
  if (is.null(max.size)) {
    ## 100 is a magic number
    max.size <- 100*total.reads
  }

  ## allocate memory to store extrapolation results
  ## first "c.extrapolate.distinct" stores the observed number of distinct
  ## molecules into estimate, then it stores the extrapolation values
  ## thus the allocated memory size is 1 plus the size of extrapolation values,
  ## which is (max.size - start.size) / step.size) + 1
  extrap.size <- as.integer( (max.size - start.size)/step.size ) + 1

  out <- .C("c_extrapolate_distinct", cf.coeffs, cf.coeffs.l, offset.coeffs,
            di, de, as.double(start.size), 
            as.double(step.size), as.double(max.size), 
            estimate = as.double(vector(mode = 'numeric', extrap.size + 1)),
            estimate.l = as.integer(0));
  
  ## return to R-coded hist.count
  hist.count <- hist.count[-1]
  if (k > length(hist.count)) {
    initial_sum = 0.0
  } else {
    initial_sum = sum(hist.count[k:length(hist.count)])
  }
  extrapolation <- out$estimate[ 1:out$estimate.l ] + initial_sum

  ## sample size vector for extrapolation
  sample.size <- start.size + step.size*( (1:length(extrapolation)) - 1 )

  ## put sample.size and extrapolation results together into a matrix
  result <- matrix(c(sample.size, extrapolation), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'extrapolation')
  return(result)
}


### construct a rational function approximation given a histogram
### di = diagonal, mt = max_terms, 
### step.adjust is an indicator for whether or not to adjust step.size
preseqR.mincount.rfa.curve <- function(hist, k, di = 0, mt = 100, ss = NULL,
                              max.extrapolation = NULL, step.adjust=TRUE,
                              header = FALSE, seed = NULL)
{
  ## set seed to reproduce the results
  if ( !is.null(seed) ) set.seed(seed)

  hist.count <- read.hist(hist, header)

  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  ## calculate total number of sample
  freq <- 1:length(hist.count)
  total.sample <- freq %*% hist.count

  ## set step.size as the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- floor(total.sample)
    step.size <- ss
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

      ## adjust step.size when it is too small
      if (step.adjust == TRUE && step.size < (total.sample / 20)) {
        step.size <- max(step.size,step.size*floor(total.sample/(20*step.size)))

        ## output the adjusted step size to stderr
        m <- paste("adjust step size to", toString(step.size), '\n', sep = ' ')
        write(m, stderr())
      }

      ## interpolate and set the size of sample for initial extrapolation
      out <- preseqR.mincount.interpolate(hist.count, step.size, k)
      yield.estimates <- out[, 2]

      ## starting sample size for extrapolation
      starting.size <- ( as.integer(total.sample/step.size) + 1 )*step.size
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  PS.coeffs <- mincount.ps(hist.count, k)
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
  ## pre check whether the power series is valid at argument 2
  if(sum(PS.coeffs) < 0.0)
  {
    m <- paste("Library expected to saturate in doubling of size",
               " unable to extrapolate", sep = ',')
    write(m, stderr())
    return(NULL)
  }

  ## allocate spaces to store constructed continued fraction approximation
  ## construct a continued fraction approximation with minimum degree
  out <- .C('c_powerseries_to_cont_frac', as.integer(di), as.integer(mt),
            as.double(PS.coeffs), as.integer(length(PS.coeffs)),
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

  ## pass results into R variables
  length(out$ps.coeffs) <- out$ps.coeffs.l
  length(out$cf.coeffs) <- out$cf.coeffs.l
  length(out$offset.coeffs) <- as.integer(abs(out$diagonal.idx))
  CF <- list(out$ps.coeffs, out$cf.coeffs, out$offset.coeffs, out$diagonal.idx,
             out$degree)
  names(CF) <- c('ps.coeffs', 'cf.coeffs', 'offset.coeffs', 'diagonal.idx',
                 'degree')
  class(CF) <- 'RFA'

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
  start <- ( starting.size - total.sample )/total.sample
  end <- ( max.extrapolation + MINOR.correction - total.sample )/total.sample
  step <- step.size/total.sample
  res <- preseqR.mincount.extrapolate(hist.count, CF, k, start, step, end)

  ## combine results from interpolation/extrapolation
  yield.estimates <- c(yield.estimates, res[, 2])
  index <- as.double(step.size) * (1: length(yield.estimates))

  ## put index and estimated yields together into a two-colunm matrix
  yield.estimates <- matrix(c(index, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(yield.estimates) <- c('sample.size', 'yield.estimate')

  result <- list(continued.fraction = CF, estimates = yield.estimates)
  return(result)
}
