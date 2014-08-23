### MAXLENGTH the allocate size for storing complexity curve
MAXLENGTH <- 10000000

### MULTINOMIAL.SAMPLE.TIMES number of random vectors to draw each time
MULTINOMIAL.SAMPLE.TIMES <- 11

### MINOR.correction a very small number to correct comparison result between
### two double type numbers when precisions can bias the result
MINOR.correction <- 1e-1

### BOOTSTRAP.factor the cut off ratio of success times / total bootstrap times
BOOTSTRAP.factor <- 0.1


### read a histogram file or a variable of a two-column table
### return the count vector of the histogram
### the ith coordinate x in a count vector means there are x distinct items,
### each of which appears i times.
### Coordinates can be zero.
read.hist <- function(hist.file, header = FALSE)
{
  if (class(hist.file) == "character") {
    hist.table <- read.table(hist.file, header = header)
  } else if (!is.null(ncol(hist.file)) && ncol(hist.file) == 2) {
      if (is.numeric(hist.file[, 1]) && is.numeric(hist.file[, 2])) {
        hist.table <- hist.file
      } else {
        stop("All items in the input should be numeric values")
      }
  } else {
    stop("input should be a variable or a file of a two-column table")
  }

  ## the first column is the frequencies of observed items
  freq <- hist.table[, 1]

  ## the second column is the number of observed distinct items for each
  ## frequency
  number.items <- hist.table[, 2]

  ## check whether frequencies are at least one and the histogram is sorted
  for (i in 1:length(freq))
    if (freq[i] <= 0 || freq[i] != floor(freq[i]) || number.items[i] < 0) {
      stop("frequencies should not be positive integers!")
    } else {
        if (i > 1 && freq[i - 1] >= freq[i])
          stop("The input histogram is not sorted in increasing order")
    }

  ## hist.count is the count vector of the histogram
  hist.count <- vector(mode = 'numeric', length = max(freq))
  hist.count[freq] <- number.items

  return(hist.count)
}


### calculate the value of the continued fraction approximation CF given the
### argument t
preseqR.rfa.estimate <- function(CF, t)
{
  if (class(CF) != "RFA")
    return()

  ## call a c-encoded function c.calculate.continued.fraction
  out <- .C("c_calculate_continued_fraction",
            cf = as.double(CF$cf.coeffs),
            cf.l = as.integer(length(CF$cf.coeffs)),
            off = as.double(CF$offset.coeffs),
            di = as.integer(CF$diagonal.idx),
            de = as.integer(CF$degree),
            coordinate = as.double(t),
            result = as.double(0));

  ## return the calculated function value
  return(out$result)
}


### extrapolate given a histogram and a continued fraction
preseqR.extrapolate.distinct <- function(hist.count, CF, start.size = NULL,
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
            di, de, hist.count, hist.count.l, as.double(start.size), 
            as.double(step.size), as.double(max.size), 
            estimate = as.double(vector(mode = 'numeric', extrap.size + 1)),
            estimate.l = as.integer(0));

  extrapolation <- out$estimate[ 1:out$estimate.l ]

  ## sample size vector for extrapolation
  sample.size <- start.size + step.size*( (1:length(extrapolation)) - 1 )

  ## put sample.size and extrapolation results together into a matrix
  result <- matrix(c(sample.size, extrapolation), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'extrapolation')
  return(result)
}


### do without replacement of random sampling given a count vector of a
### histogram; size is a user defined sample size
nonreplace.sampling <- function(size, hist.count)
{
  total.sample <- 0
  i <- 1

  ## calculate total number of sample
  freq <- 1:length(hist.count)
  total.sample <- freq %*% hist.count

  ## the number of distinct items
  distinct.sample <- sum(hist.count)

  ## identities for each distinct read
  ind <- 1:as.integer(distinct.sample)

  ## the size of each read in the library
  n <- rep(freq, as.integer(hist.count))

  ## construct a sample space X 
  ## the whole library represents by its indexes. If a read presents t
  ## times in the library, its indexes presents t times in X
  X <- rep(ind, n)

  return(sample(X, size, replace = FALSE))
}


### the function samples n histograms given a count vector of the histogram
### it is based on sampling with replacement (multinomial distribution)
replace.sampling <- function(n, hist.count)
{
  ## nonzero indexes in the hist.count
  nonzero.index <- which(hist.count != 0)

  ## nonzero values in the hist.count
  nonzero.hist.count <- hist.count[nonzero.index]

  ## calculate the distinct number of sample
  distinct.sample <- sum(hist.count)

  ## returning sampling matrix
  ## each column of the matrix represent the second column of a histogram
  ## all histograms use nonzero.index as the first column
  return( rmultinom(n, as.integer(distinct.sample), nonzero.hist.count) )
}


### given a sample vector, the function counts the number of distinct molecules
count.distinct <- function(sample)
{
  max.value <- max(sample)
  sample.table <- vector(mode = "numeric", length = max.value)
  sample.table[sample] <- 1
  return(sum(sample.table))
}


### interpolate when the sample size is no more than the size of
### the initial experiment
preseqR.interpolate.distinct <- function(hist.count, ss)
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

  ## calculate the number of distinct reads based on each sample size
  dim(s) <- length(s)
  yield.estimates <- sapply(s, function(x) count.distinct(x))

  ## put sample.size and yield.estimates together into a matrix
  result <- matrix(c(x, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'interpolation')

  return(result)
}


### check the goodness of the sample based on good Good & Toulmin's model
goodtoulmin.2x.extrap <- function(hist.count)
{
  two_fold_extrap <- 0.0
  for ( i in 1:length(hist.count) )
    two_fold_extrap <- two_fold_extrap + (-1.0)^(i + 1) * hist.count[i]
  return(two_fold_extrap)
}


### construct a rational function approximation given a frequencies of count
### data
### di = diagonal, mt = max_terms, 
### step.adjust is an indicator for whether or not to adjust step.size
preseqR.rfa.curve <- function(hist, di = 0, mt = 100, ss = NULL,
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
      out <- preseqR.interpolate.distinct(hist.count, step.size)
      yield.estimates <- out[, 2]

      ## starting sample size for extrapolation
      starting.size <- ( as.integer(total.sample/step.size) + 1 )*step.size
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

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
  if(goodtoulmin.2x.extrap(hist.count) < 0.0)
  {
    m <- paste("Library expected to saturate in doubling of size",
               " unable to extrapolate", sep = ',')
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
  res <- preseqR.extrapolate.distinct(hist.count, CF, start, step, end)

  ## combine results from interpolation/extrapolation
  yield.estimates <- c(yield.estimates, res[, 2])
  index <- as.double(step.size) * (1: length(yield.estimates))

  ## put index and estimated yields together into a two-colunm matrix
  yield.estimates <- matrix(c(index, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(yield.estimates) <- c('sample.size', 'yield.estimate')

  result <- list(continued.fraction = CF, estimates = yield.estimates)
  return(result)
}


### generate complexity curve through bootstrapping the histogram
preseqR.rfa.species.accum.curve <- function(
    hist, bootstrap.times = 100, di = 0, mt = 100, ss = NULL,
    max.extrapolation = NULL, step.adjust=TRUE, header = FALSE,
    ci = 0.95, seed = NULL)
{
  ## set seed to reproduce the results
  if ( !is.null(seed) ) set.seed(seed)

  hist.count <- read.hist(hist, header)

  ## calculate the total number of sample
  freq <- 1:length(hist.count)
  total.sample <- freq %*% hist.count

  ## calculate the distinct number of sample
  distinct.sample <- sum(hist.count)

  ## set the step.size to the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- floor(total.sample)
    step.size <- ss
  } else {
    step.size <- floor(ss)
  }

  ## adjust step.size for sampling complexity curve
  if ( step.adjust == TRUE && step.size < total.sample/20 )
  {
    step.size <- max(step.size, step.size*round(total.sample / (20*step.size)))

    ## output the adjusted step size to stderr
    m <- paste("adjust step size to", toString(step.size), '\n', sep = ' ')
    write(m, stderr())
  }

  ## set the maximum extrapolation size if it is undefined
  if (is.null(max.extrapolation)) {

    ## extrapolation 100 times; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  ## nonzero indexes in the hist.count
  nonzero.index <- which(hist.count != 0)

  ## nonzero values in the hist.count
  nonzero.hist.count <- hist.count[nonzero.index]

  ## record second columns of resampled histograms
  re.hist.second.col <- matrix(data = 0, nrow = length(nonzero.index),
                               ncol = MULTINOMIAL.SAMPLE.TIMES)

  ## the number of resampling times
  counter <- 0

  yield.estimates <- vector(mode = "numeric", length = 0)

  ## upperbound of times of iterations for bootstrapping
  upper.limit <- bootstrap.times/BOOTSTRAP.factor

  f <- function(x)
  {
    ## combine nonzero.index column and the second column to build a histogram
    ## table
    hist.table <- matrix(c(nonzero.index, x), ncol = 2, byrow = FALSE)
    preseqR.rfa.curve(
        hist.table, di, mt, step.size, max.extrapolation, step.adjust = FALSE)
  }

  while (bootstrap.times > 0) {

    ## do sampling with replacement
    ## re.hist.second.col saves the second columns of each resampled histogram
    re.hist.second.col <- replace.sampling(MULTINOMIAL.SAMPLE.TIMES, hist.count)

    ## estimate for each histogram
    out <- apply(re.hist.second.col, 2, f)

    ## eliminate NULL items in results
    out[sapply(out, is.null)] <- NULL
    ## extract yields estimation from each estimation result.
    yields <- sapply(out, function(x) x$estimates[, 2])

    if ( !is.null( dim(yields) ) )
    {
      ## update sampling status
      success.times <- dim(yields)[2]
      bootstrap.times <- bootstrap.times - success.times
      yield.estimates <- cbind(yield.estimates, yields)
    }

    ## update sampling tmes
    counter <- counter + MULTINOMIAL.SAMPLE.TIMES
    if (counter > upper.limit)
      break;
  }

  ## enough successful sampling
  if (bootstrap.times <= 0) {

    ## the number of sampled points for complexity curve
    n <- dim(yield.estimates)[1]

    ## sample sizes
    index <- as.double(step.size) * ( 1:n )

    # median values are used as complexity curve
    median.estimate <- apply(yield.estimates, 1, median)
    variance <- apply(yield.estimates, 1, var)

    # confidence interval based on lognormal
    if (ci <= 0 && ci >= 1)
      ci = 0.95
    C <- exp(qnorm((1 + ci) / 2.0) * sqrt(log(1.0 + variance / (median.estimate^2))))
    left.interval <- median.estimate/C
    right.interval <- median.estimate*C

    ## combine results and output a matrix
    result <- matrix(c(index, median.estimate, left.interval, right.interval),
                    ncol = 4, byrow = FALSE)
    lower.ci = sprintf('lower.%.2fCI', ci)
    upper.ci = sprintf('uppper.%.2fCI', ci)
    colnames(result) <- c('sample.size', 'yield.estimate', lower.ci, upper.ci)
    return(result)
  } else {
      write("fail to bootstrap!", stderr())
      return(NULL)
  }
}

print.RFA <- function(x, digit = 4, ...)
{
  s <- "CONTINUED FRACTION APPROXIMATION :\n\n"

  ## print the degree of the continued fraction approximation
  s <- paste(s, "DEGREE\t", toString(x$degree), "\n\n", sep = '')

  ## print the diagonal value
  s <- paste(s, "DIAGONAL VALUE\t", toString(x$diagonal.idx), "\n\n", sep = '')

  ## print the coefficients depending on the value of diagonal value
  s <- paste(s, "COEFFICIENTS:\n\n", sep = '')

  di <- abs(x$diagonal.idx)
    
  ## the function to print a coefficient
  ## S is the set of coefficients
  ## shift is the difference between the index of S and the index of 
  ## coefficients of a continued fraction approximation
  f <- function(index, S, shift)
  {
    s <- formatC(S[index], digit, format = 'f')
    s <- paste('a_', toString(index + shift), ' = ', s, sep = '')
  }
 
  ## print offset values if any
  if (di > 0)
  {
    index <- 1:di
    dim(index) <- di

    tmp <- apply( index, 1, function(t) f(t, x$offset.coeffs, -1) )
    tmp <- paste(tmp, collapse = '\n')
    s <- paste(s, tmp, '\n', sep = '')
  }

  ## print coeffients if any
  if (length(x$cf.coeffs) > 0)
  {
    index <- 1:length(x$cf.coeffs)
    dim(index) <- length(x$cf.coeffs)

    tmp <- apply( index, 1, function(t) f(t, x$cf.coeffs, di - 1) )
    tmp <- paste(tmp, collapse = '\n')
    s <- paste(s, tmp, '\n', sep = '')
  }
  cat(s)
  invisible(x)
}
