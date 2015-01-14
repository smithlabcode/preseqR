### calculate expected number of species appearing at least k times
### (k >= 2) through non-parameteric emprical estimator

### map the histogram into a power series
### read the calculated histogram through the power series;
### the parameter hist is the original histogram; k is the minimum required freq
get.hist.freq <- function(hist, j) {
  if (j %in% hist[, 1]) {
    return(hist[j, 2])
  } else {
    return(0)
  }
}

### calculate the coefficients of the power series with freq >= k given a histogram
mincount.ps <- function(hist, k) {
  max.freq <- max(hist[, 1])
  res <- vector(length = max.freq, mode = 'numeric')
  for (i in 1:max.freq) {
    co.eff <- 0
    for ( l in 0:(k-1) ) {
      for (j in 0:i) {
        index <- j + l
        if (i <= index) {
          co.eff <- co.eff + (-1)^(j + 1) * get.hist.freq(hist, index) *
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

