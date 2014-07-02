\name{preseqR.continued.fraction.estimate}
\alias{preseqR.continued.fraction.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to predict library complexity
}
\description{
preseqR.continued.fraction.estimate creates an continued fraction to estimate 
the number of distinct molecules, reads, or species given its histogram. It also
provides a complexity curve to describe the complexity of the capture-recapture
data.
}
\usage{
preseqR.continued.fraction.estimate(hist, di = 0, mt = 100, ss = 1e+06, 
									mv = 1e+10, max.extrapolation = 1e+10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
	  A histogram. It can be either a file name or a count vector of the 
	  histogram. For histogram file, it contains two columns. The first column
	  is frequencies of molecules. Values of frequency should be at least one.
	  For each given frequency, the second column is the number of molecules
	  with that frequency.
}
  \item{di}{
	  Diagonal value for a constructed continued fraction.
}
  \item{mt}{
	  Maximum number of parameters in a continued fraction.
}
  \item{ss}{
	  Step size of sampling points along a library.
}
  \item{mv}{
	  The maximum value to train a continued fraction.
}
  \item{max.extrapolation}{
	  The maximum possible value to extrapolate
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{CF}{All components of a continued fraction. ps.coeffs is the 
	  coeffiencies of a power series, which the continued fraction estimates;
	  co.coeffs is the coeffiencies of the continued fraction; offset.coeffs is
  	  offset coeffiencies of the contined fraction; di is the diagonal value of
      the continued fraction; de is the degree of the contined fracton}
  \item{yield.estimates}{Yields of distinct molecules given a sample size list
  as the first column}
}
\references{
	\url{http://smithlabresearch.org/software/preseq/}
}
\author{
	Chao Deng
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (hist, di = 0, mt = 100, ss = 1e+06, mv = 1e+10, max.extrapolation = 1e+10) 
{
    if (mode(hist) == "character") {
        hist.count = preseqR.read.hist(hist)
    }
    else {
        hist.count = hist
    }
    MIN_REQUIRED_TERMS = 4
    total.sample = 0
    for (i in 1:length(hist.count)) total.sample <- total.sample + 
        i * hist.count[i]
    step.size = ss
    if (step.size > total.sample) {
        yield.estimates = vector(mode = "numeric", length = 0)
        starting.size = step.size
    }
    else {
        if (step.size < (total.sample/20)) {
            step.size = max(step.size, step.size * round(total.sample/(20 * 
                step.size)))
            m = paste("adjust step size to", toString(step.size), 
                "\n", sep = " ")
            write(m, stderr())
        }
        out = preseqR.interpolate.distinct(hist.count, step.size)
        yield.estimates = out$yield.estimates
        starting.size = out$sample.size
    }
    counts.before.first.zero = 1
    while (as.integer(counts.before.first.zero) <= length(hist.count) && 
        hist.count[counts.before.first.zero] != 0) 
		counts.before.first.zero <- counts.before.first.zero + 1
    mt = min(mt, counts.before.first.zero - 1)
    mt = mt - (mt\%\%2)
    if (mt < MIN_REQUIRED_TERMS) {
        m = paste("max count before zero is les than min required count (4), ", 
            "sample not sufficiently deep or duplicates removed", 
            sep = "")
        write(m, stderr())
        return()
    }
    if (goodtoulmin.2x.extrap(hist.count) < 0) {
        m = paste("Library expected to saturate in doubling of size, ", 
            "unable to extrapolate", sep = "")
        write(m, stderr())
        return()
    }
    hist.count = c(0, hist.count)
    out <- .C("c_continued_fraction_estimate", as.double(hist.count), 
        as.integer(length(hist.count)), as.integer(di), as.integer(mt), 
        step.size = as.double(step.size), as.double(mv), 
		ps.coeffs = as.double(vector(mode = "numeric", length = MAXLENGTH)), 
		ps.coeffs.l = as.integer(0), 
        cf.coeffs = as.double(vector(mode = "numeric", length = MAXLENGTH)), 
        cf.coeffs.l = as.integer(0), 
		offset.coeffs = as.double(vector(mode = "numeric", length = MAXLENGTH)),
	   	diagonal.idx = as.integer(0), 
        degree = as.integer(0), is.valid = as.integer(0))
    if (!out$is.valid) {
        write("Fail to construct and need to bootstrap to obtain estimates", 
            stderr())
        return()
    }
    length(out$ps.coeffs) = out$ps.coeffs.l
    length(out$cf.coeffs) = out$cf.coeffs.l
    length(out$offset.coeffs) = as.integer(abs(out$diagonal.idx))
    CF = list(out$ps.coeffs, out$cf.coeffs, out$offset.coeffs, 
        out$diagonal.idx, out$degree)
    names(CF) = c("ps.coeffs", "cf.coeffs", "offset.coeffs", 
        "diagonal.idx", "degree")
    if (starting.size > max.extrapolation) {
		index = as.integer(out$step.size) * (1: length(yield.estimates));
		yield.estimates = list(sample.size = index, yields = yield.estimates);
        result = list(CF, yield.estimates, out$step.size)
        names(result) = c("continued.fraction", "yield.estimates", 
            "step.size")
        return(result)
    }
    est <- preseqR.extrapolate.distinct(hist.count, CF, (starting.size - 
        total.sample)/total.sample, out$step.size/total.sample, 
        (max.extrapolation - total.sample)/total.sample)
    est = est[-1]
    yield.estimates = c(yield.estimates, est)
	index = as.integer(out$step.size) * (1: length(yield.estimates));
	yield.estimates = list(sample.size = index, yields = yield.estimates);
    result = list(CF, yield.estimates, out$step.size)
    names(result) = c("continued.fraction", "yield.estimates", 
        "step.size")
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ approximation }
