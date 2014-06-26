\name{bootstrap.complex.curve}
\alias{bootstrap.complex.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Complexity curve
}
\description{
	The function estimates the complexity curve of a library when its histogram
	is provided. Bootstrap is add to help the approximated rational function 
	estimate the curve and make a confidence interval.
}
\usage{
bootstrap.complex.curve(hist, times = 100, di = 0, mt = 100, ss = 1e+06, 
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
  \item{times}{
     The number of resampling times as a bootstrap process.
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
     The maximum possible value to  extrapolate
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{yield.estimates}{Yields of distinct molecules given a sample size list
    as the first column}
  \item{LOWER_0.95CI}{Lower bound for a 95\% confident interval}
  \item{UPPER_0.95CI}{Upper bound for a 96\% confident interval}
}
\references{
http://smithlabresearch.org/software/preseq/
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
function (hist, times = 100, di = 0, mt = 100, ss = 1e+06, 
    mv = 1e+10, max.extrapolation = 1e+10) 
{
   	if (mode(hist) == "character") {
		hist.count = preseqR.read.hist(hist);
	}   
	else {
		hist.count = hist;
	}  
    hist.count = preseqR.read.hist(hist.file)
    total.sample = 0
    for (i in 1:length(hist.count)) total.sample <- total.sample + 
        i * hist.count[i]
    if (times == 1) {
        out <- preseqR.continued.fraction.estimate(hist.count, 
            di, mt, ss, mv, max.extrapolation)
        if (!is.null(out)) {
            return(out$yield.estimates)
        }
        else {
            return()
        }
    }
    else if (times > 1) {
		WER_0.95CI
        N = 0
        count = 0
        step.size = 0
        estimates = matrix(data = NA, nrow = max.extrapolation/ss, 
            ncol = times, byrow = FALSE)
        for (i in 1:as.integer(times)) {
            sample = preseqR.hist.sample(hist.count, as.integer(total.sample), 
                replace = TRUE)
            hist = preseqR.sample2hist.count(sample, replace = TRUE)
            out <- preseqR.continued.fraction.estimate(hist, 
                di, mt, ss, mv, max.extrapolation)
            if (!is.null(out)) {
                count <- count + 1
                N = length(out$yield.estimates$yields)
                step.size = out$step.size
                estimates[, i][1:N] = out$yield.estimates$yields;
            }
        }
        if (N == 0) {
            write("can not make prediction based on the given histogram", 
                stderr())
            return()
        }
        if (count < BOOTSTRAP.factor * times) {
            write("fail to bootstrap since the histogram is poor", 
                stderr())
            return()
        }
        index = step.size * (1:N)
        mean = apply(estimates[1:N, ], 1, mean, na.rm = TRUE)
        variance = apply(estimates[1:N, ], 1, var, na.rm = TRUE)
        n = as.vector(apply(estimates, 1, function(x) length(which(!is.na(x)))))
        n = n[1:N]
        left.interval = mean - qnorm(0.975) * sqrt(variance / n); 
        right.interval = mean + qnorm(0.975) * sqrt(variance / n); 
        yield.estimates = list(sample.size = index, yields = yield.estimates)
        result = list(yield.estimates, left.interval, right.interval);
        names(result) = c("yield.estimates", "LOWER_0.95CI", "UPPER_0.95CI");
		return(result);
   }
    else {
        write("the paramter times should be at least one", stderr())
        return()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Complexity }
\keyword{ Library }% __ONLY ONE__ keyword per line
