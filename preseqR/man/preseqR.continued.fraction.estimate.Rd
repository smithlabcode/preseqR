\name{preseqR.continued.fraction.estimate}
\alias{preseqR.continued.fraction.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to predict the number of distinct items
}
\description{
preseqR.continued.fraction.estimate constructs a continued fraction to estimate 
the number of distinct items, which are not found in an intial experiment, will
be covered if an additional experiment conducts.}
\usage{
preseqR.continued.fraction.estimate(hist, di = 0, mt = 100, ss = 1e+06, 
mv = 1e+10, max.extrapolation = 1e+10, step.adjust = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist}{
	  The histogram of an experiment. It can be either the file name of a 
	  histogram or its count vector. For a histogram file, it should have two 
	  columns. Each row (x, y) represents there are y distinct items, each of
	  which occurs x times in the experiment. Both x, y are positive integers
	  and The first columns must be an increasing order. 
	  
	  A count vector is a non-negative integral vector. The ith coordinate x in
	  the vector represents that x distinct items appear in the experiment, each
	  of which occurs i times.
}
  \item{di}{
	  An integer represents diagonal value for a constructed continued fraction.
}
  \item{mt}{
	  An integer represents maximum number of parameters in a continued fraction.
}
  \item{ss}{
	  A positive real number represents the step size. It defines the size of 
	  an additional experiment.
}
  \item{mv}{
	  A positive real number represents the maximum total size of an experiment
	  that the constructed continued fraction can apply. The total size of an 
	  experiment is calculated by the size of an initial experiment plus the 
	  size of additional experiment.
}
  \item{max.extrapolation}{
	  A positive real number represents the upperbound for extrapolation. In
	  other word, it defines the maximum possible size of experiment, in which 
	  the constructed continued fraction would be called to estimate the total
	  number of distinct items.   
}
  \item{step.adjust}{
	  A logic value represents whether or not the step size could change. If
	  step.adjust equals TRUE, the value of step size (ss) could adjust in 
	  respond to the size of an initial experiment.
  }
}
\details{
Items could be anything depending on the experiment. For example, items could be
reads or molecules if it is an DNA or RNA sequencing experiment. Items could be
species if it is a population experiment. Generally, it applies to all 
experiments, which follows  capture-recapture model, although performances may
vary. 
}
\value{
  \item{continued.fraction}{
	  A list contains five components of a continued fraction: ps.coeffs is the 
	  coeffiencies of a power series, which the continued fraction tries to 
	  estimates; cf.coeffs is coeffiencies of the continued fraction; 
	  offset.coeffs is offset coeffiencies of the continued fraction; di is the 
	  diagonal value of the continued fraction; de is the degree of the 
	  continued fracton.
  }
  \item{yield.estimates}{
	  A two-columns table represents estimated number of distinct items given 
	  various sizes of experiments.
  }
}
\references{
	\url{http://smithlabresearch.org/software/preseq/}
}
\author{
	Chao Deng
}
\note{
	  The format of a constructed continued fraction is totally defined by its 
	  coefficients, its diagonal value and offset coefficients. For example,
	  when the diagonal value is zero, there are no offset coefficients. The 
	  format of a continued fraction should like this:
	  
	  \deqn{\frac{a_0}{1 + \frac{a_1}{1 + frac{a_2}{\ldots}}}}

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	  See Baker, G. \& Graves-Morris, P. \emph{Pade Approximants} about the
	  form of a continued fraction
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## import data
data(SRR1041830_5M_hist)
## hist.count is a count vector of a histogram
hist.count = vector(mode="numeric", length=max(SRR1041830_5M_hist$frequencies))
hist.count[SRR1041830_5M_hist$frequencies] = SRR1041830_5M_hist$distinct.reads
result = preseqR.continued.fraction.estimate(hist.count)
## print the constructed continued fraction
result$continued.fraction
## print the estimated yields given various sequencing sample size
yields = result$yield.estimates
names(yields)
matrix(as.numeric(unlist(yields)), ncol = 2, byrow = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ approximation }
\keyword{ continued}
\keyword{ fraction}
