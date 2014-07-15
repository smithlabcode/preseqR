\name{preseqR.interpolate.distinct}
\alias{preseqR.interpolate.distinct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to estimate the number of distinct items
}
\description{
preseqR.interpolate.distinct estimates the number of distinct items when small
size experiments conduct.
}
\usage{
preseqR.interpolate.distinct(hist.count, ss)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist.count}{
	  The count vector of a histogram. A count vector is a non-negative integral
	  vector. The ith coordinate x in the vector represents that x distinct 
	  items appear in the experiment, each of which occurs i times.
}
  \item{ss}{
	  A positive real number represents the step size. It defines the size of 
	  an additional experiment.
}
}

\details{
preseqR.interpolate.distinct estimate the number of distinct items given various
sizes of experiments. All experiments are conducted under the same condition as
the initial experiment except their sizes must be less than the size of an 
initial experiment. The function simulates small size experiments through 
sampling the initial experiment without replacement. Then it counts the number 
of distinct items in the sample as estimation.
}
\value{
  \item{yield.estimates}{
	  A two-columns table represents estimated number of distinct items given 
	  various sizes of experiments.
  }
}
\author{
	Chao Deng
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## import data
data(SRR1041830_5M_hist)
## hist.count is a count vector of a histogram
hist.count = SRR1041830_5M_hist$distinct.reads
## interpolate with step size 1e6
result = preseqR.interpolate.distinct(hist.count, ss = 1e6)
## print returning value
names(result)
matrix(as.numeric(unlist(result)), ncol = 2, byrow = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ interpolate }
\keyword{ nonreplacemnt }
\keyword{ sampling }
