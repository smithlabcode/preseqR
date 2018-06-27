\name{Dickens}
\alias{Dickens}
\docType{data}
\title{Dickens' vocabulary}
\description{Words frequencies of a collection of Charles Dickens from Project Gutenberg}


\details{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{N_j}, the number of unique words appeared exactly \eqn{j}
    times in a collection of Charles Dickens.
}

\references{
\url{http://zipfr.r-forge.r-project.org/}
}

\examples{
##load library
library(preseqR)

##load data
data(Dickens)
}

\keyword{ datasets }
