\name{Twitter}
\alias{Twitter}
\docType{data}
\title{Social network}
\description{Following relationships of Twitter's social network}

\details{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of users with exactly \eqn{j} followers.
}

\references{
Zafarani R, Liu H (2009) Social computing data repository at ASU.
}

\examples{
##load library
library(preseqR)

##load data
data(Twitter)
}

\keyword{ datasets }
