\name{Shakespeare}
\alias{Shakespeare}
\docType{data}
\title{Shakespeare's word type frequencies}
\description{The Shakespeare's word type frequencies data was from 
Efron, B., & Thisted, R. (1976).}
\references{
Efron, B., & Thisted, R. (1976). Estimating the number of unseen species:
How many words did Shakespeare know?. Biometrika, 63(3), 435-447.
}

\details{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of unique words appeared \eqn{j}
    times in Shakespeare's work.
}

\examples{
##load library
library(preseqR)

##load data
data(Shakespeare)
}

\keyword{ datasets }
