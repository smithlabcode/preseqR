\name{SRR1301329_base}
\alias{SRR1301329_base}
\docType{data}
\title{Coverage histogram of a WES data}
\description{
The coverage histogram is based on a whole-exome sequencing (WES) 
data from Simons Foundation Autism Research Initiative with the accession number
SRR1301329. Only forward reads are used to generate the coverage histogram.
}
\references{
  Simons Foundation Autism Research Initiative (\url{https://www.sfari.org/}).
}

\details{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{N_j}, the number of nucleotides in the genome covered exactly
    \eqn{j} times.
}  

\examples{
##load library
library(preseqR)

##load data
data(SRR1301329_base)
}

\keyword{ datasets }
