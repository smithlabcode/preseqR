\name{SRR1301329_1M_read}
\alias{SRR1301329_1M_read}
\docType{data}
\title{Read counts of a WES data}
\description{The read counts are based on a whole-exome sequencing (WES) data
from Simons Foundation Autism Research Initiative with the accession number
SRR1301329. One million reads are randomly sampled from the raw data to
generate the read counts.}
\references{
  Simons Foundation Autism Research Initiative (\url{https://www.sfari.org/}).
}

\details{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{N_j}, the number of reads observed exactly \eqn{j} times in the
    data.
}  

\examples{
##load library
library(preseqR)

##load data
data(SRR1301329_1M_read)
}

\keyword{ datasets }
