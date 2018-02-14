\name{SRR611492}
\alias{SRR611492}
\docType{data}
\title{Coverage histogram of a scWGS data}
\description{
The coverage histogram is based on a single-cell whole-genome
sequencing data (scWGS) through MALBAK protocol. The accession number of the raw
data is SRR1301329. Only forward reads are used to generate the coverage
histogram.
}
\references{
  Zong, C., Lu, S., Chapman, A. R., & Xie, X. S. (2012). Genome-wide detection
  of single-nucleotide and copy-number variations of a single human cell.
  Science, 338(6114), 1622-1626.
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
data(SRR611492)
}

\keyword{ datasets }
