\name{SRR061157_k31}
\alias{SRR061157_k31}
\docType{data}
\title{\eqn{k}-mer counts of a metagenomic data}
\description{The \eqn{k}-mer counts are based on a metagenome sequencing data from
Human Microbiome Project with the accession number
SRR061157. Only forward reads are used to generate the \eqn{k}-mer counts.}
\references{
  Human Microbiome Project (\url{https://hmpdacc.org/}).
}

\details{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{N_j}, the number of 31-mers observed exactly \eqn{j} times.
}  

\examples{
##load library
library(preseqR)

##load data
data(SRR061157_k31)
}

\keyword{ datasets }
