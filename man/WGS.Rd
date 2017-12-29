\name{Whole-genome sequencing datasets}
\alias{Whole-genome sequencing datasets}
\docType{data}
\title{Whole-genome sequencing datasets}
\description{Whole-genome sequencing datasets from Chen et al. (2017).
}
\references{
Chen, C., Xing, D., Tan, L., Li, H., Zhou, G., Huang, L., & Xie, X. S. (2017). 
Single-cell whole-genome analyses by Linear Amplification via Transposon
Insertion (LIANTI). Science, 356(6334), 189-194.
}

\details{
    A two-column matrix.  
    The first column is the frequency \eqn{j = 1,2,\dots}; and the second column
    is \eqn{n_j}, the number of unique words appeared \eqn{j}
    times in Shakespeare's work. The original raw data is from Chen et al (2017).
    The raw data is then aligned to the human genome to obtain the coverage
    histograms.
}

\examples{
##load library
library(preseqR)

##load data
data(WGS)
}

\keyword{ data }
