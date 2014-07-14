% --- Source file: SRR1003759_5M_hist.txt --- 
\name{SRR1003759_5M_hist.txt} 
\alias{SRR1003759_5M_hist.txt} 
\title{A Histogram of Single-cell RNA-seq Experiment} 

\description{ 
SRR1003759_5M_hist is a histogram, which has 606 rows and 2 columns. Each row 
(x, y) represents there are y distinct reads, each of which occurs x times in 
the Single-cell RNA-seq experiment. The raw sequencing data is from NCBI with
accession No. SRX360948. The RNA-seq aligner STAR is used to map reads into 
hg19. Then around 5M mapped reads are sampled by Picard tool
DownsampleSam. SRR1003759_5M_hist is the histogram of the sample.
} 

\usage{SRR1003759_5M_hist.txt} 

\format{
This data contains the following variables:
    \itemize{
    \item frequencies: int
    \item distinct reads: int 
}
}
\keyword{datasets}
\keyword{histogram}
