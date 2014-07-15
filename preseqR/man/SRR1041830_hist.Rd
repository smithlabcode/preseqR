% --- Source file: SRR1041830_hist.txt --- 
\name{SRR1041830_hist} 
\alias{SRR1041830_hist}
\title{Histogram of A ChIP-Seq Experiment} 

\description{ 
SRR1041830_5M_hist is a histogram, which has 20 rows and 2 columns. The first
row is a header. For the rest, each row (x, y) represents there are y distinct 
reads, each of which occurs x times in the ChIP-Seq experiment. The raw 
sequencing data is from NCBI with accession No. SRX386200. Bowtie 2 is used to
map reads into hg19. SRR1041830_5M_hist is the histogram of all mapped reads.
} 

\usage{data(SRR1041830_hist)} 

\format{
This data contains the following variables:
    \itemize{
    \item frequencies: int
    \item distinct reads: int 
}
}
%\source{sequencing data see NCBI Accession: SRX386200}
\keyword{datasets}
\keyword{histogram}
