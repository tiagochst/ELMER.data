\name{human.TF}
\alias{human.TF}
\docType{data}
\title{All human transcription factors}
\description{
    A data frame contains symbols and gene ids for 1982 human trancription factors.
    This data is used for get.TFs function in \pkg{ELMER}.
}
\usage{data("human.TF")}
\format{A data frame contains symbols and gene ids for 1982 human transcription factors}
\value{A data frame contains symbols and gene ids for 1982 human transcription factors}

\source{
Ravasi T, Suzuki H, Cannistraci CV, Katayama S, Bajic VB, Tan K, Akalin A, 
Schmeier S, Kanamori-Katayama M, Bertin N, et al: An atlas of combinatorial 
transcriptional regulation in mouse and man. Cell 2010, 140:744-752.
}
\references{
1. Ravasi T, Suzuki H, Cannistraci CV, Katayama S, Bajic VB, Tan K, Akalin A, 
Schmeier S, Kanamori-Katayama M, Bertin N, et al: An atlas of combinatorial 
transcriptional regulation in mouse and man. Cell 2010, 140:744-752.

2. Yao L, Shen H, Laird PW, Farnham PJ,Berman BP: Inferring Regulatory Element 
Landscapes and Transcription Factor Networks from Cancer Methylomes. in revision of
Genome Biology
}

\examples{
data("human.TF")
dim(human.TF)
head(human.TF)
}
\keyword{datasets}
