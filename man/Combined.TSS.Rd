\name{Combined.TSS}
\alias{Combined.TSS}
\docType{data}
\title{
   Human TSS regions consist of TSS from hg19 UCSC gene and GENCODE V15
}
\description{
  This data set contains all TSS from hg19 UCSC gene annotation and hg19 GENCODE V15.
  This is used for identification of distal elements. 
}
\usage{data("Combined.TSS")}
\format{
  A GRanges object contains coordinates of combined TSS.
}
\value{
  A GRanges object contains coordinates of combined TSS.
}
\source{
GENCODE : \url{ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_15/gencode.v15.annotation.gtf.gz}
}
\references{
Yao L, Shen H, Laird PW, Farnham PJ,Berman BP: Inferring Regulatory Element 
Landscapes and Transcription Factor Networks from Cancer Methylomes. in revision of
Genome Biology
}
\examples{
data("Combined.TSS")
Combined.TSS
}
\keyword{datasets}
