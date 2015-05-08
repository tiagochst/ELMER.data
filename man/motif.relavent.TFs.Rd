\name{motif.relavent.TFs}
\alias{motif.relavent.TFs}
\docType{data}
\title{
  TFs may recognize the same motif.
}
\description{
  The transcription factor from the same family may recognize the similar motif. We
  collect the TFs that recognize a similar motif for each motif.
}
\usage{data("motif.relavent.TFs")}
\format{
  A list with names (motif) and contents (TF which recognize a similar motif) for
  each list.
}
\value{
  A list with names (motif) and contents (TF which recognize a similar motif) for
  each list.
}
\references{
Yao L, Shen H, Laird PW, Farnham PJ,Berman BP: Inferring Regulatory Element 
Landscapes and Transcription Factor Networks from Cancer Methylomes. in revision of
Genome Biology
}
\examples{
data("motif.relavent.TFs")
str(motif.relavent.TFs)
}
\keyword{datasets}
