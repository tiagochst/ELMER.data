\name{Union.enhancer}
\alias{Union.enhancer}
\docType{data}
\title{
    A comprehensive list of genomic strong enhancers
}
\description{
    Human comprehensive genomic enhancer set came from a combination of enhancers 
    from the Roadmap Epigenomics Mapping Consortium (REMC) and the Encyclopedia of 
    DNA Elements (ENCODE) Project, in which enhancers were identified using ChromHMM 
    for 98 tissues or cell lines. We used the union of genomic elements labeled as 
    EnhG1, EnhG2, EnhA1 or EnhA2 (representing intergenic and intragenic active enhancers) 
    in any of the 98 cell types, resulting in a total of 389,967 non-overlapping 
    enhancer regions. FANTOM5 enhancers (43,011) identified by eRNAs for 400 
    distinct cell types were added to this list.
}
\usage{data("Union.enhancer")}
\format{
  A GRanges object contains coordinates of this human comprehensive genomic enhancer
  set.
}
\value{
  A GRanges object contains coordinates of this human comprehensive genomic enhancer
  set.
}
\source{
REMC ChromHMM : \url{http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final}
FANTOM5 eRNA : \url{http://enhancer.binf.ku.dk/Welcome.html}
}
\references{
1. RoadmapEpigenomicsConsortium: Integrative analysis of 111 reference human epigenomes. 
Nature 2015, 19:317-330.
2. Bernstein BE, Stamatoyannopoulos JA, Costello JF, Ren B, Milosavljevic A, Meissner A, 
Kellis M, Marra MA, Beaudet AL, Ecker JR, et al: The NIH Roadmap Epigenomics Mapping 
Consortium. Nat Biotechnol 2010, 28:1045-1048.
3. ENCODE_Project_Consortium: An integrated encyclopedia of DNA elements in the human genome. 
Nature 2012, 489:57-74.
4. Andersson R, Gebhard C, Miguel-Escalada I, Hoof I, Bornholdt J, Boyd M, Chen Y,
Zhao X, Schmidl C, Suzuki T, et al: An atlas of active enhancers across human 
cell types and tissues. Nature 2014, 507:455-461.
5. Yao L, Shen H, Laird PW, Farnham PJ,Berman BP: Inferring Regulatory Element 
Landscapes and Transcription Factor Networks from Cancer Methylomes. in revision of
Genome Biology
}
\examples{
data("Union.enhancer")
Union.enhancer
}
\keyword{datasets}
