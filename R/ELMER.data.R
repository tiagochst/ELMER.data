#' @title Data for ELMER package
#' @description
#' ELMER is package using DNA methylation to 
#' identify enhancers, and correlates enhancer state with expression of nearby genes 
#' to identify one or more transcriptional targets. Transcription factor (TF) binding 
#' site analysis of enhancers is coupled with expression analysis of all TFs to 
#' infer upstream regulators. ELMER.data provide the necessary data for 
#' ELMER analysis:
#' \itemize{
#'   \item Probes.motif: motif occurences within -/+250bp of probe sites on HM450K/EPIC array aligned against hg19/hg38.
#'   \item DNA methylation platform manifest: from http://zwdzwd.github.io/InfiniumAnnotation
#'   \item TF.family TFs family from TFClass  
#'   \item TF.subfamily TFs subfamily from TFClass  
#' }
#' For more information how to create these objects please read the vignette of this package with the
#' follwing command: \code{browseVignettes("ELMER.data")}
#' @docType package
#' @seealso \code{\link[ELMER.data]{EPIC.hg19.manifest}}, \code{\link[ELMER.data]{EPIC.hg38.manifest}}, 
#' \code{\link[ELMER.data]{hm450.hg19.manifest}}, \code{\link[ELMER.data]{hm450.hg38.manifest}}, 
#' \code{\link[ELMER.data]{Probes.motif.hg19.450K}}, \code{\link[ELMER.data]{Probes.motif.hg38.450K}}, 
#' \code{\link[ELMER.data]{Probes.motif.hg38.EPIC}}, \code{\link[ELMER.data]{Probes.motif.hg19.EPIC}}, 
#' \code{\link[ELMER.data]{Human_genes__GRCh37_p13__tss}}, \code{\link[ELMER.data]{Human_genes__GRCh37_p13}},
#' \code{\link[ELMER.data]{Human_genes__GRCh38_p12}}, \code{\link[ELMER.data]{Human_genes__GRCh38_p12__tss}},
#' \code{\link[ELMER.data]{TF.subfamily}}, \code{\link[ELMER.data]{TF.family}}, and \code{\link[ELMER.data]{hocomoco.table}}
#' @name ELMER.data
#' @exportPattern ^[[:alpha:]]+
#' @keywords utilities
#' @examples
#' # Please see the datasets
NULL

#' A GRanges containing hg19 annotation with suggested overall masking for EPIC platform
#' @docType data
#' @keywords internal
#' @name EPIC.hg19.manifest
#' @import GenomicRanges
#' @format A GRanges with 866895 elements
#' @examples
#' \dontrun{
#' data("EPIC.hg19.manifest")
#' }
"EPIC.hg19.manifest"

#' A GRanges containing hg38 annotation with suggested overall masking for EPIC platform
#' @docType data
#' @keywords internal
#' @name EPIC.hg38.manifest
#' @import GenomicRanges
#' @format A GRanges with 866895 elements
#' @examples
#' \dontrun{
#' data("EPIC.hg38.manifest")
#' }
"EPIC.hg38.manifest"

#' A GRanges containing hg38 annotation with suggested overall masking for hm450 platform
#' @docType data
#' @keywords internal
#' @name hm450.hg38.manifest
#' @import GenomicRanges
#' @format A GRanges with 485577 elements
#' @examples
#' \dontrun{
#' data("hm450.hg38.manifest")
#' }
"hm450.hg38.manifest"

#' A GRanges containing hg19 annotation with suggested overall masking for hm450 platform
#' @docType data
#' @keywords internal
#' @name hm450.hg19.manifest
#' @import GenomicRanges
#' @format A GRanges with 485577 elements
#' @examples
#' \dontrun{
#'  data("hm450.hg19.manifest")
#'}
"hm450.hg19.manifest"

#' A matrix with 1 if the probe (row) has a motif (column)
#' @description 
#' It was generated using HOMER with a p-value < 1e-4 
#' to scan a +/- 250bp region around each probe using HOmo sapiens 
#' COmprehensive MOdel COllection [http://hocomoco.autosome.ru/](HOCOMOCO) v10  position
#' weight matrices (PWMs). HOCOMOCO offers 640 PMWs each has a quality rating from A to D where 
#' A represents motifs with the highest confidence, and D motifs only weakly describe the pattern with a 
#' limited applications for quantitative analyses. By default only quality A and B will be used for the 
#' Motif enrichment analysis, but the minimun quality score can be seleected by the user.
#' (Addtional information [http://hocomoco.autosome.ru/help](Source)  
#' [http://nar.oxfordjournals.org/content/44/D1/D116.full](More information)). 
#' The DNA methylation information was retrieved from: http://zwdzwd.github.io/InfiniumAnnotation
#' For more information check the vignette.
#' @docType data
#' @keywords internal
#' @name Probes.motif.hg38.EPIC
#' @format A matrix with  838881 rows and 640 columns
#' @examples
#' \dontrun{
#'   data("Probes.motif.hg38.EPIC")
#' }
"Probes.motif.hg38.EPIC"

#' A matrix with 1 if the probe (row) has a motif (column)
#' @description 
#' It was generated using HOMER with a p-value < 1e-4 
#' to scan a +/- 250bp region around each probe using HOmo sapiens 
#' COmprehensive MOdel COllection [http://hocomoco.autosome.ru/](HOCOMOCO) v10  position
#' weight matrices (PWMs). HOCOMOCO offers 640 PMWs each has a quality rating from A to D where 
#' A represents motifs with the highest confidence, and D motifs only weakly describe the pattern with a 
#' limited applications for quantitative analyses. By default only quality A and B will be used for the 
#' Motif enrichment analysis, but the minimun quality score can be seleected by the user.
#' (Addtional information [http://hocomoco.autosome.ru/help](Source)  
#' [http://nar.oxfordjournals.org/content/44/D1/D116.full](More information)). 
#' The DNA methylation information was retrieved from: http://zwdzwd.github.io/InfiniumAnnotation
#' For more information check the vignette.
#' @docType data
#' @keywords internal
#' @name Probes.motif.hg19.EPIC
#' @format A matrix with  838881 rows and 640 columns
#' @examples
#' \dontrun{
#'  data("Probes.motif.hg19.EPIC")
#' }
"Probes.motif.hg19.EPIC"

#' A matrix with 1 if the probe (row) has a motif (column)
#' @description 
#' It was generated using HOMER with a p-value < 1e-4 
#' to scan a +/- 250bp region around each probe using HOmo sapiens 
#' COmprehensive MOdel COllection [http://hocomoco.autosome.ru/](HOCOMOCO) v10  position
#' weight matrices (PWMs). HOCOMOCO offers 640 PMWs each has a quality rating from A to D where 
#' A represents motifs with the highest confidence, and D motifs only weakly describe the pattern with a 
#' limited applications for quantitative analyses. By default only quality A and B will be used for the 
#' Motif enrichment analysis, but the minimun quality score can be selected by the user.
#' (Addtional information [http://hocomoco.autosome.ru/help](Source)  
#' [http://nar.oxfordjournals.org/content/44/D1/D116.full](More information)). 
#' The DNA methylation information was retrieved from: http://zwdzwd.github.io/InfiniumAnnotation
#' For more information check the vignette.
#' @docType data
#' @keywords internal
#' @name Probes.motif.hg19.450K
#' @format A matrix with  466007 rows and 640 columns
#' @examples
#' \dontrun{
#' data("Probes.motif.hg19.450K")
#' }
"Probes.motif.hg19.450K"

#' A matrix with 1 if the probe (row) has a motif (column)
#' @description 
#' It was generated using HOMER with a p-value < 1e-4 
#' to scan a +/- 250bp region around each probe using HOmo sapiens 
#' COmprehensive MOdel COllection [http://hocomoco.autosome.ru/](HOCOMOCO) v10  position
#' weight matrices (PWMs). HOCOMOCO offers 640 PMWs each has a quality rating from A to D where 
#' A represents motifs with the highest confidence, and D motifs only weakly describe the pattern with a 
#' limited applications for quantitative analyses. By default only quality A and B will be used for the 
#' Motif enrichment analysis, but the minimun quality score can be selected by the user.
#' (Addtional information [http://hocomoco.autosome.ru/help](Source)  
#' [http://nar.oxfordjournals.org/content/44/D1/D116.full](More information)). 
#' The DNA methylation information was retrieved from: http://zwdzwd.github.io/InfiniumAnnotation
#' For more information check the vignette.
#' @docType data
#' @keywords internal
#' @name Probes.motif.hg38.450K
#' @format A matrix with  466007 rows and 640 columns
#' @examples
#' \dontrun{
#' data("Probes.motif.hg38.450K")
#' }
"Probes.motif.hg38.450K"


#' A MultiAssayExperiment containing
#' DNA methylation data: 101 probes from platform 450K 
#' Gene Expression data: 1026 genes
#' for 234 samples from TCGA-LUSC.
#' This data is used in the examples of ELMER package
#' @docType data
#' @keywords internal
#' @name data
#' @format A MultiAssayExperiment for 234 Samples (8 normal samples, 226 Primary solid tumor)
#' @examples
#' \dontrun{
#' data("elmer.data.example")
#' }
"data"

#' A MultiAssayExperiment containing
#' DNA methylation data: 16 promoter probes from platform 450K 
#' Gene Expression data: 3808 genes
#' for 234 samples from TCGA-LUSC.
#' This data is used in the examples of ELMER package
#' @docType data
#' @keywords internal
#' @name mae.promoter
#' @format A MultiAssayExperiment for 234 Samples (8 normal samples, 226 Primary solid tumor)
#' @examples
#' \dontrun{
#' data("elmer.data.example.promoter")
#' }
"mae.promoter"

#' A matrix containing DNA methylation beta-values from TCGA
#' DNA methylation data: 1728 probes
#' This data is used in the examples of ELMER package
#' @docType data
#' @keywords internal
#' @name Meth
#' @format A MultiAssayExperiment for 268 Samples  and 1728 probes
"Meth"

#' A matrix containing gene expression data from TCGA
#' Gene Expression data: 3842 genes
#' This data is used in the examples of ELMER package
#' @docType data
#' @keywords internal
#' @name GeneExp
#' @format A gene expression matrix for 234 Samples and 3842 genes
"GeneExp"


#' A list of 641 motifs with TF families (with similar bidings) from TFClass
#' Created with the following function from ELMER pacakge 
#' TF.family <-  createMotifRelevantTfs()
#' @docType data
#' @keywords internal
#' @name TF.family
#' @format A list of 641 motifs with TF families (with similar bidings)
"TF.family"

#' A list of 641 motifs with TF subfamilies (with similar bidings) from TFClass
#' Created with the following function from ELMER pacakge 
#' TF.family <-  createMotifRelevantTfs("subfamily")
#' @docType data
#' @keywords internal
#' @name TF.subfamily
#' @format A list of 641 motifs with TF subfamilies (with similar bidings)
"TF.subfamily"


#' Table parsed from hocomoco v11
#' @docType data
#' @keywords internal
#' @name hocomoco.table
#' @format A dataframe with 771 rows (motifs) and 20 columns
"hocomoco.table"

#' Table parsed from Lambert, Samuel A., et al. "The human transcription factors." Cell 172.4 (2018): 650-665.
#' @docType data
#' @keywords internal
#' @name human.TF
#' @format A dataframe with 1639 rows (motifs) and 27 columns
"human.TF"


#' A matrix containing ENSEMBL hg19 transcripts metadata accessed using biomart
#' This data is used if ensembl cannot be reached
#' @docType data
#' @name Human_genes__GRCh37_p13__tss
#' @keywords internal
#' @format A matrix with metadata for 196317 transcripts
"Human_genes__GRCh37_p13__tss"


#' A matrix containing ENSEMBL hg38 transcripts metadata accessed using biomart
#' This data is used if ensembl cannot be reached
#' @docType data
#' @keywords internal
#' @name Human_genes__GRCh38_p12__tss
#' @format A matrix with metadata for 208423 transcripts
"Human_genes__GRCh38_p12__tss"

#' A matrix containing ENSEMBL hg19 gene metadata accessed using biomart
#' This data is used if ensembl cannot be reached
#' @docType data
#' @keywords internal
#' @name Human_genes__GRCh37_p13
#' @format A matrix with metadata for 60482 genes
"Human_genes__GRCh37_p13"

#' A matrix containing ENSEMBL hg38 gene metadata accessed using biomart
#' This data is used if ensembl cannot be reached
#' @docType data
#' @name Human_genes__GRCh38_p12
#' @keywords internal
#' @format A matrix with metadata for 58639 genes
"Human_genes__GRCh38_p12"
