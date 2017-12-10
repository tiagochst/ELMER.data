
#' @title Data for ELMER package
#' @description
#' ELMER is package using DNA methylation to 
#' identify enhancers, and correlates enhancer state with expression of nearby genes 
#' to identify one or more transcriptional targets. Transcription factor (TF) binding 
#' site analysis of enhancers is coupled with expression analysis of all TFs to 
#' infer upstream regulators. ELMER.data provide 2 necessary data for 
#' ELMER analysis:
#' \itemize{
#'   \item Probes.motif: motif occurences within -/+250bp of probe sites on HM450K/EPIC array aligned against hg19/hg38.
#'   \item DNA methylation platform manifest: from http://zwdzwd.github.io/InfiniumAnnotation
#'   \item TF.family TFs family from TFClass  
#'   \item TF.subfamily TFs subfamily from TFClass  
#' }
#'     For more information how to create these objects please read the vignette of this package with the
#'     follwing command: \code{browseVignettes("ELMER.data")}
#'
#' @docType package
#' @name ELMER.data
#' @exportPattern ^[[:alpha:]]+
#' @keywords utilities
#' @examples
#' data("Probes.motif.hg38.EPIC")
#' data("Probes.motif.hg19.EPIC")
#' data("Probes.motif.hg38.450K")
#' data("Probes.motif.hg19.450K")
#' data("EPIC.manifest")
#' data("EPIC.manifest.hg38")
#' data("hm450.manifest")
#' data("hm450.manifest.hg38")
#' data("elmer.data.example")
NULL

#' A GRanges containing hg19 annotation with suggested overall masking for EPIC platform
#' @docType data
#' @keywords internal
#' @name EPIC.manifest
#' @import GenomicRanges
#' @format A GRanges with 866895 elements
#' @examples
#' data("EPIC.manifest")
NULL

#' A GRanges containing hg38 annotation with suggested overall masking for EPIC platform
#' @docType data
#' @keywords internal
#' @name EPIC.manifest.hg38
#' @import GenomicRanges
#' @format A GRanges with 866895 elements
#' @examples
#' data("EPIC.manifest.hg38")
NULL

#' A GRanges containing hg38 annotation with suggested overall masking for hm450 platform
#' @docType data
#' @keywords internal
#' @name hm450.manifest.hg38
#' @import GenomicRanges
#' @format A GRanges with 485577 elements
#' @examples
#' data("hm450.manifest.hg38")
NULL

#' A GRanges containing hg19 annotation with suggested overall masking for hm450 platform
#' @docType data
#' @keywords internal
#' @name hm450.manifest
#' @import GenomicRanges
#' @format A GRanges with 485577 elements
#' @examples
#' data("hm450.manifest")
NULL

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
#' data("Probes.motif.hg38.EPIC")
NULL

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
#' data("Probes.motif.hg19.EPIC")
NULL

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
#' data("Probes.motif.hg19.450K")
NULL

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
#' data("Probes.motif.hg38.450K")
NULL


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
#' data("elmer.data.example")
NULL

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
#' data("elmer.data.example.promoter")
NULL

#' A matrix containing DNA methylation beta-values from TCGA
#' DNA methylation data: 1728 probes
#' This data is used in the examples of ELMER package
#' @docType data
#' @keywords internal
#' @name Meth
#' @format A MultiAssayExperiment for 268 Samples  and 1728 probes
NULL

#' A matrix containing gene expression data from TCGA
#' Gene Expression data: 3842 genes
#' This data is used in the examples of ELMER package
#' @docType data
#' @keywords internal
#' @name GeneExp
#' @format A gene expression matrix for 234 Samples and 3842 genes
NULL


#' A matrix containing the relashionship between the DNA methylation level at the probes
#' for a given motif and the gene expression of TF.
#' This data is used in the examples of ELMER package
#' @docType data
#' @keywords internal
#' @name TF.meth.cor
#' @format A matrix with 1968 rows (TFs) and 12 columns (motifs)
NULL

#' A list of 641 motifs with TF families (with similar bidings) from TFClass
#' Created with the following function from ELMER pacakge 
#' TF.family <-  createMotifRelevantTfs()
#' @docType data
#' @keywords internal
#' @name TF.family
#' @format A list of 641 motifs with TF families (with similar bidings)
NULL

#' A list of 641 motifs with TF subfamilies (with similar bidings) from TFClass
#' Created with the following function from ELMER pacakge 
#' TF.family <-  createMotifRelevantTfs("subfamily")
#' @docType data
#' @keywords internal
#' @name TF.subfamily
#' @format A list of 641 motifs with TF subfamilies (with similar bidings)
NULL
