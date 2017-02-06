
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
#'   \item Union.enhancer: A comprehensive list of genomic strong enhancers.
#'   \item DNA methylation platform manifest: from http://zwdzwd.github.io/InfiniumAnnotation
#'   \item motif.relavent.TFs TFs family from TFClass  
#' }
#'     For more information how to create these objects please read the vignette of this package with the
#'     follwing command: \code{browseVignettes("ELMER.data")}
#'
#' @docType package
#' @name ELMER.data
#' @keywords utilities
#' @examples
#' data("Union.enhancer.hg38")
#' data("Union.enhancer.hg19")
#' data("motif.relavent.TFs")
#' data("Probes.motif.hg38.EPIC")
#' data("Probes.motif.hg19.EPIC")
#' data("Probes.motif.hg38.450K")
#' data("Probes.motif.hg19.450K")
#' data("EPIC.manifest")
#' data("EPIC.manifest.hg38")
#' data("hm450.manifest")
#' data("hm450.manifest.hg38")
#' data("motif.relevant.TFs")
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

#' A list with TF that might bind to each motif
#' @docType data
#' @keywords internal
#' @name motif.relevant.TFs
#' @format A List with 641 elements
#' @examples
#' data("motif.relevant.TFs")
NULL

#' A matrix with 1 if the probe (row) has a motif (column)
#' @docType data
#' @keywords internal
#' @name Probes.motif.hg38.EPIC
#' @format A matrix with  838881 rows and 640 columns
#' @examples
#' data("Probes.motif.hg38.EPIC")
NULL

#' A matrix with 1 if the probe (row) has a motif (column)
#' @docType data
#' @keywords internal
#' @name Probes.motif.hg19.EPIC
#' @format A matrix with  838881 rows and 640 columns
#' @examples
#' data("Probes.motif.hg19.EPIC")
NULL

#' A matrix with 1 if the probe (row) has a motif (column)
#' @docType data
#' @keywords internal
#' @name Probes.motif.hg19.450K
#' @format A matrix with  466007 rows and 640 columns
#' @examples
#' data("Probes.motif.hg19.450K")
NULL

#' A matrix with 1 if the probe (row) has a motif (column)
#' @docType data
#' @keywords internal
#' @name Probes.motif.hg38.450K
#' @format A matrix with  466007 rows and 640 columns
#' @examples
#' data("Probes.motif.hg38.450K")
NULL


#' A GRanges with enhancer regions
#' @docType data
#' @keywords internal
#' @name Union.enhancer.hg19
#' @format A GRRanges with 571084 elements
#' @examples
#' data("Union.enhancer.hg19")
NULL

#' A GRanges with enhancer regions
#' @docType data
#' @keywords internal
#' @name Union.enhancer.hg38
#' @format A GRRanges with 571084 elements
#' @examples
#' data("Union.enhancer.hg38")
NULL