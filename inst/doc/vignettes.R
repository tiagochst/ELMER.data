## ---- echo = FALSE,hide=TRUE, message=FALSE,warning=FALSE------------------
library(ELMER.data)
library(DT)
library(dplyr)

## ---- eval = FALSE---------------------------------------------------------
#  devtools::install_github(repo = "tiagochst/ELMER.data")
#  library("ELMER.data")
#  library("GenomicRanges")

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  for(plat in c("450K","EPIC")) {
#    for(genome in c("hg38","hg19")) {
#      base <- "http://zwdzwd.io/InfiniumAnnotation/current/"
#      path <- file.path(base,plat,paste(plat,"hg19.manifest.rds", sep ="."))
#      if (grepl("hg38", genome)) path <- gsub("hg19","hg38",path)
#      if(plat == "EPIC") {
#        annotation <- paste0(base,"EPIC/EPIC.hg19.manifest.rds")
#      } else {
#        annotation <- paste0(base,"hm450/hm450.hg19.manifest.rds")
#      }
#      if(grepl("hg38", genome)) annotation <- gsub("hg19","hg38",annotation)
#      if(!file.exists(basename(annotation))) {
#        if(Sys.info()["sysname"] == "Windows") mode <- "wb" else  mode <- "w"
#        downloader::download(annotation, basename(annotation), mode = mode)
#      }
#    }
#  }

## ---- message=FALSE--------------------------------------------------------
data("EPIC.hg19.manifest")
as.data.frame(EPIC.hg19.manifest)[1:5,] %>% datatable(options = list(scrollX = TRUE,pageLength = 5)) 
data("EPIC.hg38.manifest")
as.data.frame(EPIC.hg38.manifest)[1:5,] %>% datatable(options = list(scrollX = TRUE,pageLength = 5)) 
data("hm450.hg19.manifest")
as.data.frame(hm450.hg19.manifest)[1:5,] %>% datatable(options = list(scrollX = TRUE,pageLength = 5)) 
data("hm450.hg38.manifest")
as.data.frame(hm450.hg38.manifest)[1:5,] %>% datatable(options = list(scrollX = TRUE,pageLength = 5)) 

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  library(xml2)
#  library(httr)
#  library(dplyr)
#  library(rvest)
#  createMotifRelevantTfs <- function(classification = "family"){
#    message("Accessing hocomoco to get last version of TFs ", classification)
#    file <- paste0(classification,".motif.relevant.TFs.rda")
#  
#    # Download from http://hocomoco.autosome.ru/human/mono
#    tf.family <- "http://hocomoco11.autosome.ru/human/mono?full=true" %>% read_html()  %>%  html_table()
#    tf.family <- tf.family[[1]]
#    # Split TF for each family, this will help us map for each motif which are the some ones in the family
#    # basicaly: for a TF get its family then get all TF in that family
#    col <- ifelse(classification == "family", "TF family","TF subfamily")
#    family <- split(tf.family,f = tf.family[[col]])
#  
#    motif.relevant.TFs <- plyr::alply(tf.family,1, function(x){
#      f <- x[[col]]
#      if(f == "") return(x$`Transcription factor`) # Case without family, we will get only the same object
#      return(unique(family[as.character(f)][[1]]$`Transcription factor`))
#    },.progress = "text")
#    #names(motif.relevant.TFs) <- tf.family$`Transcription factor`
#    names(motif.relevant.TFs) <- tf.family$Model
#    # Cleaning object
#    attr(motif.relevant.TFs,which="split_type") <- NULL
#    attr(motif.relevant.TFs,which="split_labels") <- NULL
#    return(motif.relevant.TFs)
#  }
#  TF.family <- createMotifRelevantTfs("family")
#  TF.subfamily <- createMotifRelevantTfs("subfamily")
#  save(TF.family,file = "data/TF.family.rda", compress = "xz")
#  save(TF.subfamily,file = "data/TF.subfamily.rda", compress = "xz")

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  hocomoco.table <- "http://hocomoco11.autosome.ru/human/mono?full=true" %>% read_html()  %>%  html_table()
#  hocomoco.table <- hocomoco.table[[1]]
#  save(hocomoco.table,file = "data/hocomoco.table.rda", compress = "xz")

## --------------------------------------------------------------------------
data("Probes.motif.hg19.450K")
dim(Probes.motif.hg19.450K)
str(Probes.motif.hg19.450K)

## --------------------------------------------------------------------------
data("Probes.motif.hg38.450K")
dim(Probes.motif.hg38.450K)
str(Probes.motif.hg38.450K)

## --------------------------------------------------------------------------
data("Probes.motif.hg19.EPIC")
dim(Probes.motif.hg19.EPIC)
str(Probes.motif.hg19.EPIC)

## --------------------------------------------------------------------------
data("Probes.motif.hg38.EPIC")
dim(Probes.motif.hg38.EPIC)
str(Probes.motif.hg38.EPIC)

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  getInfiniumAnnotation <- function(plat = "450K", genome = "hg38"){
#    message("Loading object: ",file)
#    newenv <- new.env()
#    if(plat == "EPIC" & genome == "hg19") data("EPIC.hg19.manifest", package = "ELMER.data",envir=newenv)
#    if(plat == "EPIC" & genome == "hg38") data("EPIC.hg38.manifest", package = "ELMER.data",envir=newenv)
#    if(plat == "450K" & genome == "hg19") data("hm450.hg19.manifest", package = "ELMER.data",envir=newenv)
#    if(plat == "450K" & genome == "hg38") data("hm450.hg38.manifest", package = "ELMER.data",envir=newenv)
#    annotation <- get(ls(newenv)[1],envir=newenv)
#    return(annotation)
#  }
#  # To find for each probe the know motif we will use HOMER software (http://homer.salk.edu/homer/)
#  # Step:
#  # 1 - get DNA methylation probes annotation with the regions
#  # 2 - Make a bed file from it
#  # 3 - Execute section: Finding Instance of Specific Motifs from http://homer.salk.edu/homer/ngs/peakMotifs.html to the HOCOMOCO TF motifs
#  # Also, As HOMER is using more RAM than the available we will split the files in to 100k probes.
#  # Obs: for each probe we create a winddow of 500 bp (-size 500) around it. This might lead to false positives, but will not have false negatives.
#  # The false posives will be removed latter with some statistical tests.
#  TFBS.motif <- "http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif"
#  if(!file.exists(basename(TFBS.motif))) downloader::download(TFBS.motif,basename(TFBS.motif))
#  for(plat in c("EPIC","450K")){
#    for(gen in c("hg38","hg19")){
#  
#      file <- paste0(plat,gen,".txt")
#      print(file)
#      if(!file.exists(file)){
#        # STEP 1
#        gr <- getInfiniumAnnotation(plat = plat,genome =  gen)
#  
#        # This will remove masked probes. They have poor quality and might be arbitrarily positioned (Wanding Zhou)
#        gr <- gr[!gr$MASK_general]
#  
#        df <- data.frame(seqnames=seqnames(gr),
#                         starts=as.integer(start(gr)),
#                         ends=end(gr),
#                         names=names(gr),
#                         scores=c(rep(".", length(gr))),
#                         strands=strand(gr))
#        step <- 10000 # nb of lines in each file. 10K was selected to not explode RAM
#        n <- nrow(df)
#        pb <- txtProgressBar(max = floor(n/step), style = 3)
#  
#        for(j in 0:floor(n/step)){
#          setTxtProgressBar(pb, j)
#          # STEP 2
#          file.aux <- paste0(plat,gen,"_",j,".bed")
#          if(!file.exists(gsub(".bed",".txt",file.aux))){
#            end <- ifelse(((j + 1) * step) > n, n,((j + 1) * step))
#            write.table(df[((j * step) + 1):end,], file = file.aux, col.names = F, quote = F,row.names = F,sep = "\t")
#  
#            # STEP 3 use -mscore to get scores
#            cmd <- paste("source ~/.bash_rc; annotatePeaks.pl" ,file.aux, gen, "-m", basename(TFBS.motif), "-size 500 -cpu 12 >", gsub(".bed",".txt",file.aux))
#            system(cmd)
#          }
#        }
#      }
#      close(pb)
#      # We will merge the results from each file into one
#      peaks <- NULL
#      pb <- txtProgressBar(max = floor(n/step), style = 3)
#      for(j in 0:floor(n/step)){
#        setTxtProgressBar(pb, j)
#        aux <-  readr::read_tsv(paste0(plat,gen,"_",j,".txt"))
#        colnames(aux)[1] <- "PeakID"
#        if(is.null(peaks)) {
#          peaks <- aux
#        } else {
#          peaks <- rbind(peaks, aux)
#        }
#      }
#      close(pb)
#      print("Writing file...")
#      readr::write_tsv(peaks,path=file,col_names = TRUE)
#      print("DONE!")
#      gc()
#    }
#  }
#  
#  getMatrix <- function(filename) {
#    motifs <- readr::read_tsv(file)
#    # From 1 to 21 we have annotations
#    matrix <- Matrix::Matrix(0, nrow = nrow(motifs), ncol = ncol(motifs) - 21 ,sparse = TRUE)
#    colnames(matrix) <- gsub(" Distance From Peak\\(sequence,strand,conservation\\)","",colnames(motifs)[-c(1:21)])
#    rownames(matrix) <- motifs$PeakID
#    matrix[!is.na(motifs[,-c(1:21)])] <- 1
#    matrix <- as(matrix, "nsparseMatrix")
#    return(matrix)
#  }
#  
#  for(plat in c("EPIC","450K")){
#    for(gen in c("hg19","hg38")){
#      file <- paste0(plat,gen,".txt")
#  
#      if(file == "450Khg19.txt"){
#        if(file.exists("Probes.motif.hg19.450K.rda")) next
#        Probes.motif.hg19.450K <- getMatrix(file)
#        save(Probes.motif.hg19.450K, file = "Probes.motif.hg19.450K.rda", compress = "xz")
#        rm(Probes.motif.hg19.450K)
#      }
#      if(file == "450Khg38.txt"){
#        if(file.exists("Probes.motif.hg38.450K.rda")) next
#        Probes.motif.hg38.450K <- getMatrix(file)
#        save(Probes.motif.hg38.450K, file = "Probes.motif.hg38.450K.rda", compress = "xz")
#        rm(Probes.motif.hg38.450K)
#      }
#  
#      if(file == "EPIChg19.txt"){
#        if(file.exists("Probes.motif.hg19.EPIC.rda")) next
#        Probes.motif.hg19.EPIC <- getMatrix(file)
#        save(Probes.motif.hg19.EPIC, file = "Probes.motif.hg19.EPIC.rda", compress = "xz")
#        rm(Probes.motif.hg19.EPIC)
#      }
#  
#      if(file == "EPIChg38.txt"){
#        if(file.exists("Probes.motif.hg38.EPIC.rda")) next
#  
#        Probes.motif.hg38.EPIC <- getMatrix(file)
#        save(Probes.motif.hg38.EPIC, file = "Probes.motif.hg38.EPIC.rda", compress = "xz")
#        rm(Probes.motif.hg38.EPIC)
#      }
#    }
#  }

## --------------------------------------------------------------------------
data("Probes.motif.hg19.450K")
as.data.frame(as.matrix(Probes.motif.hg19.450K[1:20,1:20])) %>% 
  datatable(options = list(scrollX = TRUE,pageLength = 5)) 

## ----sessionInfo-----------------------------------------------------------
sessionInfo()

