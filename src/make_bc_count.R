library(IRanges)
library(Biostrings)
library(ggplot2)
library(SummarizedExperiment)



#rsync -av --include='*.umi.class.viral.bam.tsv.gz' --include='*/' --exclude='*' /Volumes/C1HT/Adrien-Philipp-Esther/Julien/data/2020-07_MAPseq4_MOS1_2_3_5_8_MO1.sindbis/ tmp/data/

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 1. Retreive UMI-corrected counts
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


#' Get umi corrected counts from a single .bam.tsv.gz file for each injection-site mapped barcode
#' @param bam_tsv path to the file
#' @return a vector of count for each mapped barcode found in the file. 
#'         Attribute "meta" additionaly contains summary statistics for spikes (useful for normalization)
read_single_bam_tsv <- function(bam_tsv) {
  # load the TSV file
  x <- read.table(sep="\t",bam_tsv,header=TRUE,stringsAsFactors=FALSE,quote="",comment.char="")
  
  # summary statistics by class and specific to spikes
  x$class <- factor(x$class,c("viral","unknown","spike"))
  meta <- local({
    N <- tabulate(x$class)
    names(N) <- levels(x$class)
    n <- tabulate(selfmatch(x$umi[x$class=="spike"]))
    c(N,spike_uniq_umi=sum(n==1),spike_multi_umi=sum(n>1))
  })

  # count each barcode
  n <- local({
    bc <- x$barcode_mapping[x$class=="viral" & x$barcode_mapping!=""]
    n <- tabulate(selfmatch(bc),length(bc))
    names(n) <- bc
    n[n>0]
  })
  attr(n,"meta") <- meta
  n
}

#' Get umi corrected counts from multiple .bam.tsv.gz file for all injection-site mapped barcode
#' @param paths paths to the files to read counts from
#' @param bc_suffix a character vector of the same length that path that specify a string to add at the end of the barcodes.
#'                  This is typically equal to "pup/injection-site" and is used to distinguish the same barcodes sequenced in different pup
#' @return a SummarizedExperiment object with a count of the barcodes in each file
read_multi_bam_tsv <- function(paths,bc_prefix="") {
  n <- mclapply(mc.cores=4,paths,read_single_bam_tsv)
  bc <- unlist(mapply(paste0,bc_prefix,lapply(n,names)),use.names = FALSE)
  N <- Matrix::sparseMatrix(
    i = selfmatch(bc),
    j = rep(seq_along(n),lengths(n)),
    x = unlist(n,use.names = FALSE),
    dims = c(length(bc),length(n)),
    dimnames = list(bc,NULL)
  )
  N <- N[Matrix::rowSums(N)>0,]
  N <- SummarizedExperiment(list(count=N),colData=as.data.frame(t(sapply(n,attr,"meta"))))
  
  N$path <- paths
  N
}

#' Read all tsv files and save the resulting SummarizedExperiment
make_umi_corrected_count_matrix <- function(dir) {
  pat <- "^([^_]*)_(.*).umi.class.viral.bam.tsv.gz$"
  paths <- list.files(dir,pat,recursive=TRUE,full.name=TRUE)
  N <- read_multi_bam_tsv(paths,sub("\\..*","/",basename(dirname(paths))))
  N$pup <- sub(pat,"\\1",basename(N$path))
  N$target <- sub(pat,"\\2",basename(N$path))
  N$injection <- sub("^([^_]*)_(.*)\\.umi.*$","\\2",basename(dirname(N$path)))
  N
}

  
