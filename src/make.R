



#-#-#-#-#-#-#-#-#-#
# WHICH REGION IS EACH READ COMING FROM?
# build the index mapping file that associates to any sequence an index (corresponding to a region)
#-#-#-#-#-#-#-#-#-#
make_demuxmap <- function(index_file) {
  library(IRanges)
  library(Biostrings)
  library(Biobase)
  
  index <- as(read.table(index_file,sep="\t",header=TRUE,stringsAsFactors = FALSE),"DataFrame")
  index$sequence <- DNAStringSet(gsub(" *","",index$sequence))
  index <- index[grepl("(.*)_(.*)",index$name),]
  
  # determine minimal barcode distance
  local({
    d <- as.matrix(stringDist(index$sequence))
    write.table(d,sep=" ")
    diag(d) <- +Inf
    cat("--------------------------------------------\n")
    cat("Min barcode distance:",min(rowMin(d)),"\n")
    cat("--------------------------------------------\n")
  })
  
  map <- DataFrame(seq=DNAStringSet(mkAllStrings(c("A","C","G","T","N"),6)))
  n <- sapply(reverseComplement(index$sequence),neditAt,map$seq,at=1)
  map$closest_index_id <- max.col(-n)
  map$closest_index_seq <- index$sequence[map$closest_index_id]
  map$closest_index_name <- index$name[map$closest_index_id]
  map$closest_index_mismatch <- n[cbind(seq_along(map$closest_index_id),map$closest_index_id)]
  
  map
}


#-#-#-#-#-#-#-#-#-#
# map a FASTA file onto the given bowtie index (Justus library sequences???), and generate a BAM
# The method uses Bowtie and find all possible alignments with 3 mismatches allowed
#-#-#-#-#-#-#-#-#-#
bowtie_map <- function(fa.file,bwt.index,bam.file) {
  library(Rbowtie)
  library(Rsamtools)
  
  tmp.sam <- tempfile(fileext=".sam")
  bowtie(sequences=fa.file,index=bwt.index,outfile=tmp.sam,f=TRUE,p=3,n=3,v=3,l=32,sam=TRUE,all=TRUE,best=TRUE,force=TRUE) # v = allowed mismatches; n = max mismatches in seed (0 to 3, if seed = 28 nt); p = number of CPU for analysis (keep 3) ; l = length of the seed (28 by default)
  asBam(tmp.sam,sub(".bam$","",bam.file),overwrite=TRUE)
  return(bam.file)
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Find clusters in a self-mapped BAM graph
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
make_barcode_cluster <- function(self.bam) {
  library(igraph)
  library(GenomicAlignments)
  library(Biostrings)
  library(S4Vectors)
  
  # Load the barcode graph
  aln <- readGAlignments(self.bam,param=ScanBamParam(what="qname",tag="NM",flag = scanBamFlag(isMinusStrand = FALSE)))
  g <- graph_from_edgelist(cbind(mcols(aln)$qname,as.character(seqnames(aln))),directed=FALSE)
  #E(g)$weight <- mcols(aln)$NM
  
  #-#-#-#-#-#-#-#-#-#
  # Cluster the barcodes
  # the cluster center is the most abundant sequence within a cluster, considered as the "real" sequence of this barcode
  # 3 mismatches are allowed between barcodes to define them as bellonging to a same cluster
  #-#-#-#-#-#-#-#-#-#
  k <- components(g)
  k <- DataFrame(vertex=names(k$membership),cluster_id=unname(k$membership),cluster_size=k$csize[k$membership])
  pat <- "([A-Z]*)_([0-9]*)"
  k$barcode_count <- as.integer(sub(pat,"\\2",k$vertex))
  k$barcode <- DNAStringSet(sub(pat,"\\1",k$vertex))
  o <- order(k$barcode_count,decreasing=TRUE)
  k$is_cluster_center[o] <- !duplicated(k$cluster_id[o])
  k$cluster_barcode_count <- as.vector(tapply(k$barcode_count,k$cluster_id,sum)[k$cluster_id])
  
  k
}


#-#-#-#-#-#-#-#-#-#
# Map all barcodes of a pup onto filtered S1 clusters identified
#-#-#-#-#-#-#-#-#-#
make_cluster_assignment <- function(in.tsv,in.bam,out.tsv.gz) {
  library(GenomicAlignments)
  aln <- readGAlignments(in.bam,param=ScanBamParam(what="seq",tag="NM",flag=scanBamFlag(isMinusStrand=FALSE)))
  x <- read.table(in.tsv,sep="\t",stringsAsFactors=FALSE,header=TRUE)
  i <- match(x$bc32,mcols(aln)$seq)
  x$barcode_mapping <- as.character(seqnames(aln))[i]
  x$barcode_mapping_err <- mcols(aln)$NM[i]
  con <- gzfile(out.tsv.gz,"w")
  on.exit(close(con))
  write.table(x,file=con,sep="\t",quote=FALSE,row.names=FALSE,na="")
}





