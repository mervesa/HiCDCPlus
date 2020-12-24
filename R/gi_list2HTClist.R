#'gi_list2HTClist
#'
#'This function converts a gi_list instance into a HTClist instance
#'compatible for use with the R Bioconductor package HiTC
#'https://bioconductor.org/packages/HiTC/
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param gi_list List of \code{GenomicInteractions} objects with a counts
#'column where each object
#'named with chromosomes contains intra-chromosomal interaction information
#'(minimally containing counts and genomic distance in \code{mcols(gi_list)}---
#'see \code{?gi_list_validate} for a detailed explanation of valid
#'\code{gi_list} instances).
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to chromosomes in \code{gi_list}.
#'@return a HTClist instance compatible for use with HiTC
#'@examples gi_list<-generate_binned_gi_list(50e3,chrs=c('chr22'))
#'gi_list<-add_hic_counts(gi_list,
#'hic_path<-system.file("extdata", "GSE63525_HMEC_combined_example.hic",
#' package = "HiCDCPlus"))
#'htc_list<-gi_list2HTClist(gi_list)
#'@export

gi_list2HTClist<-function(gi_list,chrs=NULL){
  gi_list_validate(gi_list)
  if (is.null(chrs)){
    chrs<-names(gi_list)
  }
  gi_list<-gi_list[chrs]
  if(!all(vapply(gi_list,function(x) sum(names(mcols(x))=='counts')>0,TRUE))){
    stop("Not all chromosomes have a 'counts' column in their mcols(.)")
  }
  gi2htcexp<-function(gi){
    xgi<-gi@regions
    names(xgi)<-paste0(GenomicRanges::seqnames(xgi),':',
                       GenomicRanges::start(xgi),'-',
                       GenomicRanges::end(xgi))
    htc<-HiTC::forceSymmetric(HiTC::HTCexp(intdata=Matrix::sparseMatrix(
      i=InteractionSet::anchors(gi,type="first",id=TRUE),
      j=InteractionSet::anchors(gi,type="second",id=TRUE),
      x=mcols(gi)$counts,symmetric=TRUE),
      xgi=xgi,
      ygi=xgi))
    return(htc)
  }
  htc_list<-HiTC::HTClist(lapply(gi_list,gi2htcexp))
  return(htc_list)
}

