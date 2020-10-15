#'HTClist2gi.list
#'
#'This function converts a HTClist instance into a gi_list instance with counts
#'for further use with this package, HiCDCPlus
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param htc_list A valid HTClist instance (see \code{vignette("HiTC")})
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to chromosomes in \code{htc_list}.
#'@param Dthreshold maximum distance (included) to check for significant
#'interactions, defaults to 2e6 or maximum in the data; whichever is smaller.
#'@return a thresholded gi_list instance with intra-chromosomal counts for
#'further use with HiCDCPlus
#'@examples gi_list<-generate.binned.gi.list(50e3,chrs=c('chr22'))
#'gi_list<-add.hic.counts(gi_list,
#'hic_path=system.file("extdata", "GSE63525_HMEC_combined_example.hic",
#' package = "HiCDCPlus"))
#'htc_list<-gi.list2HTClist(gi_list)
#'gi_list2<-HTClist2gi.list(htc_list,Dthreshold=Inf)
#'@export

HTClist2gi.list<-function(htc_list,chrs=NULL,Dthreshold=2e6){
  #convert htc_list back to gi_list
  if (!is.null(chrs)){
    filter<-paste0(chrs,chrs)
  }else{
    filter<-names(htc_list)
  }
  htc_list<-htc_list[HiTC::isIntraChrom(htc_list)&
                       names(htc_list)%in%filter]
  gi_list_names<-gsub("\\b(\\S+?)\\1\\S*\\b","\\1",names(htc_list),perl=TRUE)
  htcexp2gi<-function(htcexp,Dthreshold=2e6){
    anchors<-htcexp@xgi
    df<-tryCatch({as.data.frame(Matrix::summary(methods::as(htcexp@intdata,'dgTMatrix')))},
                 error=function(e){as.data.frame(Matrix::Matrix(htcexp@intdata,sparse=TRUE))})
    gi<-InteractionSet::GInteractions(anchors[df$i],anchors[df$j])
    mcols(gi)$counts<-df$x
    mcols(gi)$D<-InteractionSet::pairdist(gi)
    gi<-gi[mcols(gi)$D<=Dthreshold]
    return(gi)
  }
  gi_list<-lapply(htc_list,htcexp2gi)
  names(gi_list)<-gi_list_names
  return(gi_list)
}
  