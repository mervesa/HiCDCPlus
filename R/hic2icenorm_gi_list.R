#'hic2icenorm_gi_list
#'
#' This function converts a .hic file into a gi_list instance with ICE 
#' normalized counts on the counts column for TAD annotation using  a copy of
#' TopDom (see ?TopDom_0.0.2) as well as an (optional) .hic file with ICE
#' normalized counts for visualization with Juicebox. This function requires
#' installing the Bioconductor package \code{HiTC}. 
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param hic_path Path to the .hic file.
#'@param binsize Desired bin size in bp (default 50000).
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to chromosomes in \code{gen} and \code{gen_ver}
#'except 'chrY' and 'chrM'.
#'@param hic_output If TRUE, a .hic file with the name
#'\code{gsub("\\.hic$","_icenorm.hic",hic_path)} is generated containing the ICE
#'normalized counts under 'NONE' normalization.
#'@param gen name of the species: e.g., default \code{'Hsapiens'}
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}
#'@param Dthreshold maximum distance (included) to check for significant
#'interactions, defaults to maximum in the data.
#'@return a thresholded gi_list instance with ICE normalized intra-chromosomal
#'counts for further use with this package, HiCDCPlus.
#'@examples hic_path<-system.file("extdata", "GSE63525_HMEC_combined_example.hic",
#' package = "HiCDCPlus")
#'gi_list=hic2icenorm_gi_list(hic_path,binsize=50e3,chrs=c('chr22'))
#'@export

hic2icenorm_gi_list<-function(hic_path,binsize=50000,chrs=NULL,hic_output=FALSE,
                              gen="Hsapiens",gen_ver="hg19",Dthreshold=Inf){
#convert .hic files to HiTC object with ICE normalized counts
if (is.null(chrs)){
  chrs<-get_chrs(gen,gen_ver)[chrs%in%c("chrY","chrM")]
}

htc_list<-list()
for (chrom in chrs){
  htc_list[[paste0(chrom,chrom)]]<-.hic2htcexp(chrom,chrom,binsize,hic_path,gen,gen_ver)
}
htc_list<-HiTC::HTClist(htc_list)
gi_list<-HTClist2gi_list(htc_list,Dthreshold=Dthreshold)
rm(htc_list)
if(hic_output){
  hicoutput_path<-path.expand(gsub("\\.hic$","_icenorm.hic",hic_path))
  hicoutput_pathdir<-gsub("/[^/]+$", "",hicoutput_path)
  if (hicoutput_pathdir==hicoutput_path){
    hicoutput_pathdir<-gsub("\\[^\\]+$", "",hicoutput_path)
  }
  if (hicoutput_pathdir==hicoutput_path){
    hicoutput_pathdir<-gsub("\\\\[^\\\\]+$", "",hicoutput_path)
  }
  if (!hicoutput_pathdir==hicoutput_path&!dir.exists(hicoutput_pathdir)){
    dir.create(hicoutput_pathdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
  }
  hicdc2hic(gi_list,
            hicfile=hicoutput_path,
            mode='raw',
            chrs=chrs,
            gen_ver=gen_ver)
}
return(gi_list)
}
