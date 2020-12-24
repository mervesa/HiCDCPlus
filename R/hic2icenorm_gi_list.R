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
hic2htcexp<-function(chrom1,chrom2,binsize,hic_path,gen,gen_ver){
chr_select1 <- gsub("chr", "", chrom1)
chr_select2 <- gsub("chr", "", chrom2)
if (!.Platform$OS.type=="windows"){
count_matrix <- tryCatch(
  straw(norm = "NONE", fn = path.expand(hic_path), bs = binsize, ch1 = chr_select1, ch2 = chr_select2, 
        u = "BP"),
  error=function(e){
    tryCatch(straw(norm = "NONE", fn = path.expand(hic_path), bs = binsize, ch1 = chrom1, ch2 = chrom2, 
                   u = "BP"),
             error=function(e){
               straw_dump(norm = "NONE",fn=path.expand(hic_path),bs=binsize,ch1=chr_select1,ch2=chr_select2,u="BP")   
             })
  })
}else{
    count_matrix<-straw_dump(norm = "NONE",fn=path.expand(hic_path),bs=binsize,ch1=chr_select1,ch2=chr_select2,u="BP")   
}
xgi<-generate_binned_gi_list(binsize,chrs=chrom1,gen=gen,gen_ver=gen_ver)[[chrom1]]@regions
ygi<-generate_binned_gi_list(binsize,chrs=chrom2,gen=gen,gen_ver=gen_ver)[[chrom2]]@regions
names(xgi)<-paste0(GenomicRanges::seqnames(xgi),':',
                   GenomicRanges::start(xgi),'-',
                   GenomicRanges::end(xgi))
names(ygi)<-paste0(GenomicRanges::seqnames(ygi),':',
                   GenomicRanges::start(ygi),'-',
                   GenomicRanges::end(ygi))
minrow1<-min(GenomicRanges::start(xgi))
minrow2<-min(GenomicRanges::start(ygi))
maxrow1<-max(GenomicRanges::start(xgi))
maxrow2<-max(GenomicRanges::start(ygi))
count_matrix<-count_matrix%>%dplyr::filter(.data$x<=maxrow1&.data$y<=maxrow2)
if (sum(count_matrix$x==maxrow1&count_matrix$y==maxrow2)<=0){
  count_matrix<-dplyr::bind_rows(count_matrix,data.frame(x=maxrow1,y=maxrow2,
                                                         counts=0))
}
if (sum(count_matrix$x==minrow1&count_matrix$y==minrow2)<=0){
  count_matrix<-dplyr::bind_rows(count_matrix,data.frame(x=minrow1,y=minrow2,
                                                         counts=0))
}
count_matrix$x<-subjectHits(GenomicRanges::findOverlaps(
  GenomicRanges::GRanges(chrom1,
                         IRanges::IRanges(count_matrix$x,count_matrix$x+1)),xgi,minoverlap=2))
count_matrix$y<-subjectHits(GenomicRanges::findOverlaps(
  GenomicRanges::GRanges(chrom2,
                         IRanges::IRanges(count_matrix$y,count_matrix$y+1)),ygi,minoverlap=2))
count_matrix<-count_matrix%>%dplyr::arrange(.data$x,.data$y)
htc<-HiTC::normICE(HiTC::forceSymmetric(HiTC::HTCexp(
  intdata=Matrix::sparseMatrix(
    j=count_matrix$x,
    i=count_matrix$y,
    x=count_matrix$counts,symmetric = TRUE),
    ygi=ygi,
    xgi=xgi)),eps=1e-2, sparse.filter=0.025, max_iter=10)
return(htc)
}

htc_list<-list()
for (chr in chrs){
  htc_list[[paste0(chr,chr)]]<-hic2htcexp(chr,chr,binsize,hic_path,gen,gen_ver)
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
