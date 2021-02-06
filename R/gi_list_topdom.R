#'gi_list_topdom
#'
#'This function converts a gi_list instance with ICE normalized counts
#'into TAD annotations through an implementation of TopDom v0.0.2
#'(https://github.com/HenrikBengtsson/TopDom) adapted as 
#'\code{TopDom} at this package. If you're using this function, please cite 
#'TopDom according to the documentation at 
#'https://github.com/HenrikBengtsson/TopDom/blob/0.0.2/docs/
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param gi_list List of \code{GenomicInteractions} objects where each object
#' named with chromosomes contains intrachromosomal interaction information
#' (see
#' \code{?gi_list_validate} for a detailed explanation of valid \code{gi_list}
#' instances). 
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to chromosomes in \code{gi_list}.
#'@param file_out If true, outputs TAD annotations into files with paths
#'beginning with \code{fpath}. Defaults to FALSE
#'@param fpath Outputs TAD annotations into files with paths beginning
#'in \code{fpath}.
#'@param window.size integer, number of bins to extend. Defaults to 5.
#'@param verbose TRUE if you would like to troubleshoot TopDom.
#'@return a list instance with TAD annotation reporting for each chromosome
#'@examples 
#'hic_path<-system.file("extdata", "GSE63525_HMEC_combined_example.hic",
#' package = "HiCDCPlus")
#'gi_list<-hic2icenorm_gi_list(hic_path,binsize=50e3,chrs='chr22')
#'tads<-gi_list_topdom(gi_list)
#'@export


gi_list_topdom<-function(gi_list,chrs=NULL,file_out=FALSE,fpath=NULL,window.size=5, verbose=FALSE){
if(is.null(chrs)){
    chrs<-names(gi_list)
}
gi_list<-gi_list[chrs]
out<-list()
for (chrom in chrs){
m<-InteractionSet::inflate(gi_list[[chrom]],
           rows=gi_list[[chrom]]@regions,
           columns=gi_list[[chrom]]@regions,
           fill=mcols(gi_list[[chrom]])$counts)@matrix
m[is.na(m)]<-0
rownames(m) <- NULL
chr<-unique(GenomicRanges::seqnames(gi_list[[chrom]]@regions))
sub_mat <- cbind.data.frame(GenomicRanges::seqnames(gi_list[[chrom]]@regions), 
                            GenomicRanges::start(gi_list[[chrom]]@regions), 
                            GenomicRanges::end(gi_list[[chrom]]@regions), m)
colnames(sub_mat) = NULL
tmpfile <- paste0(base::tempfile(), ".txt")
data.table::fwrite(sub_mat, tmpfile, sep="\t",
                   col.names=FALSE, row.names=FALSE,quote=FALSE)
outfile<-NULL
if (file_out) outfile<-path.expand(paste0(fpath, "_", chr))
out[[chrom]]<-.TopDom(matrix.file=tmpfile,
             window.size=5, outFile=outfile,verbose=verbose)
system2("rm", args = path.expand(tmpfile))
}
return(out)
}
