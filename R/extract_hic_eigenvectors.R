#' extract_hic_eigenvectors
#'
#' This function uses Juicer command line tools to extract first eigenvectors
#' across chromosomes from counts data in a .hic file and outputs them to text 
#' file of the structure chr start end score where the score column contains 
#' the eigenvector elements.
#'@import BSgenome splines
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param hicfile path to the input .hic file.
#'@param mode Normalization mode to extract first eigenvectors from 
#'Allowable options are: 'NONE' for raw (normalized counts if .hic file is 
#'written using \code{hicdc2hic} or \code{hic2icenorm_gi_list}), 'KR' for
#'Knight-Ruiz normalization, 'VC' for Vanilla-Coverage normalization and
#''VC_SQRT' for square root vanilla coverage. Defaults to 'KR'.
#'@param binsize the uniform binning size for compartment scores in bp.
#' Defaults to 100e3.
#'@param chrs a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes except "Y", and "M" for
#'the specified \code{gen} and \code{gen_ver}.
#'@param gen name of the species: e.g., default \code{'Hsapiens'}.
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}.
#'@return path to the eigenvector text files for each chromosome containing
#'chromosome, start, end
#'and compartment score values that may need to be flipped signs for each
#'chromosome. File paths 
#'follow \code{gsub('.hic','_<chromosome>_eigenvectors.txt',hicfile)}
#'@examples eigenvector_filepaths<-extract_hic_eigenvectors(
#'hicfile=system.file("extdata", "GSE63525_HMEC_combined_example.hic",
#'package = "HiCDCPlus"),
#'chrs=c("chr22"),binsize=50e3)
#'@export

extract_hic_eigenvectors<-function(hicfile,mode='KR',binsize=100e3,
                                   chrs=NULL,gen="Hsapiens",gen_ver="hg19"){
    options(scipen=9999)
    #set memory limit to max if i386
    if (.Platform$OS.type=='windows'&Sys.getenv("R_ARCH")=="/i386") {
        gc(reset=TRUE,full=TRUE)
        utils::memory.limit(size=4095)
    }
    filepaths<-c()
    if (is.null(chrs)){chrs<-get_chrs(gen,gen_ver)}
    jarpath<-.download_juicer()
    for(chrom in chrs){
    chr_select<-gsub("chr","",chrom)
    tmpfile <- paste0(base::tempfile(), ".txt")
    tryCatch(system2("java", args = c(ifelse(.Platform$OS.type=='windows'&Sys.getenv("R_ARCH")=="/i386","-Xmx2g","-Xmx8g"), "-jar", path.expand(jarpath), 
                             "eigenvector", "-p", mode, path.expand(hicfile),
                             chr_select, "BP", binsize, path.expand(tmpfile))),
             error=function(e){
            system2("java", args = c(ifelse(.Platform$OS.type=='windows'&Sys.getenv("R_ARCH")=="/i386","-Xmx2g","-Xmx8g"), "-jar", path.expand(jarpath), 
                            "eigenvector", "-p", mode, path.expand(hicfile),
                             chrom, "BP", binsize, path.expand(tmpfile)))                 
             })
    out.df<-data.table::fread(tmpfile)
    system2("rm", args = path.expand(tmpfile))
    out.df<-out.df%>%dplyr::mutate(chr=chrom,
                            start=binsize*seq(0,nrow(out.df)-1,1),
                            end=binsize*seq(1,nrow(out.df),1))%>%
        dplyr::rename("score"=.data$V1)%>%dplyr::select(.data$chr,.data$start,
                                                        .data$end,.data$score)
    out.df[1,"start"]=1
    outpath<-gsub('\\.hic$',paste0('_',chrom,'_eigenvectors.txt'),hicfile)
    data.table::fwrite(out.df,outpath,row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
    filepaths<-c(filepaths,path.expand(outpath))
    }
    return(filepaths)
}
    
