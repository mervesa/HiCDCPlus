#' straw_dump
#'
#' Interface for Juicer's dump in case C++ straw fails (known to fail on Windows
#' due to zlib compression not being OS agnostic and particularly not preserving
#' null bytes, which .hic files are delimited with). This function reads the 
#' .hic file, finds the appropriate matrix and slice of data, writes it to a 
#' temp file, reads and modifies it, and outputs as an R DataFrame (and also
#' deletes the temp file).
#' 
#' Usage: straw_dump <oe/observed> <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize> <outfile>
#' 
#' @param norm Normalization to apply. Must be one of NONE/VC/VC_SQRT/KR.
#'     VC is vanilla coverage, VC_SQRT is square root of vanilla coverage,
#'     and KR is Knight-Ruiz or Balanced normalization.
#' @param fn path to the .hic file
#' @param ch1 first chromosome location (e.g., "1")
#' @param ch2 second chromosome location (e.g., "8")
#' @param u BP (BasePair) or FRAG (restriction enzyme FRAGment)
#' @param bs The bin size. By default, for BP, this is one of 
#'     <2500000, 1000000, 500000,
#'     250000, 100000, 50000, 25000, 10000, 5000> and for FRAG this is one of 
#'     <500, 200,
#'     100, 50, 20, 5, 2, 1>.
#' @return Data.frame of a sparse matrix of data from hic file. x,y,counts
#' @export
straw_dump<-function(norm,fn,ch1,ch2,u,bs){
    options(scipen=9999)
    outdir<-tempdir(check=TRUE)
    jarpath<-paste0(outdir,'/juicer_tools.jar')
    if(!file.exists(jarpath)) {
    if(.Platform$OS.type=="windows"){
    utils::download.file(url='https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar',
        destfile=jarpath,quiet=TRUE, mode="wb", method=as.character(ifelse(capabilities("libcurl"),"libcurl","wininet")))
    } else {
    utils::download.file(url='https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar',
        destfile=jarpath,quiet=TRUE)
        
    }
    }
    tmpfile <- path.expand(base::tempfile())
    chrflag1<-grepl("^chr",ch1)
    chrflag2<-grepl("^chr",ch1)
    tryCatch({
        system2("java", args = c("-Xmx8g", "-jar", path.expand(jarpath), "dump",
                                 "observed", 
                                 norm,
                                 path.expand(fn),
                                 ifelse(chrflag1,gsub("chr","",ch1),ch1),
                                 ifelse(chrflag2,gsub("chr","",ch2),ch2),
                                 u,
                                 bs,
                                 tmpfile))
        df<-as.data.frame(data.table::fread(tmpfile))},
        error=function(e){
            system2("java", args = c("-Xmx8g", "-jar", path.expand(jarpath), "dump",
                                     "observed", 
                                     "NONE",
                                     path.expand(fn),
                                     ifelse(chrflag1,ch1,paste0("chr",ch1)),
                                     ifelse(chrflag2,ch2,paste0("chr",ch2)),
                                     u,
                                     bs,
                                     tmpfile))
            df<-as.data.frame(data.table::fread(tmpfile))
        })
    colnames(df)<-c('x','y','counts')
    system2("rm", args = tmpfile)
    return(df)
}