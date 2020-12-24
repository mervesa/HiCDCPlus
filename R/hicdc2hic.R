#' hicdc2hic
#'
#' This function converts various modes from HiCDCPlus gi_list (uniformly binned)
#' instance back into a .hic file with the \code{mode} passed
#' as counts that can be retrieved using Juicer Dump
#' (https://github.com/aidenlab/juicer/wiki/Data-Extraction)
#' with 'NONE' normalization. 
#'@import BSgenome splines
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param gi_list List of \code{GenomicInteractions} objects where each object
#'named with chromosomes contains intra-chromosomal interaction information
#'(minimally containing counts and genomic distance in \code{mcols(gi_list)}---
#'see \code{?gi_list_validate} for a detailed explanation of valid
#'\code{gi_list} instances).
#'@param hicfile the path to the .hic file
#'@param mode What to put to the .hic file
#'as score. Allowable options are: 'pvalue' for -log10 significance p-value,
#''qvalue' for -log10 FDR corrected p-value,
#''normcounts' for raw counts/expected counts, and
#''zvalue' for standardized counts (raw counts-expected counts)/modeled
#'standard deviation of expected counts
#' and 'raw' to pass-through 'raw counts. Defaults to 'normcounts'.
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to chromosomes in \code{gi_list}.
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}
#'@return path of the .hic file.
#'@examples 
#'outdir<-paste0(tempdir(check=TRUE),'/')
#'gi_list<-generate_binned_gi_list(50e3,chrs='chr22')
#'gi_list<-add_hic_counts(gi_list,
#'hic_path=system.file("extdata", "GSE63525_HMEC_combined_example.hic",
#' package = "HiCDCPlus"))
#'hicdc2hic(gi_list,hicfile=paste0(outdir,'out.hic'),
#'mode='raw')
#'@export

hicdc2hic <- function(gi_list, hicfile, mode = "normcounts", chrs = NULL, gen_ver = "hg19") {
    options(scipen = 9999)
    gi_list_validate(gi_list)
    binsize<-gi_list_binsize_detect(gi_list)
    if (is.null(chrs)) 
        chrs <- names(gi_list)
    tmpfile <- paste0(base::tempfile(), ".txt")
    gi_list_write(gi_list, tmpfile, columns = "minimal_plus_score", score = mode)
    #generate path to the file if not exists
    hicdc2hicoutput <- path.expand(hicfile)
    hicdc2hicoutputdir<-gsub("/[^/]+$", "",hicdc2hicoutput)
    if (hicdc2hicoutputdir==hicdc2hicoutput){
        hicdc2hicoutputdir<-gsub("\\[^\\]+$", "",hicdc2hicoutput)
    }
    if (hicdc2hicoutputdir==hicdc2hicoutput){
        hicdc2hicoutputdir<-gsub("\\\\[^\\\\]+$", "",hicdc2hicoutput)
    }
    if (!hicdc2hicoutputdir==hicdc2hicoutput&!dir.exists(hicdc2hicoutputdir)){
        dir.create(hicdc2hicoutputdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
    }
    # run pre
    outdir<-tempdir(check=TRUE)
    jarpath<-paste0(outdir,'/juicer_tools.jar')
    if(!file.exists(jarpath)) {
        if(.Platform$OS.type=="windows"){
            utils::download.file(url='https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar',
                                 destfile=jarpath,quiet=TRUE, mode="wb", method="wininet")
        } else {
            utils::download.file(url='https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar',
                                 destfile=jarpath,quiet=TRUE)
            
        }
    }
    system2("java", args = c("-Xmx8g", "-jar", path.expand(jarpath), "pre", "-v", "-d", "-r", binsize, path.expand(tmpfile), path.expand(hicdc2hicoutput), 
        gen_ver))
    # remove file
    system2("rm", args = path.expand(tmpfile))
    return(hicdc2hicoutput)
}
