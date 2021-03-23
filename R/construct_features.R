#' construct_features
#'
#'This function lists all restriction enzyme cutsites of a given genome and
#'genome version with genomic features outlined in Carty et al. (2017)
#'https://www.nature.com/articles/ncomms15454; GC content, mappability,
#'and effective length
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param output_path the path to the folder and name prefix you want to place
#'feature files into. The feature file will have the suffix '_bintolen.txt.gz'.
#'@param gen name of the species: e.g., default \code{'Hsapiens'}.
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}.
#'@param sig restriction enzyme cut pattern (or a vector of patterns; e.g.,
#''GATC' or c('GATC','GANTC')).
#'@param bin_type 'Bins-uniform' if uniformly binned by binsize in
#'bp, or 'Bins-RE-sites' if binned by number of
#'restriction enzyme fragments.
#'@param binsize binsize in bp if bin_type='Bins-uniform' (or number of
#'RE fragment cut sites if bin_type='Bins-RE-sites'), defaults to 5000.
#'@param wg_file path to the bigWig file containing mappability values across
#'the genome of interest.
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes (except Y and M)
#'in the genome specified.
#'@param feature_type 'RE-based' if features are to be computed based on
#'restriction enzyme fragments. 'RE-agnostic' ignores restriction enzyme cutsite 
#'information and computes features gc and map based on binwide averages. bin_type
#'has to be 'Bins-uniform' if \code{feature_type='RE-agnostic'}.
#'@return a features 'bintolen' file that contains GC, mappability and length
#'features.
#'@examples 
#'outdir<-paste0(tempdir(check=TRUE),'/')
#'construct_features(output_path=outdir,gen='Hsapiens',
#'gen_ver='hg19',sig=c('GATC','GANTC'),bin_type='Bins-uniform',binsize=100000,
#'wg_file=NULL,chrs=c('chr21'))
#'@export

construct_features <- function(output_path, gen = "Hsapiens", gen_ver = "hg19",
                      sig = "GATC", bin_type = "Bins-uniform", binsize = 5000, 
                      wg_file = NULL, chrs = NULL, feature_type = "RE-based") {
    
    if (is.null(chrs)) {
        chrs <- get_chrs(gen, gen_ver)
    }
    bintolen<-lapply(chrs,function(x) data.frame(stringsAsFactors = FALSE))
    names(bintolen)<-chrs
    for (chrom in names(bintolen)){
    bintolen[[chrom]] <- construct_features_chr(chrom=chrom,
                                          gen = gen, 
                                          gen_ver = gen_ver,
                                          sig = sig, 
                                          bin_type = bin_type, 
                                          binsize = binsize, 
                                          wg_file = wg_file,
                                          feature_type = feature_type)
    }
    bintolen<-suppressWarnings(dplyr::bind_rows(bintolen))
    bintolenoutput <- path.expand(paste0(output_path, "_bintolen.txt.gz"))
    bintolenoutputdir<-gsub("/[^/]+$", "",bintolenoutput)
    if (bintolenoutputdir==bintolenoutput){
        bintolenoutputdir<-gsub("\\[^\\]+$", "",bintolenoutput)
    }
    if (bintolenoutputdir==bintolenoutput){
        bintolenoutputdir<-gsub("\\\\[^\\\\]+$", "",bintolenoutput)
    }
    if (!bintolenoutputdir==bintolenoutput&!dir.exists(bintolenoutputdir)){
        dir.create(bintolenoutputdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
    }
    data.table::fwrite(bintolen, bintolenoutput, row.names = FALSE, quote = FALSE, sep = "\t")
    return(bintolenoutput)
}
