#' generate_bintolen_gi_list
#'
#'Generates a gi_list instance from a bintolen file generated by
#' \code{generate.features} (see \code{?generate.features}) for details).
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param bintolen_path path to the flat file containing columns named bins and
#'features
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes
#'specified in the bintolen file.
#'@param Dthreshold maximum distance (included) to check for significant
#'interactions, defaults to 2e6 or maximum in the data; whichever is smaller.
#'@param binned TRUE if the bintolen file is uniformly binned. Defaults to TRUE.
#'@param binsize bin size in bp to be generated for the object. Defaults to the
#'binsize in the bintolen file, if exists.
#'@param gen name of the species: e.g., default \code{'Hsapiens'}
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}
#'@return a valid gi_list instance with genomic features derived from specified
#'restriction enzyme cut patterns when generating the bintolen file using
#'\code{construct_features} (see \code{?construct_features} for help).  
#'Genomic 1D features are stored in the regions metadata handle
#'of each list element (e.g., \code{gi_list[[1]]@regions@elementMetadata}).
#'@examples 
#'chrs<-'chr22'
#'bintolen_path<-system.file("extdata", "test_bintolen.txt.gz",
#' package = "HiCDCPlus")
#'gi_list<-generate_bintolen_gi_list(bintolen_path,chrs)
#'@export

generate_bintolen_gi_list <- function(bintolen_path, chrs = NULL, Dthreshold = 2e+06, binned = TRUE, binsize = NULL, gen = "Hsapiens", 
    gen_ver = "hg19") {
    input.file.read <- function(filepath) {
        # reads files of RDS, .txt, or .txt.gz filepath: valid path ending in .txt, .txt.gz, or .rds
        if (grepl("\\.txt.gz$", filepath) | grepl("\\.txt$", filepath)) {
            return(data.table::fread(filepath))
        } else if (grepl("\\.rds$", filepath)) {
            return(readRDS(filepath))
        } else {
            stop("Can only read in paths ending with .txt,.txt.gz, or .rds")
        }
    }
    if (binned & is.null(binsize) & !is.null(bintolen_path)) {
        bintolen <- input.file.read(bintolen_path)
        # get starts to the nearest 1000 if binned
        bintolen <- bintolen %>% tidyr::separate(.data$bins, c("chr", "start", "end"), "-") %>% dplyr::mutate(start = floor(as.numeric(.data$start)/1000) * 
            1000, end = as.numeric(.data$end))
        binsize <- round(mean(abs(bintolen$end - bintolen$start))/1000) * 1000
    } else if (binned & is.null(binsize) & is.null(bintolen_path)) {
        stop("If binned, need to specify at least one of binsize and bintolen_path")
    }
    # read bintolen if not read already
    if (!is.null(bintolen_path) & !exists("bintolen")) {
        bintolen <- input.file.read(bintolen_path)
        if (binned) {
            # if binned, get starts to the nearest 1000
            bintolen <- bintolen %>% tidyr::separate(.data$bins, c("chr", "start", "end"), "-") %>% dplyr::mutate(start = floor(as.numeric(.data$start)/1000) * 
                1000, end = as.numeric(.data$end))
        } else {
            bintolen <- bintolen %>% tidyr::separate(.data$bins, c("chr", "start", "end"), "-") %>% dplyr::mutate(start = as.numeric(.data$start), 
                end = as.numeric(.data$end))
        }
    }
    if (is.null(chrs)) {
        chrs <- sort(unique(bintolen$chr))
    }
    # generate gi_list
    if (binned) {
        gi_list <- generate_binned_gi_list(binsize, chrs, Dthreshold, gen, gen_ver)
    } else {
        gi_list <- generate_df_gi_list(bintolen, chrs, Dthreshold, gen, gen_ver)
    }
    # merge gc, map, len features
    gi_list <- add_1D_features(gi_list, bintolen, chrs)
    return(gi_list)
}
