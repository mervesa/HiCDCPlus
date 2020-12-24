#' add_2D_features
#'
#'Adds 2D features to a gi_list instance. If any bin on gi_list overlaps with
#'multiple feature records, features are aggregated among matches according
#'to the univariate vector valued function agg (e.g., sum, mean). For efficient
#'use of memory, using add/expand 1D features (see \code{?add_1D_features} and
#'\code{expand_1D_features}) in sequence is recommended
#'instead of using \code{add_2D_features} directly for each chromosome.
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@importFrom S4Vectors mcols<- mcols queryHits subjectHits
#'@param gi Element of a valid \code{gi_list} instance
#'(restricted to a single chromosome e.g., \code{gi_list[['chr9']]}---see
#'\code{?gi_list_validate} for a detailed explanation of valid \code{gi_list}
#'instances). 
#'@param df data frame for a single chromosome containing columns named chr,
#'startI and startJ and features to be added with their respective names 
#'(if df contains multiple chromosomes, you can convert it into a
#'list of smaller data.frames for each chromosome and apply this function
#'with \code{sapply}).
#'@param features features to be added. Needs to be subset of 
#'\code{colnames(df)}. Defaults to all columns in \code{df} other than
#''chr','start',and 'end'.
#'@param agg any vector valued function with one data argument:
#' defaults to \code{mean}.
#'@return a gi_list element with 2D features stored in metadata handle
#'(i.e., \code{mcols(gi)}).
#'@examples 
#'df<-data.frame(chr='chr9',start=seq(1e6,10e6,1e6))
#'gi_list<-generate_df_gi_list(df,Dthreshold=500e3)
#'feats<-data.frame(chr='chr9',
#'startI=seq(1e6,10e6,1e6),startJ=seq(1e6,10e6,1e6),counts=rpois(10,lambda=5))
#'gi_list[['chr9']]<-add_2D_features(gi_list[['chr9']],feats)
#'@export

add_2D_features <- function(gi, df, features = NULL, agg = sum) {
    if (!methods::is(gi, "GInteractions")) {
        stop("gi_list has to be a list of
            InteractionSet::GInteractions objects")
    }
    if (!("startI" %in% colnames(df) & "startJ" %in% colnames(df))) {
        stop("df should have 'chr','startI',and 'startJ' columns")
    }
    if (is.null(features)) {
        features <- colnames(df)[!colnames(df) %in% c("chr", "startI", "startJ")]
    } else {
        df <- df %>% dplyr::select(c("startI", "startJ", features))
    }
    chrom <- levels(GenomicRanges::seqnames(InteractionSet::regions(gi)))
    if (length(chrom) > 1) {
        stop("Multiple chromosomes in df.")
    }
    dfGI <- GenomicInteractions::GenomicInteractions(GenomicRanges::GRanges(chrom, IRanges::IRanges(df$startI, df$startI + 
        1)), GenomicRanges::GRanges(chrom, IRanges::IRanges(df$startJ, df$startJ + 1)))
    overlaps <- GenomicRanges::findOverlaps(gi, dfGI, minoverlap = 2)
    df$queryHits <- NA
    df$queryHits[subjectHits(overlaps)] <- queryHits(overlaps)
    if (sum(is.na(df$queryHits)) > 0) {
        warning("Bins and counts mismatch. This will slow down the
            performance of counts integration.Check if genome and/or bin size of
            counts data is aligned with the GenomicInteractions object.")
        df <- df %>% dplyr::filter(!is.na(.data$queryHits))
    }
    if (!length(unique(df$queryHits)) == nrow(df)) {
        df <- df %>% dplyr::group_by(.data$queryHits) %>% dplyr::select(.data$queryHits, dplyr::all_of(features)) %>% dplyr::summarize_all(agg)
    }
    for (feature in features) {
        mcols(gi)[, feature] <- 0
        mcols(gi)[, feature][df$queryHits] <- as.numeric(unlist(df[, feature]))
    }
    return(gi)
}

