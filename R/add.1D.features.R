#' add.1D.features.R
#'
#'Adds 1D features to the gi_list instance. If any bin on gi_list
#'overlaps with multiple feature records, feature values are aggregated
#'for the bin according to the
#'vector valued function agg (e.g., sum, mean)
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@importFrom S4Vectors mcols<- mcols queryHits subjectHits
#'@param gi_list List of \code{GenomicInteractions} objects where each object
#' named with chromosomes contains intrachromosomal interaction information
#' (see
#' \code{?gi.list.validate} for a detailed explanation of valid \code{gi_list} instances). 
#'@param df DataFrame with columns named 'chr', and'start' and features to
#'be added with their respective names.
#'@param chrs a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes
#'specified in the data frame \code{df}.
#'@param features features to be added. Needs to be a subset of 
#'\code{colnames(df)}. Defaults to all columns in \code{df} other than
#''chr','start',and 'end'.
#'@param agg any vector valued function with one data argument:
#' defaults to \code{mean}.
#'@return a gi_list instance with 1D features stored in regions metadata handle
#' of each list element (e.g., \code{gi_list[[1]]@regions@elementMetadata})
#' in the instance
#'@examples 
#'df<-data.frame(chr='chr9',start=seq(1e6,10e6,1e6),end=seq(2e6,11e6,1e6))
#'gi_list<-generate.df.gi.list(df)
#'feats<-data.frame(chr='chr9',start=seq(1e6,10e6,1e6),gc=runif(10))
#'gi_list<-add.1D.features(gi_list,feats)
#'@export

add.1D.features <- function(gi_list, df, chrs = NULL, features = NULL, agg = mean) {
    gi.list.validate(gi_list)
    if (!("chr" %in% colnames(df) & "start" %in% colnames(df))) {
        stop("df should have 'chr' and 'start' columns")
    }
    if (is.null(features)) {
        features <- colnames(df)[!colnames(df) %in% c("chr", "start", "end")]
    } else {
        df <- df[c("chr", "start", features)]
    }
    if (is.null(chrs)) 
        chrs <- sort(unique(df$chr))
    for (chrom in chrs) {
        df_chr <- as.data.frame(df %>% dplyr::filter(.data$chr == chrom) %>% dplyr::select(-.data$chr))
        df_chrGR <- GenomicRanges::GRanges(seqnames = chrom, IRanges::IRanges(start = df_chr$start, end = df_chr$start + 1))
        overlaps <- GenomicRanges::findOverlaps(InteractionSet::regions(gi_list[[chrom]]), df_chrGR, minoverlap = 2)
        df_chr$queryHits <- NA
        df_chr$queryHits[subjectHits(overlaps)] <- queryHits(overlaps)
        if (!length(unique(df_chr$queryHits[!is.na(df_chr$queryHits)])) == nrow(df_chr)) {
            df_chr <- df_chr %>% dplyr::filter(!is.na(.data$queryHits)) %>% dplyr::group_by(.data$queryHits) %>% dplyr::summarize_all(agg) %>% 
                dplyr::ungroup()
        }
        for (feature in features) {
            gi_list[[chrom]]@regions@elementMetadata[[feature]] <- 0
            gi_list[[chrom]]@regions@elementMetadata[[feature]][df_chr$queryHits] <- as.numeric(unlist(df_chr[, feature]))
        }
    }
    return(gi_list)
}
