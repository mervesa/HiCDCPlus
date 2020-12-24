#' expand_1D_features
#'
#'Expands 1D features on the regions metadata handle
#'of each list element (e.g., \code{gi_list[[1]]@regions@elementMetadata})
#'to the to 2D metadata e.g., \code{mcols(gi_list[[1]])}). Two feature values
#'corresponding to each anchor is summarized as a score using a vector
#'valued function agg that takes two vector valued arguments of the same size
#'and outputs a vector of the same size as the input vectors. This defaults
#'to the \code{transform.vec} function outlined in (Carty et al., 2017). 
#'For efficient use of memory, using add/expand 1D features (see 
#'\code{?add_1D_features} and \code{expand_1D_features}) in sequence is
#'recommended instead of using \code{add_2D_features} directly
#'for each chromosome.
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param gi_list List of \code{GenomicInteractions} objects where each object
#' named with chromosomes contains intra-chromosomal interaction information
#' (see
#' \code{?gi_list_validate} for a detailed explanation of valid \code{gi_list}
#'  instances). 
#'@param chrs a subset of chromosomes' e.g., c('chr21','chr22'). Defaults
#'to all chromosomes in the \code{gi_list} instance.
#'@param features features to be added. Defaults to all 1D features in
#' elements of \code{gi_list[[1]]@regions@elementMetadata}
#'@param agg any vector valued function with two data arguments:
#' defaults to \code{transform.vec} described in HiC-DC (Carty et al., 2017).
#'@return a gi_list element with 2D features stored in metadata handle
#'(i.e., \code{mcols(gi)}).
#'@examples 
#'df<-data.frame(chr='chr9',start=seq(1e6,10e6,1e6),end=seq(2e6,11e6,1e6))
#'gi_list<-generate_df_gi_list(df)
#'feats<-data.frame(chr='chr9',start=seq(1e6,10e6,1e6),gc=runif(10))
#'gi_list<-add_1D_features(gi_list,feats)
#'gi_list<-expand_1D_features(gi_list)
#'@export

expand_1D_features <- function(gi_list, chrs = NULL, features = NULL, agg = transform.vec) {
    # gi_list: list of GenomicInteractions objects per chromosome chrs: chromosomes to be expanded (default: all chromosomes
    # in the list) features: features to be expanded (default: all 1D features) agg: any vector valued function with two data
    # arguments default function transform.vec from HiCDC (Carty et al., 2017) For efficient use of memory, add_1D_features 
    # followed by expand_1D_features where applicable is suggested instead of add_2D_features for each chromosome
    gi_list_validate(gi_list)
    Normalize <- function(x) {
        x.new <- ifelse(is.finite(x), (x - base::mean(x[is.finite(x)]))/stats::sd(x[is.finite(x)]), NA)
        return(x.new)
    }
    transform.vec <- function(x, y) {
        as.vector(ifelse(x == 0 | y == 0, 0, Normalize(log(x * y))))
    }
    if (is.null(chrs)) {
        chrs <- names(gi_list)
    }
    if (is.null(features)) {
        features <- colnames(gi_list[[1]]@regions@elementMetadata)
    }
    for (chrom in chrs) {
        for (feature in features) {
            mcols(gi_list[[chrom]])[, feature] <- agg(gi_list[[chrom]]@regions@elementMetadata[[feature]][InteractionSet::anchorIds(gi_list[[chrom]], 
                type = "first")], gi_list[[chrom]]@regions@elementMetadata[[feature]][InteractionSet::anchorIds(gi_list[[chrom]], 
                type = "second")])
        }
    }
    return(gi_list)
}
