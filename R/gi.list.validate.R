#' gi.list.validate
#'
#' This function validates a gi_list instance.
#'@param gi_list gi_list object to be verified. In order to pass without 
#'errors,  a gi_list object (1) has to be a list of
#'InteractionSet::GInteractions objects, (2) each list element has to be named
#'as chromosomes and only contain intra-chromosomal interaction information,
#'(3) \code{mcols(.)} for each list element should at least contain pairwise
#'genomic distances in a column named 'D'.
#'@return invisible value if the gi_list instance is valid. Otherwise, an error
#'is raised.
#'@examples gi_list<-generate.binned.gi.list(1e6,chrs='chr22')
#'gi.list.validate(gi_list)
#'@export
gi.list.validate <- function(gi_list) {
    if (!is.list(gi_list)) {
        stop("gi_list has to be a list of
            InteractionSet::GInteractions objects")
    }
    if (!all(vapply(gi_list, function(x) methods::is(x, "GInteractions"), TRUE))) {
        stop("gi_list has to be a list of
            InteractionSet::GInteractions objects")
    }
    if (!all(vapply(gi_list, function(x) length(unique(GenomicRanges::seqnames(InteractionSet::regions(x)))), 8) == 1)) {
        stop("Multiple chromosomes in some gi_list list elements")
    }
    if (!all(vapply(gi_list, function(x) "D" %in% colnames(mcols(x)), TRUE))) {
        stop("gi_list list elements should at least contain pairwise genomic distances")
    }
    return(invisible(0))
}

