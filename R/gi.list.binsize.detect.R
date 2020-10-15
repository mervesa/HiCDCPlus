#' gi.list.binsize.detect
#'
#' This function finds the bin size of a uniformly binned valid gi_list instance
#' in bp. It raises an error if the gi_list instance is not uniformly binned.
#'@param gi_list gi_list object to be verified. In order to pass without errors, 
#'a gi_list object (1) has to be a list of InteractionSet::GInteractions
#'objects,(2) each list element has to be named as chromosomes and only contain
#'intra-chromosomal interaction information,
#'(3) \code{mcols(.)} for each list element should at least contain pairwise
#'genomic distances in a column named 'D' and (4) each list element needs to
#'be uniformly binned
#'@return uniform binsize in base pairs or an error if the gi_list instance
#'is not uniformly binned.
#'@examples gi_list<-generate.binned.gi.list(1e6,chrs='chr22')
#'gi.list.binsize.detect(gi_list)
#'@export
gi.list.binsize.detect <- function(gi_list) {
    gi.list.validate(gi_list)
    get.range <- function(gi) {
        return(stats::quantile(GenomicRanges::end(InteractionSet::regions(gi)) - GenomicRanges::start(InteractionSet::regions(gi)), 
                               probs = c(0.025, 0.975)))
    }
    if('chrM'%in%names(gi_list)){
        sel<-names(gi_list)[!names(gi_list)%in%'chrM']
        range_candidates<-unique(c(vapply(gi_list[sel], get.range, c(1, 2))))
    }else{
        range_candidates <- unique(c(vapply(gi_list, get.range, c(1, 2))))
    }
    if (length(range_candidates) == 1 | stats::sd(range_candidates) < 1e-04) {
        binsize <- round(base::mean(range_candidates))
    } else {
        stop("gi_list has variable binsizes")
    }
    return(binsize)
}