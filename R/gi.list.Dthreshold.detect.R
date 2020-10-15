#' gi.list.Dthreshold.detect
#'
#' This function finds the maximum genomic distance in a valid gi_list object.
#'@param gi_list A valid gi_list instance. See
#' \code{?gi.list.validate} for more details about
#' the attributes of a valid gi_list instance. 
#'@return maximum genomic distance in the object
#'@examples gi_list<-generate.binned.gi.list(1e6,chrs='chr22')
#'gi.list.Dthreshold.detect(gi_list)
#'@export
gi.list.Dthreshold.detect <- function(gi_list) {
    gi.list.validate(gi_list)
    get.Dthreshold <- function(gi) {
        return(max(mcols(gi)$D, na.rm = TRUE))
    }
    range_candidates <- vapply(gi_list, get.Dthreshold, 1)
    Dthreshold <- round(max(range_candidates))
    return(Dthreshold)
}

