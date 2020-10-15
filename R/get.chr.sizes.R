#' get.chr.sizes
#'
#' This function finds all chromosome sizes of a given genome, genome version
#' and set of chromosomes.
#'@import BSgenome
#'@param gen name of the species: e.g., default \code{'Hsapiens'}
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes (except Y and M)
#'in the genome specified.
#'@return named vector containing names as chromosomes and values as chromosome
#'sizes.
#'@examples get.chr.sizes('Hsapiens','hg19',c('chr21','chr22'))
#'@export
get.chr.sizes <- function(gen = "Hsapiens", gen_ver = "hg19", chrs = NULL) {
    genome <- paste("BSgenome.", gen, ".UCSC.", gen_ver, sep = "")
    library(genome, character.only = TRUE)
    if (is.null(chrs)) {
        # get list of chromosomes
        chrs <- get.chrs(gen, gen_ver)
    }
    sizes <- GenomeInfoDb::seqlengths(get(gen))[chrs]
    return(sizes)
}
