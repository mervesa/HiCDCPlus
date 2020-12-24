#' get_chrs
#'
#' This function finds all chromosomes of a given genome and genome version
#' except for Y and M.
#'@import BSgenome
#'@param gen name of the species: e.g., default \code{'Hsapiens'}
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}
#'@return string vector of chromosomes.
#'@examples get_chrs('Hsapiens','hg19')
#'@export
get_chrs <- function(gen = "Hsapiens", gen_ver = "hg19") {
    genome <- paste("BSgenome.", gen, ".UCSC.", gen_ver, sep = "")
    library(genome, character.only = TRUE)
    chrs <- GenomeInfoDb::seqnames(get(gen))[vapply(GenomeInfoDb::seqnames(get(gen)), nchar, FUN.VALUE = 1) <= 5]
    chrs <- chrs[!chrs %in% c("chrMT","chrY","chrM")]
    return(chrs)
}
