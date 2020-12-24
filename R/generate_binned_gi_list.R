#' generate_binned_gi_list
#'
#'Generates a valid uniformly binned gi_list instance.
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param binsize Desired binsize in bp, e.g., 5000, 25000.
#'@param chrs a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes except "Y", and "M" for
#'the specified \code{gen} and \code{gen_ver}.
#'@param Dthreshold maximum distance (included) to check for significant
#'interactions, defaults to 2e6 or maximum in the data; whichever is smaller.
#'@param gen name of the species: e.g., default \code{'Hsapiens'}.
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}.
#'@return a valid, uniformly binned gi_list instance. 
#'@examples gi_list<-generate_binned_gi_list(1e6,chrs='chr22')
#'@export

generate_binned_gi_list <- function(binsize, chrs = NULL, Dthreshold = 2e+06, gen = "Hsapiens", gen_ver = "hg19") {
    # generates binned GenomicInteraction object separated by chromosomes defined in chrs. binsize: binsize in BP chrs:
    # chromosomes for this object to be generated Dthreshold: maximum pairwise genomic distance (default 2Mb)
    gi_list <- list()
    if (is.null(chrs)) {
        chrs <- get_chrs(gen, gen_ver)
    } else {
        chrs <- chrs[chrs %in% get_chrs(gen, gen_ver)]
    }
    for (chrom in chrs) {
        seqlen <- get_chr_sizes(gen, gen_ver)[chrom]
        numbins <- ceiling(seqlen/binsize)
        maxbins <- min(round(Dthreshold/binsize), numbins)
        all.regions <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = seq(0, (numbins - 1), 1) * binsize, end = pmin(seq(1, 
            numbins, 1) * binsize, seqlen)))
        index1 <- unlist(lapply(seq(1, numbins, 1), function(x) rep(x, min(maxbins + 1, numbins - x + 1))))
        index2 <- unlist(lapply(seq(1, numbins, 1), function(x) seq(x, min(x + maxbins, numbins), 1)))
        gi_list[[chrom]] <- InteractionSet::GInteractions(index1, index2, all.regions)
        mcols(gi_list[[chrom]])$D <- InteractionSet::pairdist(gi_list[[chrom]])
        gi_list[[chrom]] <- gi_list[[chrom]][mcols(gi_list[[chrom]])$D <= Dthreshold]
    }
    return(gi_list)
}
