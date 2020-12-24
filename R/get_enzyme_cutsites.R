#' get_enzyme_cutsites
#'
#' This function finds all restriction enzyme cutsites of a given genome,
#' genome version, and set of cut patterns
#'@import BSgenome
#'@param gen name of the species: e.g., default \code{'Hsapiens'}
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}
#'@param chrs a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes (except Y and M)
#'in the genome specified by \code{gen} and \code{gen_ver}.
#'@param sig a set of restriction enzyme cut patterns (e.g.,
#''GATC' or c('GATC','GANTC'))
#'@return list of chromosomes.
#'@examples get_enzyme_cutsites(gen='Hsapiens',gen_ver='hg19',
#'sig=c('GATC','GANTC'),chrs=c('chr22'))
#'@export

get_enzyme_cutsites <- function(sig, gen = "Hsapiens", gen_ver = "hg19", chrs = NULL) {
    genome <- paste("BSgenome.", gen, ".UCSC.", gen_ver, sep = "")
    library(genome, character.only = TRUE)
    if (is.null(chrs)) {
        # get list of chromosomes
        chrs <- get_chrs(gen, gen_ver)
        chrs <- chrs[!chrs %in% c("chrY", "chrM")]
    }
    genome.chromSizes <- get_chr_sizes(gen, gen_ver, chrs)
    while (any(grepl("N", sig))) {
        if (length(sig) <= 1) {
            sig <- unique(c(sub("N", "A", sig), sub("N", "C", sig), sub("N", "G", sig), sub("N", "T", sig)))
        } else {
            rng <- length(sig)
            sig_new <- NULL
            for (i in seq(rng)) {
                if (grepl("N", sig[i])) {
                sig_new_i <- unique(c(sub("N", "A", sig[i]), sub("N", "C", sig[i]), sub("N", "G", sig[i]), sub("N", "T", sig[i])))
                sig_new <- c(sig_new, sig_new_i)
                } else {
                sig_new <- unique(c(sig_new, sig[i]))
                }
            }
            sig <- unique(sig_new)
        }
    }
    if (length(sig) <= 1) {
        enzymeCuts <- lapply(chrs, function(x) {
            range_chr <- IRanges::ranges(Biostrings::matchPattern(sig, get(gen)[[x]]))
            if (length(IRanges::start(range_chr)) <= 0) {
                warning(paste0("chromosome ", x, " does not have any cutsites with
                pattern ", sig, " .Using full chromosome instead"))
                range_chr <- IRanges::IRanges(start = 1, end = genome.chromSizes[x])
            }
            sort(GenomicRanges::GRanges(seqnames = x, ranges = range_chr))
        })
    } else {
        patternobj <- Biostrings::DNAStringSet(sig)
        enzymeCuts <- lapply(chrs, function(x) {
            range_chr <- IRanges::ranges(Biostrings::extractAllMatches(get(gen)[[x]], Biostrings::matchPDict(patternobj, get(gen)[[x]])))
            if (length(IRanges::start(range_chr)) <= 0) {
                warning(paste0("chromosome ", x, " does not have any cutsites with
                pattern set ", paste(sig, collapse = ","), " .Using full chromosome instead"))
                range_chr <- IRanges::IRanges(start = 1, end = genome.chromSizes[x])
            }
            sort(GenomicRanges::GRanges(seqnames = x, ranges = range_chr))
        })
    }
    enzymeCuts <- suppressWarnings(sort(
        unlist(methods::as(enzymeCuts, "GRangesList"))))
    return(enzymeCuts)
}
