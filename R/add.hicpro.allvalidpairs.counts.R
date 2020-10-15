#'add.hicpro.allvalidpairs.counts
#'
#'This function converts HiC-Pro outputs in allValidPairs format into a
#'gi_list instance. 
#'
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@importFrom data.table .N
#'@param gi_list valid gi_list instance. 
#'See \code{?gi.list.validate} for details. You can also detect whether
#'a gi_list instance is uniformly binned, along with its bin size
#'using \code{gi.list.binsize.detect}.
#'@param allvalidpairs_path allValidPairsfile obtained from HiC-Pro (e.g.,
#''GSM2572593_con_rep1.allvalidPairs.txt')
#'@param chrs a subset of chromosomes' e.g., c('chr21','chr22'). Defaults
#'to all chromosomes in the \code{gi_list} instance.
#'@param binned TRUE if the gi_list instance is uniformly binned 
#'(helps faster execution). Defaults to TRUE.
#'@param add_inter Interchromosomal interaction counts added as a 1D feature
#'named 'inter' on regions metadata handle of each gi_list element (e.g., 
#'\code{gi_list[[1]]@regions@elementMetadata} or not;
#'default FALSE
#'@return \code{gi_list} instance with counts on the metadata (e.g., 
#'\code{mcols(gi_list[[1]])} handle on each list element, and 'inter' on 
#'regions metadata handle of each element if \code{add_inter=TRUE}.
#'@export

add.hicpro.allvalidpairs.counts <- function(gi_list, allvalidpairs_path,
                                chrs = NULL, binned = TRUE, add_inter = FALSE) {
    gi.list.validate(gi_list)
    if (is.null(chrs)) {
        chrs <- sort(names(gi_list))
    } else {
        chrs <- chrs[chrs %in% names(gi_list)]
        if (length(chrs) == 0) {
            stop("None of the chromosomes specified exists in the gi.list object")
        }
    }
    Dthreshold <- gi.list.Dthreshold.detect(gi_list)
    allvalidpairs <- data.table::fread(allvalidpairs_path, sep = "\t", header = FALSE, select = c(2, 3, 5, 6), stringsAsFactors = FALSE)
    allvalidpairs <- allvalidpairs %>% dplyr::filter(.data$V2 == .data$V5 & .data$V2 %in% chrs & abs(.data$V3 - .data$V6) <= 
        Dthreshold)
    if (binned) {
        binsize <- gi.list.binsize.detect(gi_list)
        allvalidpairs <- allvalidpairs %>% dplyr::rename(chr = "V2") %>% dplyr::mutate(startI = floor(base::pmin(.data$V3, 
            .data$V6)/binsize) * binsize, startJ = floor(base::pmax(.data$V3, .data$V6)/binsize) * binsize) %>% dplyr::select(.data$chr, 
            .data$startI, .data$startJ)
        allvalidpairs <- data.frame(data.table::data.table(allvalidpairs)[, .N, by = names(allvalidpairs)]) %>% dplyr::rename(counts = "N")
    } else {
        allvalidpairs <- allvalidpairs %>% dplyr::rename(chr = "V2") %>% dplyr::mutate(startI = base::pmin(.data$V3, .data$V6), 
            startJ = base::pmax(.data$V3, .data$V6), counts = 1) %>% dplyr::select(.data$chr, .data$startI, .data$startJ, 
            .data$counts)
    }
    for (chrom in chrs) {
        allvalidpairs_chr <- allvalidpairs %>% dplyr::filter(.data$chr == chrom) %>% dplyr::select(.data$startI, .data$startJ, 
            .data$counts)
        gi_list[[chrom]] <- add.2D.features(gi_list[[chrom]], allvalidpairs_chr)
        print(paste0("Chromosome ",chrom," intrachromosomal counts processed."))
    }
    rm(allvalidpairs_chr, allvalidpairs)
    if (add_inter) {
        # inter--allvalidpairs
        allvalidpairs <- data.table::fread(allvalidpairs_path, sep = "\t", header = FALSE, select = c(2, 3, 5, 6), stringsAsFactors = FALSE)
        allvalidpairs <- allvalidpairs %>% dplyr::filter(!.data$V2 == .data$V5 & (.data$V2 %in% chrs | .data$V5 %in% chrs))
        if (binned) {
            binsize <- gi.list.binsize.detect(gi_list)
            allvalidpairs <- allvalidpairs %>% dplyr::rename(chrI = "V2", chrJ = "V5") %>% dplyr::mutate(startI = floor(.data$V3/binsize) * 
                binsize, startJ = floor(.data$V6/binsize) * binsize) %>% dplyr::select(.data$chrI, .data$startI, .data$chrJ, 
                .data$startJ)
            allvalidpairs <- data.frame(data.table::data.table(allvalidpairs)[, .N, by = names(allvalidpairs)]) %>% dplyr::rename(counts = "N")
        } else {
            allvalidpairs <- allvalidpairs %>% dplyr::rename(chrI = "V2", chrJ = "V5") %>% dplyr::mutate(startI = base::pmin(.data$V3, 
                .data$V6), startJ = base::pmax(.data$V3, .data$V6), counts = 1) %>% dplyr::select(.data$chrI, .data$startI, 
                .data$chrJ, .data$startJ, .data$counts)
        }
        for (chrom in chrs) {
            allvalidpairs_chr <- dplyr::bind_rows(allvalidpairs %>% dplyr::filter(.data$chrI == chrom) %>% dplyr::rename(chr = "chrI", 
                start = "startI", inter = "counts") %>% dplyr::select(.data$chr, .data$start, .data$inter), allvalidpairs %>% 
                dplyr::filter(.data$chrJ == chrom) %>% dplyr::rename(chr = "chrJ", start = "startJ", inter = "counts") %>% 
                dplyr::select(.data$chr, .data$start, .data$inter))
            gi_list <- add.1D.features(gi_list, allvalidpairs_chr, chrs = chrom, agg = sum)
            print(paste0("Chromosome ",chrom," interchromosomal counts processed."))
        }
        rm(allvalidpairs_chr, allvalidpairs)
    }
    return(gi_list)
}
