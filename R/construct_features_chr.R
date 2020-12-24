#' construct_features_chr
#'
#'This function lists all restriction enzyme cutsites of a given genome and
#'genome version with genomic features outlined in Carty et al. (2017) for
#'a single chromosome.
#'https://www.nature.com/articles/ncomms15454; GC content, mappability,
#'and effective length
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param chrom select a chromosome.
#'@param gen name of the species: e.g., default \code{'Hsapiens'}.
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}.
#'@param sig restriction enzyme cut pattern (or a vector of patterns; e.g.,
#''GATC' or c('GATC','GANTC')).
#'@param bin_type 'Bins-uniform' if uniformly binned by binsize in
#'bp, or 'Bins-RE-sites' if binned by number of
#'restriction enzyme fragments.
#'@param binsize binsize in bp if bin_type='Bins-uniform' (or number of
#'RE fragment cut sites if bin_type='Bins-RE-sites'), defaults to 5000.
#'@param wg_file path to the bigWig file containing mappability values across
#'the genome of interest.
#'@return a features 'bintolen' file that contains GC, mappability and length
#'features.
#'@examples df<-construct_features_chr(chrom='chr22',
#'gen='Hsapiens', gen_ver='hg19',sig=c('GATC','GANTC'),bin_type='Bins-uniform',
#'binsize=100000,wg_file=NULL)
#'@export

construct_features_chr <- function(chrom, gen = "Hsapiens", gen_ver = "hg19",
                               sig = "GATC", bin_type = "Bins-uniform", binsize = 5000, 
                               wg_file = NULL) {
    
    getGC <- function(regions) {
        genome <- paste("BSgenome.", gen, ".UCSC.", gen_ver, sep = "")
        library(genome, character.only = TRUE)
        gc <- rowSums(Biostrings::alphabetFrequency(BSgenome::Views(BiocGenerics::get(genome), regions), as.prob = TRUE)[, 
                                                                                                                         seq(2,3,1)])
        return(gc)
    }
    
    getMap <- function(regions, data) {
        hits <- GenomicRanges::findOverlaps(data, regions, type = "within")
        DT <- data.frame(queryHits = as.data.frame(hits)[[1]], subjectHits = as.data.frame(hits)[[2]])
        DT <- DT %>% dplyr::mutate(map = data$score[.data$queryHits], len = GenomicRanges::end(data)[.data$queryHits] - GenomicRanges::start(data)[.data$queryHits] + 
                                       1) %>% dplyr::group_by(.data$subjectHits) %>% dplyr::summarize(avgmap = stats::weighted.mean(map, w = len))
        map <- rep(0, length(regions))
        map[DT$subjectHits] <- DT$avgmap
        return(map)
    }
    
    # get chromosome size
    genome.chromSize <- get_chr_sizes(gen, gen_ver, chrom)[chrom]
    genome.chromGR <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = 1, end = genome.chromSize))
    
    print(paste0("Using ", paste(chrom, collapse = " "), "and cut patterns ", paste(sig, collapse = " ")))
    
    
    enzymeCutsites <- get_enzyme_cutsites(sig, gen, gen_ver, chrom)
    RE_sites <- enzymeCutsites
    newends <- dplyr::lead(BiocGenerics::start(RE_sites)) - 1
    newends <- c(newends[-length(newends)], genome.chromSize)
    BiocGenerics::end(RE_sites) <- newends
    minstart <- min(BiocGenerics::start(RE_sites))
    if (minstart > 1) {
        RE_sites <- GenomeInfoDb::sortSeqlevels(c(GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = 1, 
                                                                                                                       end = minstart)), RE_sites))
        RE_sites <- sort(RE_sites)
    }
    
    
    if (!(is.null(wg_file))) {
        wgdata <- rtracklayer::import.bw(wg_file, which = genome.chromGR, as = "GRanges")
    } else {
        wgdata <- NULL
    }
    
    if (bin_type == "Bins-RE-sites") {
        endL <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(start = GenomicRanges::start(RE_sites), 
                                                                                                                    width = 200))
        BiocGenerics::end(endL) <- pmin(BiocGenerics::end(endL), genome.chromSize)
        BiocGenerics::start(endL) <- pmax(BiocGenerics::start(endL), 1)

        endR <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(end = GenomicRanges::end(RE_sites), 
                                                                                                                    width = 200))

        BiocGenerics::end(endR) <- pmin(BiocGenerics::end(endR), genome.chromSize)
        BiocGenerics::start(endR) <- pmax(BiocGenerics::start(endR), 1)
        gcL <- getGC(regions = endL)
        gcR <- getGC(regions = endR)
        gc <- (gcL + gcR)/2
        RE_sites$gc <- gc
        
        if (!is.null(wgdata)) {
            endL <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(start = GenomicRanges::start(RE_sites), 
                                                                                                                        width = 500))
                BiocGenerics::end(endL) <- pmin(BiocGenerics::end(endL), genome.chromSize)
                BiocGenerics::start(endL) <- pmax(BiocGenerics::start(endL), 1)
            endR <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(end = GenomicRanges::end(RE_sites), 
                                                                                                                        width = 500))
                BiocGenerics::end(endR) <- pmin(BiocGenerics::end(endR), genome.chromSize)
                BiocGenerics::start(endR) <- pmax(BiocGenerics::start(endR), 1)
            mapL <- getMap(regions = endL, data = wgdata)
            mapR <- getMap(regions = endR, data = wgdata)
            map <- (mapL + mapR)/2
            RE_sites$map <- map
        } else {
            RE_sites$map <- 0
        }
        # generate bintolen from RE_sites
        bintolen <- as.data.frame(RE_sites, stringsAsFactors = FALSE)
        bintolen <- bintolen %>% dplyr::mutate(RE_id = dplyr::row_number())
        if (binsize > 1) {
            # We cluster enzyme cutting sites every 'binsize' consecutive RE sites per each chromosome
            bintolen <- bintolen %>% dplyr::mutate(binNumb = floor((dplyr::row_number() - 1)/binsize) + 
                                                                                     1)
            bintolen <- bintolen %>% dplyr::group_by(.data$seqnames, .data$binNumb) %>% dplyr::summarize(start = min(start), 
                                                                                                         end = max(end), gc = mean(gc), map = mean(map), RE_id = min(.data$RE_id)) %>% dplyr::mutate(width = .data$end - 
                                                                                                                                                                                                         .data$start + 1) %>% dplyr::rename(chr = "seqnames") %>% dplyr::mutate(bins = paste(.data$chr, .data$start, 
                                                                                                                                                                                                                                                                                             .data$end, sep = "-"))
        } else {
            bintolen <- bintolen %>% dplyr::rename(chr = "seqnames") %>% dplyr::mutate(bins = paste(.data$chr, .data$start, 
                                                                                                    .data$end, sep = "-"))
        }
        # construct a similar bintolen object
        bintolen <- bintolen %>% dplyr::select(.data$bins, .data$gc, .data$map, .data$width, .data$RE_id)
        bintolen <- as.data.frame(bintolen, stringsAsFactors = FALSE)
        rownames(bintolen) <- bintolen$bins
        return(bintolen)
    }
    
    if (bin_type == "Bins-uniform") {
        bins.chrom <- function(chrom, binsize) {
            cuts_start <- seq(1, genome.chromSize, by = binsize)
            cuts_end <- seq(binsize, genome.chromSize - (genome.chromSize%%binsize) + binsize, by = binsize)
            cuts_end <- c(cuts_end[-length(cuts_end)], genome.chromSize)
            bins <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = cuts_start, end = cuts_end))
            return(bins)
        }
        binsGR <- bins.chrom(chrom, binsize)
        names(binsGR) <- paste(as.character(GenomicRanges::seqnames(binsGR)), GenomicRanges::start(binsGR), GenomicRanges::end(binsGR), 
                               sep = "-")
        

        medians <- (GenomicRanges::start(enzymeCutsites) + GenomicRanges::end(enzymeCutsites))/2
        if (!is.null(wgdata)) {
                wgdata <- wgdata[GenomicRanges::seqnames(wgdata) == chrom]
        } else {
                wgdata <- NULL
        }
        # 500bp from ends for mappability and len features
        FragmentendsL <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::restrict(IRanges::IRanges(end = medians - 
                                                                                                                               1, width = 500), start = 1, end = genome.chromSize))
        FragmentendsR <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::restrict(IRanges::IRanges(start = medians, 
                                                                                                                           width = 500), start = 1, end = genome.chromSize))
        ends <- c(FragmentendsL, FragmentendsR)
        hits <- as.data.frame(GenomicRanges::findOverlaps(ends, binsGR, type = "within", select = "all"))
        LR <- data.frame(bins = names(binsGR[hits[[2]]]), start = GenomicRanges::start(ends[hits[[1]]]), end = GenomicRanges::end(ends[hits[[1]]]), 
                        stringsAsFactors = FALSE)
        LRgr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = LR$start, end = LR$end))
        if (!is.null(wgdata)) {
            LR$map <- getMap(LRgr, wgdata)
        } else {
            LR$map <- 0
        }
        map <- LR %>% dplyr::group_by(.data$bins) %>% dplyr::summarize(map = mean(map))
        # effective length
        ir <- IRanges::IRanges(LR$start, LR$end)
        LR$group <- as.data.frame(GenomicRanges::findOverlaps(ir, IRanges::reduce(ir)))[[2]]
        len <- LR %>% dplyr::group_by(.data$group, .data$bins) %>% dplyr::summarize(start = min(start), end = max(end)) %>% 
            dplyr::mutate(width = .data$end - .data$start + 1) %>% dplyr::group_by(.data$bins) %>% dplyr::summarize(len = sum(.data$width))
            
        # 200bp from ends for gc feature
        FragmentendsL <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(end = medians - 1, 
                                                                                                        width = 200))
        FragmentendsR <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = medians, 
                                                                                                         width = 200))
        ends <- c(FragmentendsL, FragmentendsR)
        hits <- as.data.frame(GenomicRanges::findOverlaps(ends, binsGR, type = "within", select = "all"))
        LR2 <- data.frame(bins = names(binsGR[hits[[2]]]), start = GenomicRanges::start(ends[hits[[1]]]), end = GenomicRanges::end(ends[hits[[1]]]), 
                              stringsAsFactors = FALSE)
        LR2gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = LR2$start, end = LR2$end))
        LR2$gc <- getGC(LR2gr)
        gc <- LR2 %>% dplyr::group_by(.data$bins) %>% dplyr::summarize(gc = mean(gc))
        LR <- dplyr::left_join(dplyr::left_join(gc, map), len)
        bintolen <- as.data.frame(LR)
        rownames(bintolen) <- bintolen$bins
        allbins <- data.frame(bins = names(binsGR), stringsAsFactors = FALSE)
        bintolen <- suppressWarnings(dplyr::left_join(allbins, bintolen) %>% 
                            tidyr::replace_na(list(gc = 0, map = 0, len = 0)))
        return(bintolen)


    }
}
    
