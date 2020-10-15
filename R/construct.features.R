#' construct.features.R
#'
#'This function lists all restriction enzyme cutsites of a given genome and
#'genome version with genomic features outlined in Carty et al. (2017)
#'https://www.nature.com/articles/ncomms15454; GC content, mappability,
#'and effective length
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param output_path the path to the folder and name prefix you want to place
#'feature files into. The feature file will have the suffix '_bintolen.txt.gz'.
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
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes (except Y and M)
#'in the genome specified.
#'@return a features 'bintolen' file that contains GC, mappability and length
#'features.
#'@examples 
#'outdir<-paste0(tempdir(check=TRUE),'/')
#'construct.features(output_path=outdir,gen='Hsapiens',
#'gen_ver='hg19',sig=c('GATC','GANTC'),bin_type='Bins-uniform',binsize=100000,
#'wg_file=NULL,chrs=c('chr21'))
#'@export

construct.features <- function(output_path, gen = "Hsapiens", gen_ver = "hg19",
                      sig = "GATC", bin_type = "Bins-uniform", binsize = 5000, 
                      wg_file = NULL, chrs = NULL) {
    
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
    
    # Check if chrom is given as an argument.  If it is entered, converted it to vector:
    if (is.null(chrs)) {
        chrs <- get.chrs(gen, gen_ver)
    }
    # ger chromosome sizes
    genome.chromSizes <- get.chr.sizes(gen, gen_ver, chrs)
    genome.chromGR <- GenomicRanges::GRanges(seqnames = chrs, ranges = IRanges::IRanges(start = 1, end = genome.chromSizes[chrs]))
    
    print(paste0("Using ", paste(chrs, collapse = " "), "and cut patterns ", paste(sig, collapse = " ")))
    
    
    enzymeCutsites <- get.enzyme.cutsites(sig, gen, gen_ver, chrs)
    RE_sites <- enzymeCutsites
    for (chr in chrs) {
        newends <- dplyr::lead(BiocGenerics::start(RE_sites[seqnames(RE_sites) == chr])) - 1
        newends <- c(newends[-length(newends)], genome.chromSizes[chr])
        BiocGenerics::end(RE_sites[seqnames(RE_sites) == chr]) <- newends
        minstart <- min(BiocGenerics::start(RE_sites[seqnames(RE_sites) == chr]))
        if (minstart > 1) {
            RE_sites <- GenomeInfoDb::sortSeqlevels(c(GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = 1, 
                end = minstart)), RE_sites))
            RE_sites <- sort(RE_sites)
        }
    }
    
    
    if (!(is.null(wg_file))) {
        wgdata <- rtracklayer::import.bw(wg_file, which = genome.chromGR, as = "GRanges")
    } else {
        wgdata <- NULL
    }
    
    if (bin_type == "Bins-RE-sites") {
        endL <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(start = GenomicRanges::start(RE_sites), 
            width = 200))
        for (chr in chrs) {
            BiocGenerics::end(endL[seqnames(endL) == chr]) <- pmin(BiocGenerics::end(endL[seqnames(endL) == chr]), genome.chromSizes[chr])
            BiocGenerics::start(endL[seqnames(endL) == chr]) <- pmax(BiocGenerics::start(endL[seqnames(endL) == chr]), 1)
        }
        endR <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(end = GenomicRanges::end(RE_sites), 
            width = 200))
        for (chr in chrs) {
            BiocGenerics::end(endR[seqnames(endL) == chr]) <- pmin(BiocGenerics::end(endR[seqnames(endR) == chr]), genome.chromSizes[chr])
            BiocGenerics::start(endR[seqnames(endL) == chr]) <- pmax(BiocGenerics::start(endR[seqnames(endR) == chr]), 1)
        }
        gcL <- getGC(regions = endL)
        gcR <- getGC(regions = endR)
        gc <- (gcL + gcR)/2
        RE_sites$gc <- gc
        
        if (!is.null(wgdata)) {
            endL <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(start = GenomicRanges::start(RE_sites), 
                width = 500))
            for (chr in chrs) {
                BiocGenerics::end(endL[seqnames(endL) == chr]) <- pmin(BiocGenerics::end(endL[seqnames(endL) == chr]), genome.chromSizes[chr])
                BiocGenerics::start(endL[seqnames(endL) == chr]) <- pmax(BiocGenerics::start(endL[seqnames(endL) == chr]), 1)
            }
            endR <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(end = GenomicRanges::end(RE_sites), 
                width = 500))
            for (chr in chrs) {
                BiocGenerics::end(endR[seqnames(endL) == chr]) <- pmin(BiocGenerics::end(endR[seqnames(endR) == chr]), genome.chromSizes[chr])
                BiocGenerics::start(endR[seqnames(endL) == chr]) <- pmax(BiocGenerics::start(endR[seqnames(endR) == chr]), 1)
            }
            mapL <- getMap(regions = endL, data = wgdata)
            mapR <- getMap(regions = endR, data = wgdata)
            map <- (mapL + mapR)/2
            RE_sites$map <- map
        } else {
            RE_sites$map <- 0
        }
        # generate bintolen from RE_sites
        bintolen <- as.data.frame(RE_sites, stringsAsFactors = FALSE)
        bintolen <- bintolen %>% dplyr::group_by(seqnames) %>% dplyr::mutate(RE_id = dplyr::row_number())
        if (binsize > 1) {
            # We cluster enzyme cutting sites every 'binsize' consecutive RE sites per each chromosome
            bintolen <- bintolen %>% dplyr::group_by(seqnames) %>% dplyr::mutate(binNumb = floor((dplyr::row_number() - 1)/binsize) + 
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
        gc()
        gc(reset = TRUE)
    }

    if (bin_type == "Bins-uniform") {
        bins.chrom <- function(chr, binsize) {
            cuts_start <- seq(1, genome.chromSizes[chr], by = binsize)
            cuts_end <- seq(binsize, genome.chromSizes[chr] - (genome.chromSizes[chr]%%binsize) + binsize, by = binsize)
            cuts_end <- c(cuts_end[-length(cuts_end)], genome.chromSizes[chr])
            bins <- GenomicRanges::GRanges(chr, IRanges::IRanges(start = cuts_start, end = cuts_end))
            return(bins)
        }
        binsGR <- BiocGenerics::unlist(GenomicRanges::GRangesList(lapply(chrs, function(x) {
            bins.chrom(x, binsize)
        })))
        names(binsGR) <- paste(as.character(GenomicRanges::seqnames(binsGR)), GenomicRanges::start(binsGR), GenomicRanges::end(binsGR), 
            sep = "-")
        
        bintolen <- list()
        i <- 1
        for (chr in chrs) {
            regions <- binsGR[GenomicRanges::seqnames(binsGR) == chr]
            cuts <- enzymeCutsites[GenomicRanges::seqnames(enzymeCutsites) == chr]
            medians <- (GenomicRanges::start(cuts) + GenomicRanges::end(cuts))/2
            if (!is.null(wgdata)) {
                data <- wgdata[GenomicRanges::seqnames(wgdata) == chr]
            } else {
                data <- NULL
            }
            # 500bp from ends for mappability and len features
            FragmentendsL <- GenomicRanges::GRanges(seqnames = seqnames(cuts), ranges = IRanges::restrict(IRanges::IRanges(end = medians - 
                1, width = 500), start = 1, end = genome.chromSizes[chr]))
            FragmentendsR <- GenomicRanges::GRanges(seqnames = seqnames(cuts), ranges = IRanges::restrict(IRanges::IRanges(start = medians, 
                width = 500), start = 1, end = genome.chromSizes[chr]))
            ends <- c(FragmentendsL, FragmentendsR)
            hits <- as.data.frame(GenomicRanges::findOverlaps(ends, regions, type = "within", select = "all"))
            LR <- data.frame(bins = names(regions[hits[[2]]]), start = GenomicRanges::start(ends[hits[[1]]]), end = GenomicRanges::end(ends[hits[[1]]]), 
                stringsAsFactors = FALSE)
            LRgr <- GenomicRanges::GRanges(chr, IRanges::IRanges(start = LR$start, end = LR$end))
            if (!is.null(wgdata)) {
                LR$map <- getMap(LRgr, data)
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
            FragmentendsL <- GenomicRanges::GRanges(seqnames = seqnames(cuts), ranges = IRanges::IRanges(end = medians - 1, 
                width = 200))
            FragmentendsR <- GenomicRanges::GRanges(seqnames = seqnames(cuts), ranges = IRanges::IRanges(start = medians, 
                width = 200))
            ends <- c(FragmentendsL, FragmentendsR)
            hits <- as.data.frame(GenomicRanges::findOverlaps(ends, regions, type = "within", select = "all"))
            LR2 <- data.frame(bins = names(regions[hits[[2]]]), start = GenomicRanges::start(ends[hits[[1]]]), end = GenomicRanges::end(ends[hits[[1]]]), 
                stringsAsFactors = FALSE)
            LR2gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(start = LR2$start, end = LR2$end))
            LR2$gc <- getGC(LR2gr)
            gc <- LR2 %>% dplyr::group_by(.data$bins) %>% dplyr::summarize(gc = mean(gc))
            LR <- dplyr::left_join(dplyr::left_join(gc, map), len)
            bintolen[[i]] <- as.data.frame(LR)
            i <- i + 1
            print(chr)
        }
        bintolen <- suppressWarnings(dplyr::bind_rows(bintolen))
        rownames(bintolen) <- bintolen$bins
        allbins <- data.frame(bins = names(binsGR), stringsAsFactors = FALSE)
        bintolen <- suppressWarnings(dplyr::left_join(allbins, bintolen) %>% 
                           tidyr::replace_na(list(gc = 0, map = 0, len = 0)))
    }
    
    bintolenoutput <- path.expand(paste0(output_path, "_bintolen.txt.gz"))
    bintolenoutputdir<-gsub("/[^/]+$", "",bintolenoutput)
    if (bintolenoutputdir==bintolenoutput){
        bintolenoutputdir<-gsub("\\[^\\]+$", "",bintolenoutput)
    }
    if (bintolenoutputdir==bintolenoutput){
        bintolenoutputdir<-gsub("\\\\[^\\\\]+$", "",bintolenoutput)
    }
    if (!bintolenoutputdir==bintolenoutput&!dir.exists(bintolenoutputdir)){
        dir.create(bintolenoutputdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
    }
    data.table::fwrite(bintolen, bintolenoutput, row.names = FALSE, quote = FALSE, sep = "\t")
    return(bintolenoutput)
}
