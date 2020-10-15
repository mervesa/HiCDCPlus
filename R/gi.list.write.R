#' gi.list.write.R
#'
#'Writes a valid gi_list instance into a file.
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param gi_list List of \code{GenomicInteractions} objects where each object
#' named with chromosomes contains intra-chromosomal interaction information
#'(see
#'\code{?gi.list.validate} for a detailed explanation of valid \code{gi_list}
#'instances). 
#'@param fname path to the file to write to (can end with .txt, or .txt.gz).
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes
#'in the \code{gi_list}.
#'@param columns Can be 'minimal', which is
#'just distance and counts (and \code{HiCDCPlus} result columns
#''qvalue','pvalue','mu',and 'sdev', if exists; see \code{?HiCDCPlus})
#'information, 'minimal_plus_features', which is
# distance, counts, and other calculated 2D features,
#''minimal_plus_score', which generates a .hic pre compatible text file,
#'or 'all', which is
#'distance, counts, calculated 2D features, as well as all 1D features.
#'Defaults to 'minimal'.
#'@param rows Can be 'all' or 'significant', which filters rows according to 
#'FDR adjusted pvalue column 'qvalue' (this has to exist in \code{mcols(.)})
#'at \code{significance_threshold}. Defaults to 'all'.
#'@param significance_threshold Row filtering threshold on 'qvalue'.
#'Defaults to 0.05.
#'@param score Score column to extract to .hic pre compatible file.
#'See \code{mode} options in \code{?hicdc2hic} for more details.
#'@return a tab separated flat file concatenating all intra-chromosomal
#'interaction information.
#'@examples 
#'outputdir<-paste0(tempdir(check=TRUE),'/')
#'gi_list<-generate.binned.gi.list(1e6,chrs='chr22')
#'gi.list.write(gi_list,paste0(outputdir,'test.txt'))
#'@export

gi.list.write <- function(gi_list, fname, chrs = NULL, columns = "minimal", rows = "all", significance_threshold = 0.05, score = NULL) {
    gi.list.validate(gi_list)
    if (is.null(chrs)) {
        chrs <- sort(names(gi_list))
    } else {
        chrs <- chrs[chrs %in% names(gi_list)]
        if (length(chrs) == 0) {
            stop("None of the chromosomes specified exists in the gi.list object")
        }
    }
    if (!(rows %in% c("all", "significant") & columns %in% c("minimal", "minimal_plus_features", "minimal_plus_score", "all") & 
        is.finite(significance_threshold) & significance_threshold >= 0 & significance_threshold <= 1)) {
        stop("Invalid option selected in one of the fields. Check ?gi.list.write
        for allowed options")
    }
    if (rows == "significant" & !all(vapply(gi_list[chrs], function(x) "qvalue" %in% names(mcols(x)), TRUE))) {
        stop("No FDR-adjusted p-value 'qvalue' column to extract significant
    interactions.")
    }
    if (!is.null(score)) {
        if (!score %in% c("raw", "pvalue", "qvalue", "zvalue", "normcounts")) {
            stop("Invalid option for mode/score. Allowable options are
            'raw','pvalue','qvalue','normcounts'")
        }
        if (score == "raw" & !all(vapply(gi_list[chrs], function(x) "counts" %in% names(mcols(x)), TRUE))) {
            stop("No counts column to extract scores for export.")
        }
        if (!score == "raw" & !all(vapply(gi_list[chrs], function(x) c("qvalue", "pvalue", "mu", "sdev") %in% names(mcols(x)), 
            c(TRUE, TRUE, TRUE, TRUE)))) {
            stop("No HiCDCPlus result columns to extract scores for export.")
        }
    }
    if (is.null("score") & columns == "minimal_plus_score") {
        stop("No score option specified for the column selection.")
    }
    write.flag = TRUE
    for (i in seq(length(chrs))) {
        chrom <- chrs[i]
        if (rows == "significant") {
            ix <- mcols(gi_list[[chrom]])$qvalue < significance_threshold
        }
        if (columns %in% c("minimal", "minimal_plus_features", "minimal_plus_score")) {
            infoI <- data.frame(InteractionSet::anchors(gi_list[[chrom]], type = "first"), stringsAsFactors = FALSE) %>% dplyr::select(.data$seqnames, 
                .data$start,.data$end) %>% dplyr::rename(chr = "seqnames") %>% dplyr::rename_all(function(x) paste0(x, "I"))
            infoJ <- data.frame(InteractionSet::anchors(gi_list[[chrom]], type = "second"), stringsAsFactors = FALSE) %>% 
                dplyr::select(.data$seqnames, .data$start,.data$end) %>% dplyr::rename(chr = "seqnames") %>% dplyr::rename_all(function(x) paste0(x, 
                "J"))
            if (columns == "minimal") {
                res_features <- c("counts", "pvalue", "qvalue", "mu", "sdev")
                res_features <- res_features[res_features %in% colnames(mcols(gi_list[[chrom]]))]
                if (length(res_features) == 0) {
                features <- (mcols(gi_list[[chrom]]))[c("D")]
                } else {
                features <- (mcols(gi_list[[chrom]]))[c("D", res_features)]
                }
                gi.out <- cbind(infoI, infoJ, features)
            }
            if (columns == "minimal_plus_features") {
                features <- mcols(gi_list[[chrom]])
                gi.out <- cbind(infoI, infoJ, features)
            }
            if (columns == "minimal_plus_score") {
                binsize <- gi.list.binsize.detect(gi_list)
                if (score == "raw") {
                scoredata <- mcols(gi_list[[chrom]])$counts
                }
                if (score == "pvalue") {
                scoredata <- pmax(-log10(mcols(gi_list[[chrom]])$pvalue), 0)
                }
                if (score == "qvalue") {
                scoredata <- pmax(-log10(mcols(gi_list[[chrom]])$qvalue), 0)
                }
                if (score == "zvalue") {
                scoredata <- (mcols(gi_list[[chrom]])$counts - mcols(gi_list[[chrom]])$mu)/mcols(gi_list[[chrom]])$sdev
                }
                if (score == "normcounts") {
                scoredata <- mcols(gi_list[[chrom]])$counts/mcols(gi_list[[chrom]])$mu
                }
                gi.out <- cbind(infoI, infoJ) %>% dplyr::mutate(score = scoredata, strI = 0, strJ = 0, fragI = 0, fragJ = 1, 
                startI = .data$startI + binsize/2, startJ = .data$startJ + binsize/2) %>% dplyr::select(.data$strI, .data$chrI, 
                .data$startI, .data$fragI, .data$strJ, .data$chrJ, .data$startJ, .data$fragJ, .data$score) %>% dplyr::filter(is.finite(.data$score) & 
                .data$score > 0)
            }
        } else {
            gi.out <- data.frame(gi_list[[chrom]]) %>% dplyr::rename(chr1 = "seqnames1", chr2 = "seqnames2") %>% dplyr::rename_all(function(x) gsub("2", 
                "J", gsub("1", "I", x))) %>% dplyr::select( -.data$strandI, -.data$strandJ, -.data$widthI, -.data$widthJ)
        }
        if (rows == "significant") {
            gi.out <- gi.out[ix, ]
        } 
        if (write.flag) {
            # first chromosome to be printed
            fname<-path.expand(fname)
            fnamedir<-gsub("/[^/]+$", "",fname)
            if (fnamedir==fname){
                fnamedir<-gsub("\\[^\\]+$", "",fname)
            }
            if (fnamedir==fname){
                fnamedir<-gsub("\\\\[^\\\\]+$", "",fname)
            }
            if (!fnamedir==fname&!dir.exists(fnamedir)){
                dir.create(fnamedir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            }
            if (.Platform$OS.type=="windows"){
                outdf<-list()
                outdf[[chrom]]<-data.frame(gi.out)
            }else{
                data.table::fwrite(data.frame(gi.out), fname, append = FALSE, sep = "\t", quote = FALSE, row.names = FALSE, verbose = FALSE)
            }
            write.flag=FALSE
        } else {
            if (.Platform$OS.type=="windows"){
                outdf[[chrom]]<-data.frame(gi.out)
            }else{
            data.table::fwrite(data.frame(gi.out), fname, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, verbose = FALSE)
            }
        }
    }
    if (.Platform$OS.type=="windows"){
        data.table::fwrite(dplyr::bind_rows(outdf), fname, append = FALSE, sep = "\t", quote = FALSE, row.names = FALSE, verbose = FALSE)
    }
    return(fname)
}
