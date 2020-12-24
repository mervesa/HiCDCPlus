#'add_hicpro_matrix.counts
#'
#'This function converts HiC-Pro matrix and bed outputs into a gi_list instance.
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param gi_list valid, uniformly binned gi_list instance. 
#'See \code{?gi_list_validate} and \code{gi_list_binsize_detect} for details.
#'@param absfile_path absfile BED out of HiC-Pro (e.g., 'rawdata_10000_abs.bed')
#'@param matrixfile_path matrix count file out of HiC-Pro (e.g.,
#''rawdata_10000.matrix')
#'@param chrs a subset of chromosomes' e.g., c('chr21','chr22'). Defaults
#'to all chromosomes in the \code{gi_list} instance.
#'@param add_inter Interchromosomal interaction counts added as a 1D feature
#'named 'inter' on regions metadata handle of each gi_list element (e.g., 
#'\code{gi_list[[1]]@regions@elementMetadata} or not;
#'default FALSE
#'@return \code{gi_list} instance with counts on the metadata (e.g., 
#'\code{mcols(gi_list[[1]])} handle on each list element, and 'inter' on 
#'regions metadata handle of each element if \code{add_inter=TRUE}.
#'@export

add_hicpro_matrix_counts <- function(gi_list, absfile_path, matrixfile_path, chrs = NULL, add_inter = FALSE) {
    gi_list_validate(gi_list)
    if (is.null(chrs)) {
        chrs <- sort(names(gi_list))
    } else {
        chrs <- chrs[chrs %in% names(gi_list)]
        if (length(chrs) == 0) {
            stop("None of the chromosomes specified exists in the gi.list object")
        }
    }
    Dthreshold <- gi_list_Dthreshold.detect(gi_list)
    absI <- data.table::fread(absfile_path, sep = "\t", header = FALSE) %>% dplyr::rename(chrI = "V1", startI = "V2", end = "V3", 
        first = "V4") %>% dplyr::select(-.data$end)
    absJ <- data.table::fread(absfile_path, sep = "\t", header = FALSE) %>% dplyr::rename(chrJ = "V1", startJ = "V2", end = "V3", 
        second = "V4") %>% dplyr::select(-.data$end)
    count_matrix <- data.table::fread(matrixfile_path, sep = "\t", stringsAsFactors = FALSE) %>% dplyr::rename(first = "V1", 
        second = "V2", counts = "V3")
    count_matrix <- dplyr::left_join((dplyr::left_join(as.data.frame(count_matrix), absI) %>% dplyr::filter(!is.na(.data$chrI))), 
        absJ) %>% dplyr::filter(!is.na(.data$chrJ)) %>% dplyr::filter(.data$chrI == .data$chrJ & abs(.data$startJ - .data$startI) <= 
        Dthreshold) %>% dplyr::select(.data$chrI, .data$startI, .data$startJ, .data$counts) %>% dplyr::rename(chr = "chrI")
    rm(absI, absJ)
    for (chrom in chrs) {
        count_matrix_chr <- count_matrix %>% dplyr::filter(.data$chr == chrom)
        gi_list[[chrom]] <- add_2D_features(gi_list[[chrom]], count_matrix_chr)
        print(paste0("Chromosome ",chrom," intrachromosomal counts processed."))
    }
    rm(count_matrix_chr, count_matrix)
    # intercounts--hicpro files
    if (add_inter) {
        absI <- data.table::fread(absfile_path, sep = "\t", header = FALSE) %>% dplyr::rename(chrI = "V1", startI = "V2", 
            end = "V3", first = "V4") %>% dplyr::select(-.data$end)
        absJ <- data.table::fread(absfile_path, sep = "\t", header = FALSE) %>% dplyr::rename(chrJ = "V1", startJ = "V2", 
            end = "V3", second = "V4") %>% dplyr::select(-.data$end)
        count_matrix <- data.table::fread(matrixfile_path, sep = "\t", stringsAsFactors = FALSE) %>% dplyr::rename(first = "V1", 
            second = "V2", counts = "V3")
        count_matrix <- dplyr::left_join((dplyr::left_join(as.data.frame(count_matrix), absI) %>% dplyr::filter(!is.na(.data$chrI))), 
            absJ) %>% dplyr::filter(!is.na(.data$chrJ)) %>% dplyr::filter(!.data$chrI == .data$chrJ)
        count_matrix <- dplyr::bind_rows(count_matrix %>% dplyr::group_by(.data$chrI, .data$startI) %>% dplyr::summarize(counts = sum(.data$counts)) %>% 
            dplyr::rename(chr = "chrI", start = "startI"), count_matrix %>% dplyr::group_by(.data$chrJ, .data$startJ) %>% 
            dplyr::summarize(counts = sum(.data$counts)) %>% dplyr::rename(chr = "chrJ", start = "startJ"))
        count_matrix <- count_matrix %>% dplyr::group_by(.data$chr, .data$start) %>% dplyr::summarize_all(sum) %>% dplyr::ungroup() %>% 
            dplyr::rename(inter = "counts")
        rm(absI, absJ)
        gi_list <- add_1D_features(gi_list, count_matrix, chrs, agg = sum)
        rm(count_matrix)
    }
    return(gi_list)
}
