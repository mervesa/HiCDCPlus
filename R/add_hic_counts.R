#' add_hic_counts
#'
#'This function adds counts from a .hic file into a valid, binned, gi_list
#'instance.
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param gi_list valid, uniformly binned gi_list instance. 
#'See \code{?gi_list_validate} and \code{gi_list_binsize_detect} for details.
#'@param hic_path path to the .hic file
#'@param chrs a subset of chromosomes' e.g., c('chr21','chr22'). Defaults
#'to all chromosomes in the \code{gi_list} instance.
#'@param add_inter Interchromosomal interaction counts added as a 1D feature
#'named 'inter' on regions metadata handle of each gi_list element (e.g., 
#'\code{gi_list[[1]]@regions@elementMetadata} or not;
#'default FALSE
#'@return \code{gi_list} instance with counts on the metadata (e.g., 
#'\code{mcols(gi_list[[1]])} handle on each list element, and 'inter' on 
#'regions metadata handle of each element if \code{add_inter=TRUE}.
#'@examples gi_list<-generate_binned_gi_list(50e3,chrs='chr22')
#'gi_list<-add_hic_counts(gi_list,
#'hic_path=system.file("extdata", "GSE63525_HMEC_combined_example.hic",
#' package = "HiCDCPlus"))
#'@export

add_hic_counts <- function(gi_list, hic_path, chrs = NULL, add_inter = FALSE) {
    gi_list_validate(gi_list)
    binsize <- gi_list_binsize_detect(gi_list)
    Dthreshold <- gi_list_Dthreshold.detect(gi_list)
    if (is.null(chrs)) {
        chrs <- sort(names(gi_list))
    } else {
        chrs <- chrs[chrs %in% names(gi_list)]
        if (length(chrs) == 0) {
            stop("None of the chromosomes specified exists in the
                 gi_list object")
        }
    }
    for (chrom in chrs) {
        chr_select <- gsub("chr", "", chrom)
        if (!.Platform$OS.type=="windows"){
        count_matrix <- tryCatch(
            straw(norm = "NONE", fn = path.expand(hic_path), bs = binsize, ch1 = chr_select, ch2 = chr_select, 
            u = "BP"),
            error=function(e){
            tryCatch(straw(norm = "NONE", fn = path.expand(hic_path), bs = binsize, ch1 = chrom, ch2 = chrom, 
                      u = "BP"),
                     error=function(e){
                      straw_dump(norm = "NONE",fn=path.expand(hic_path),bs=binsize,ch1=chr_select,ch2=chr_select,u="BP")   
                     })
            })
        }else{
            count_matrix<-straw_dump(norm = "NONE",fn=path.expand(hic_path),bs=binsize,ch1=chr_select,ch2=chr_select,u="BP")   
        }
        count_matrix <- count_matrix %>% dplyr::filter(abs(.data$y - .data$x) <= Dthreshold) %>% dplyr::rename(startI = "x", 
            startJ = "y")
        if (nrow(count_matrix)==0){
            msg<-paste0(chrom, "does not have any counts in this file. Dropping from gi_list.")
            warning(msg)
            gi_list[[chrom]]<-NULL
            next
        }
        gi_list[[chrom]] <- add_2D_features(gi_list[[chrom]], count_matrix)
        msg<-paste0("Chromosome ",chrom," intrachromosomal counts processed.")
        message(msg)
    }
    rm(count_matrix, chr_select)
    if (add_inter) {
        # intercounts--hic files
        inter_matrix <- data.frame(stringsAsFactors = FALSE)
        for (chrom in chrs) {
            inter_matrix_add <- data.frame(stringsAsFactors = FALSE)
            for (chr2 in chrs) {
                if (chr2 == chrom) next
                chr_select <- gsub("chr", "", chrom)
                chr_other <- gsub("chr", "", chr2)
                if(!.Platform$OS.type=="windows"){
                count_matrix <- tryCatch(
                    straw(norm = "NONE", fn = path.expand(hic_path), bs = binsize, ch1 = chr_select, ch2 = chr_other, 
                          u = "BP"),
                    error=function(e){
                        tryCatch(straw(norm = "NONE", fn = path.expand(hic_path), bs = binsize, ch1 = chrom, ch2 = chr2, 
                              u = "BP"),
                              error=function(e){
                                  straw_dump(norm = "NONE",fn=path.expand(hic_path),bs=binsize,ch1=chr_select,ch2=chr_other,u="BP")   
                              })
                    })
                }else{
                    count_matrix<-straw_dump(norm = "NONE",fn=path.expand(hic_path),bs=binsize,ch1=chr_select,ch2=chr_other,u="BP")   
                    
                }
                count_matrix$chr <- chrom
                inter_matrix_add <- dplyr::bind_rows(inter_matrix_add, count_matrix %>% dplyr::group_by(.data$chr, .data$x) %>% 
                dplyr::summarize(inter = sum(.data$counts)))
                rm(count_matrix)
            }
            inter_matrix <- dplyr::bind_rows(inter_matrix, inter_matrix_add %>% dplyr::group_by(.data$chr, .data$x) %>% dplyr::summarize(inter = sum(.data$inter)))
            msg<-paste0("Chromosome ",chrom," interchromosomal counts processed.")
            message(msg)
        }
        inter_matrix <- inter_matrix %>% dplyr::ungroup() %>% dplyr::rename(start = "x")
        gi_list <- add_1D_features(gi_list, inter_matrix, chrs, agg = sum)
        rm(inter_matrix, inter_matrix_add)
    }
    return(gi_list)
}
