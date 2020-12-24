#' generate_df_gi_list
#'
#'Generates a gi_list instance from a data frame object describing the regions.
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param df DataFrame with columns named 'chr', 'start', (and optionally
#''end', if the regions have gaps) and 1D features with their respective
#'column names.
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes
#'specified in \code{df}.
#'@param Dthreshold maximum distance (included) to check for significant
#'interactions, defaults to 2e6 or maximum in the data, whichever is smaller.
#'@param gen name of the species: e.g., default \code{'Hsapiens'}
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}
#'@return a valid gi_list instance with genomic features supplied from 
#'\code{df}. Genomic 1D features are stored in the regions metadata handle
#'of each list element (e.g., \code{gi_list[[1]]@regions@elementMetadata}).
#'@examples 
#'df<-data.frame(chr='chr9',start=seq(1e6,10e6,1e6))
#'gi_list<-generate_df_gi_list(df)
#'@export

generate_df_gi_list <- function(df, chrs = NULL, Dthreshold = 2e+06,
                                gen = "Hsapiens", gen_ver = "hg19") {
    gi_list <- list()
    if (!("chr" %in% colnames(df) & "start" %in% colnames(df))) {
        stop("the df should have columns named 'chr' and 'start'")
    }
    if (!sum(stats::complete.cases(df%>%
            dplyr::select(.data$chr,.data$start))) == nrow(df)) {
        stop("the df has NULL/NA values for 'chr','start'")
    }
    if (!(sum(is.finite(df$start)) == nrow(df))) {
        stop("the df has NA/NaN/Inf values for 'start'")
    }
    if (is.null(chrs)) 
        chrs <- sort(unique(df$chr))
    if ('end'%in%colnames(df)){
    df <- df %>% dplyr::arrange(.data$chr, .data$start,.data$end)
        
    }else{
    df <- df %>% dplyr::arrange(.data$chr, .data$start)
    }
    for (chrom in chrs) {
        df_chr <- df %>% dplyr::filter(.data$chr == chrom)
        if (!'end'%in%colnames(df_chr)){
        df_chr<-df_chr%>%dplyr::mutate(end=c(.data$start[-1],
            get_chr_sizes(gen,gen_ver,chrom)[chrom]))
        }
        all.regions <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = df_chr$start, end = df_chr$end))
        eff_binsize <- stats::quantile(df_chr$end - df_chr$start, probs = 0.01)
        numbins <- length(all.regions)
        maxbins <- min(round(Dthreshold/eff_binsize), numbins)
        index1 <- unlist(lapply(seq(1, numbins, 1), function(x) rep(x, min(maxbins + 1, numbins - x + 1))))
        index2 <- unlist(lapply(seq(1, numbins, 1), function(x) seq(x, min(x + maxbins, numbins), 1)))
        gi_list[[chrom]] <- InteractionSet::GInteractions(index1, index2, all.regions)
        mcols(gi_list[[chrom]])$D <- InteractionSet::pairdist(gi_list[[chrom]])
        gi_list[[chrom]] <- gi_list[[chrom]][mcols(gi_list[[chrom]])$D <= Dthreshold]
    }
    return(gi_list)
}
