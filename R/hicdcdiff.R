#' hicdcdiff
#'
#'This function calculates differential interactions for a
#'set of chromosomes across conditions and replicates. You need to install
#'\code{DESeq2} from Bioconductor to use this function.
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param input_paths a list with names as condition names and values
#'as paths to \code{gi_list} RDS objects (see \code{?gi_list_validate}
#'for a detailed explanation of valid \code{gi_list} instances) 
#'saved with \code{saveRDS} or
#'paths to .hic files for each replicate. e.g.,\code{list(
#'CTCF=c('~/Downloads/GM_CTCF_rep1_MAPQ30_10kb.rds',
#''~/Downloads/GM_CTCF_rep2_MAPQ30_10kb.rds'),
#'SMC=c('~/Downloads/GM_SMC_rep1_MAPQ30_10kb.rds',
#''~/Downloads/GM_SMC_rep2_MAPQ30_10kb.rds'))}
#'@param filter_file path to the text file containing columns 
#'chr', startI, and
#'startJ denoting the name of the chromosomes and starting coordinates
#'of 2D interaction bins to be compared across conditions, respectively.
#'@param output_path the path to the folder and name prefix you want to
#'place DESeq-processed matrices (in a .txt file), plots
#'(if \code{diagnostics=TRUE}) and DESeq2 objects (if \code{DESeq.save=TRUE}).
#'Files will be generated for each chromosome.
#'@param bin_type 'Bins-uniform' if uniformly binned by binsize in
#'bp, or 'Bins-RE-sites' if binned by number of
#'restriction enzyme fragment cutsites!
#'@param binsize binsize in bp if bin_type='Bins-uniform' (or number of
#'RE fragments if bin_type='Bins-RE-sites'), e.g., default 5000
#'@param granularity Desired distance granularity to base dispersion parameters
#'on in bp. For uniformly binned analysis 
#'(i.e., \code{bin_type=='Bins-uniform'}),
#'this defaults to the bin size. Otherwise, it is 5000.
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes (except Y and M)
#'in the filter_file.
#'@param Dmin minimum distance (included) to check for significant interactions,
#'defaults to 0. Put Dmin=1 to ignore D=0 bins in calculating
#'normalization factors.
#'@param Dmax maximum distance (included) to check for significant interactions,
#' defaults to 2e6 or maximum in the data; whichever is minimum.
#'@param diagnostics if TRUE, generates diagnostic plots of the normalization
#'factors, geometric means of such factors by distance bin, as well as MA Plots
#'(see DESeq documentation for details about MA plots). Defaults to FALSE.
#'@param DESeq.save if TRUE, saves the DESeq objects for each chromosome
#'as an .rds file in the \code{output_path}. Defaults to FALSE.
#'@param fitType follows fitType in \code{DESeq2::estimateDispersions}.
#'Allowable options are 'parametric' (parametric regression),'local' 
#'(local regression), and 'mean' (constant across interaction bins). 
#'Default is 'local'.
#'@return paths of a list of three entities.
#'\code{outputpaths} will have differential bins among those in filter_file.
#'\code{deseq2paths} will have the DESeq2 object stored as an .rds file.
#'Available if \code{DESeq.save=TRUE}
#'\code{plotpaths} will have diagnostic plots (e.g., MA, dispersion, PCA)
#'if \code{diagnostics=TRUE}.
#'@examples 
#'outputdir<-paste0(tempdir(check=TRUE),'/')
#'hicdcdiff(input_paths=list(NSD2=c(
#'system.file("extdata", "GSE131651_NSD2_LOW_arima_example.hic",
#' package = "HiCDCPlus"),
#' system.file("extdata", "GSE131651_NSD2_HIGH_arima_example.hic",
#' package = "HiCDCPlus")),
#' TKO=c(system.file("extdata", "GSE131651_TKOCTCF_new_example.hic",
#' package = "HiCDCPlus"),
#' system.file("extdata", "GSE131651_NTKOCTCF_new_example.hic",
#' package = "HiCDCPlus"))),
#' filter_file=system.file("extdata", "GSE131651_analysis_indices.txt.gz",
#' package = "HiCDCPlus"),
#'          chrs='chr22',
#'          output_path=outputdir,
#'          fitType = 'mean',
#'          binsize=50000,
#'          diagnostics=FALSE)
#'@export

hicdcdiff <- function(input_paths, filter_file, output_path, bin_type = "Bins-uniform", binsize = 5000, granularity = 5000, 
    chrs = NULL, Dmin = 0, Dmax = 2e+06, diagnostics = FALSE, DESeq.save = FALSE, fitType = "local") {
    options(scipen = 100, digits = 4, warn = -1)
    conds <- names(input_paths)
    granularity <- ifelse(bin_type == "Bins-uniform" & granularity == 5000, binsize, granularity)
    dband <- seq(0, Dmax, granularity)
    dband[1] <- Dmin
    dband <- unique(sort(dband[dband >= Dmin]))
    deseq2paths <- NULL
    outputpaths <- NULL
    plotpaths <- NULL
    sigs <- data.table::fread(filter_file)
    #set memory limit to max if i386
    if (.Platform$OS.type=='windows'&Sys.getenv("R_ARCH")=="/i386") {
        gc(reset=TRUE,full=TRUE)
        utils::memory.limit(size=4095)
    }
    # set default chrs if need be
    if (is.null(chrs)) {
        # get list of chromosomes from sigs
        if (!"chr" %in% colnames(sigs)) 
            stop("No column named 'chr' in filter_file.")
        chrs <- sort(unique(sigs$chr))
    }

    for (chr in chrs) {
        plotpaths_add <- NULL
        outputpaths_add <- NULL
        deseq2paths_add <- NULL
        retlist <- .readDiffInputFiles(conds, input_paths, sigs, chr, Dmin, Dmax, dband, bin_type, binsize)
        if (diagnostics) {
            plotpaths_add <- .plotNormalizationFactors(retlist[["normfac.final"]], binsize, chr, dband, conds, output_path)
        }
        countcols <- c()
        for (cond in conds) {
            countcols <- c(countcols, paste0(cond, ".", seq(length(input_paths[[cond]]))))
        }
        count.mat <- retlist[["normfac"]][countcols]
        rownames(count.mat) <- paste0(retlist[["normfac"]]$startI, ":", retlist[["normfac"]]$startJ)
        normcols <- paste0(countcols, ".norm")
        normfac.mat <- retlist[["normfac"]][normcols]
        rownames(normfac.mat) <- rownames(count.mat)
        sampleinfo <- data.frame(condition = rep(conds, times = unlist(lapply(input_paths, length))))
        rownames(sampleinfo) <- colnames(count.mat)
        sampleinfo$samples <- rownames(sampleinfo)
        count.mat <- as.matrix(count.mat)
        normfac.mat <- as.matrix(normfac.mat)
        ## run deseq2
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.mat, colData = sampleinfo, design = ~condition)
        dds$condition <- stats::relevel(dds$condition, ref = conds[1])
        DESeq2::normalizationFactors(dds) <- normfac.mat
        keep <- rownames(DESeq2::counts(dds)) %in% retlist[["sigbins"]]
        dds <- dds[keep, ]
        dds <- DESeq2::estimateDispersionsGeneEst(dds)
        DESeq2::dispersions(dds) <- tryCatch({
            dds <- DESeq2::estimateDispersionsFit(dds, fitType = fitType)
            GenomicRanges::mcols(dds)$dispFit},
            error=function(e) {
                msg<-paste0("fitType failed to fit for ",chr,". Overriding with gene estimates")
                message(msg)
                return(mcols(dds)$dispGeneEst)
            })
        dds <- DESeq2::nbinomWaldTest(dds)
        for (i in seq(1, length(conds) - 1, 1)) {
            for (j in seq(i + 1, length(conds), 1)) {
                # create contrast variables
                var <- paste0("res", conds[j], "over", conds[i])
                assign(var, DESeq2::results(dds, contrast = c("condition", conds[j], conds[i])))
                if (diagnostics) {
                # MA plot
                maplot_path <- path.expand(paste0(output_path, "plotMA_", conds[j], "over", conds[i], "_", chr, ".pdf"))
                maplot_pathdir<-gsub("/[^/]+$", "",maplot_path)
                if (maplot_pathdir==maplot_path){
                    maplot_pathdir<-gsub("\\[^\\]+$", "",maplot_path)
                }
                if (maplot_pathdir==maplot_path){
                    maplot_pathdir<-gsub("\\\\[^\\\\]+$", "",maplot_path)
                }
                if (!maplot_pathdir==maplot_path&!dir.exists(maplot_pathdir)){
                    dir.create(maplot_pathdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
                }
                grDevices::pdf(maplot_path)
                DESeq2::plotMA(get(var))
                grDevices::dev.off()
                plotpaths_add <- c(plotpaths_add, maplot_path)
                }
                # write to file
                assign(var, as.data.frame(get(var)))
                assign(var, get(var) %>% tibble::rownames_to_column("bin") %>% dplyr::mutate(chr = chr) %>% tidyr::separate(col = .data$bin, 
                into = c("startI", "startJ"), sep = ":", remove = TRUE, convert = TRUE))
                assign(var, get(var)[c("chr", setdiff(names(get(var)), "chr"))])
                filepath <- paste0(output_path, "diff_", var, "_", chr, ".txt.gz")
                filepath<-path.expand(filepath)
                filepathdir<-gsub("/[^/]+$", "",filepath)
                if (filepathdir==filepath){
                    filepathdir<-gsub("\\[^\\]+$", "",filepath)
                }
                if (filepathdir==filepath){
                    filepathdir<-gsub("\\\\[^\\\\]+$", "",filepath)
                }
                if (!filepathdir==filepath&!dir.exists(filepathdir)){
                    dir.create(filepathdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
                }
                data.table::fwrite(get(var), filepath, quote = FALSE, sep = "\t", row.names = FALSE)
                outputpaths_add <- c(outputpaths_add, filepath)
            }
        }
        
        if (diagnostics) {
            # PCA plot
            rld <- DESeq2::rlog(dds, fitType = fitType)
            pcaplot_path <- path.expand(paste0(output_path, "diff_", chr, "_PCA.pdf"))
            pcaplot_pathdir<-gsub("/[^/]+$", "",pcaplot_path)
            if (pcaplot_pathdir==pcaplot_path){
                pcaplot_pathdir<-gsub("\\[^\\]+$", "",pcaplot_path)
            }
            if (pcaplot_pathdir==pcaplot_path){
                pcaplot_pathdir<-gsub("\\\\[^\\\\]+$", "",pcaplot_path)
            }
            if (!pcaplot_pathdir==pcaplot_path&!dir.exists(pcaplot_pathdir)){
                dir.create(pcaplot_pathdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            }
            grDevices::pdf(pcaplot_path)
            p <- DESeq2::plotPCA(rld, intgroup = "samples")
            print(p)
            grDevices::dev.off()
            plotpaths_add <- c(plotpaths_add, pcaplot_path)
            # Dispersion plot
            dispplot_path <- path.expand(paste0(output_path, "dispersionplot.pdf"))
            dispplot_pathdir<-gsub("/[^/]+$", "",dispplot_path)
            if (dispplot_pathdir==dispplot_path){
                dispplot_pathdir<-gsub("\\[^\\]+$", "",dispplot_path)
            }
            if (dispplot_pathdir==dispplot_path){
                dispplot_pathdir<-gsub("\\\\[^\\\\]+$", "",dispplot_path)
            }
            if (!dispplot_pathdir==dispplot_path&!dir.exists(dispplot_pathdir)){
                dir.create(dispplot_pathdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            }
            grDevices::pdf(dispplot_path)
            DESeq2::plotDispEsts(dds)
            grDevices::dev.off()
            plotpaths_add <- c(plotpaths_add, dispplot_path)
        }
        
        # save DESeq file
        if (DESeq.save) {
            deseq2path <- paste0(output_path, chr, "_DESeq2_obj.rds")
            deseq2output <- path.expand(deseq2path)
            deseq2outputdir<-gsub("/[^/]+$", "",deseq2output)
            if (deseq2outputdir==deseq2output){
            deseq2outputdir<-gsub("\\[^\\]+$", "",deseq2output)
            }
            if (deseq2outputdir==deseq2output){
            deseq2outputdir<-gsub("\\\\[^\\\\]+$", "",deseq2output)
            }
            if (!deseq2outputdir==deseq2output&!dir.exists(deseq2outputdir)){
            dir.create(deseq2outputdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            }
            saveRDS(dds, deseq2path)
            deseq2paths_add <- c(deseq2paths_add, deseq2path)
        }
        deseq2paths <- c(deseq2paths, deseq2paths_add)
        outputpaths <- c(outputpaths, outputpaths_add)
        plotpaths <- c(plotpaths, plotpaths_add)
    }
    return(list(deseq2paths = deseq2paths, outputpaths = outputpaths, plotpaths = plotpaths))
}
