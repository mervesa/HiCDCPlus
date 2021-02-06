#' HiCDCPlus_parallel
#'
#' This function finds significant interactions in a HiC-DC readable matrix
#' and expresses statistical significance of counts through the following with
#' a parallel implementation (using sockets; compatible with Windows):
#' 'pvalue': significance \emph{P}-value, 'qvalue': FDR corrected 
#' \emph{P}-value, mu': expected counts, 'sdev': modeled standard deviation
#' of expected counts.
#'@import BSgenome splines
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@importFrom S4Vectors mcols<- mcols queryHits subjectHits
#'@param gi_list List of \code{GenomicInteractions} objects where each object
#'named with chromosomes contains intrachromosomal interaction information
#'(minimally containing counts and genomic distance in 
#'\code{mcols(gi_list[[1]])}---see
#'\code{?gi_list_validate} for a detailed explanation of valid \code{gi_list}
#'instances). 
#'@param covariates covariates to be considered in addition to genomic
#'distance D. Defaults to all covariates besides 
#''D','counts','mu','sdev',pvalue','qvalue' 
#'in \code{mcols(gi)}
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes
#'in the \code{gi_list}.
#'@param distance_type distance covariate form: 'spline' or 'log'.
#'Defaults to 'spline'.
#'@param model_distribution 'nb' uses a Negative Binomial model, 
#''nb_vardisp' uses a Negative Binomial model with a distance specific 
#'dispersion parameter inferred from the data, 'nb_hurdle' uses the legacy
#'HiC-DC model.
#'@param binned TRUE if uniformly binned or FALSE if binned by
#'restriction enzyme fragment cutsites
#'@param df degrees of freedom for the genomic distance spline
#'function if \code{distance_type='spline'}. Defaults to 6, which corresponds to
#'a cubic spline as explained in Carty et al. (2017)
#'@param Dmin minimum distance (included) to check for significant interactions,
#' defaults to 0
#'@param Dmax maximum distance (included) to check for significant interactions,
#' defaults to 2e6 or maximum in the data; whichever is minimum.
#'@param ssize Distance stratified sampling size. Can decrease for
#' large chromosomes. Increase recommended if
#' model fails to converge. Defaults to 0.01.
#'@param splineknotting Spline knotting strategy. Either "uniform", uniformly
#'spaced in distance, or placed based on distance distribution of counts 
#'"count-based" (i.e., more closely spaced where counts are more dense).
#'@param ncore Number of cores to parallelize. Defaults to 
#'\code{parallel::detectCores()-1}.
#'@return A valid \code{gi_list} instance with additional \code{mcols(.)} for
#'each chromosome: pvalue': significance \emph{P}-value, 'qvalue': FDR 
#'corrected \emph{P}-value, mu': expected counts, 'sdev': modeled standard
#'deviation of expected counts.
#'@examples gi_list<-generate_binned_gi_list(50e3,chrs='chr22')
#'gi_list<-add_hic_counts(gi_list,
#'hic_path=system.file("extdata", "GSE63525_HMEC_combined_example.hic",
#' package = "HiCDCPlus"))
#'gi<-HiCDCPlus_parallel(gi_list,ncore=1)
#'@export

HiCDCPlus_parallel <- function(gi_list, covariates = NULL, chrs = NULL, distance_type = "spline", model_distribution = "nb", binned = TRUE, 
    df = 6, Dmin = 0, Dmax = 2e+06, ssize = 0.01, splineknotting = "uniform", ncore=NULL) {
    options(scipen = 9999, digits = 4)
    gi_list_validate(gi_list)
    if (is.null(chrs)) chrs <- names(gi_list)
    # check if D and counts exist on each chromosome
    if (!(all(vapply(gi_list, function(x) sum(colnames(S4Vectors::mcols(x)) %in% c("counts", "D")) == 2, TRUE)))) {
        stop("Some gi_list elements do not contain pairwise distances D
    and counts, the minimum set of features needed for HiC-DC+ modeling.")
    }
    if (!model_distribution %in% c("nb", "nb_hurdle", "nb_vardisp")) {
        stop("Allowable options for model_distribution are 'nb','nb_hurdle', and
    nb_vardisp'.")
    }
    if (is.null(ncore)) ncore<-parallel::detectCores()-1
    cl <- parallel::makeCluster(ncore)
    parallel::clusterEvalQ(cl,library("dplyr"))
    parallel::clusterEvalQ(cl,library("rlang"))
    gi_list[chrs] <- parallel::parLapply(cl, gi_list[chrs], HiCDCPlus_chr,covariates=covariates,
                         distance_type=distance_type,
                         model_distribution=model_distribution,
                         binned=binned,
                         df=df,
                         Dmin=Dmin,
                         Dmax=Dmax,
                         ssize=ssize,
                         splineknotting = splineknotting,
                         model_filepath = NULL)
    parallel::stopCluster(cl)
    return(gi_list)
}
