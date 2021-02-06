#' HiCDCPlus_chr
#'
#' This function finds significant interactions in a HiC-DC readable matrix
#' restricted to a single chromosome and expresses statistical significance of 
#' counts through the following:
#' 'pvalue': significance \emph{P}-value, 'qvalue': FDR corrected 
#' \emph{P}-value, mu': expected counts, 'sdev': modeled standard deviation
#' of expected counts.
#'@import BSgenome splines
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@importFrom S4Vectors mcols<- mcols queryHits subjectHits
#'@param gi Instance of a single chromosome \code{GenomicInteractions} object
#' containing intra-chromosomal interaction information
#' (minimally containing counts and genomic distance). 
#'@param covariates covariates to be considered in addition to genomic
#'distance D. Defaults to all covariates besides 
#''D','counts','mu','sdev',pvalue','qvalue' 
#'in \code{mcols(gi)}
#'@param distance_type distance covariate form: 'spline' or 'log'.
#'Defaults to 'spline'.
#'@param model_distribution 'nb' uses a Negative Binomial model, 
#''nb_vardisp' uses a Negative Binomial model with a distance specific 
#'dispersion parameter inferred from the data, 'nb_hurdle' uses the legacy
#'HiC-DC model.
#'@param binned TRUE if uniformly binned or FALSE if binned by
#'restriction enzyme fragment cut sites.
#'@param df degrees of freedom for the genomic distance spline
#'function if \code{distance_type='spline'}. Defaults to 6, which corresponds to
#'a cubic spline as explained in Carty et al. (2017)
#'@param Dmin minimum distance (included) to check for significant interactions,
#'defaults to 0
#'@param Dmax maximum distance (included) to check for significant interactions,
#'defaults to 2e6 or maximum in the data; whichever is minimum.
#'@param ssize Distance stratified sampling size. Can decrease for
#'large chromosomes. Increase recommended if
#'model fails to converge. Defaults to 0.01.
#'@param splineknotting Spline knotting strategy. Either "uniform", uniformly
#'spaced in distance, or placed based on distance distribution of counts 
#'"count-based" (i.e., more closely spaced where counts are more dense).
#'@param model_filepath Outputs fitted HiC-DC model object as an .rds
#'file with chromosome name indicatd on it. Defaults to NULL (no output).
#'@return A valid \code{gi} instance with additional \code{mcols(.)}: 
#'pvalue': significance \emph{P}-value, 'qvalue': FDR 
#'corrected \emph{P}-value, mu': expected counts, 'sdev': modeled standard
#'deviation of expected counts.
#'@examples gi_list<-generate_binned_gi_list(50e3,chrs='chr22')
#'gi_list<-add_hic_counts(gi_list,
#'hic_path<-system.file("extdata", "GSE63525_HMEC_combined_example.hic",
#' package = "HiCDCPlus"))
#'gi<-HiCDCPlus_chr(gi_list[[1]])
#'@export

HiCDCPlus_chr <- function(gi, covariates = NULL, distance_type = "spline", model_distribution = "nb", binned = TRUE, 
    df = 6, Dmin = 0, Dmax = 2e+06, ssize = 0.01, splineknotting = "uniform", model_filepath=NULL) {
    options(scipen = 9999, digits = 4)
    # remove logD, pvalue, qvalue, mu and sdev if they exist
    col_rem <- names(S4Vectors::mcols(gi))
    col_rem <- col_rem[col_rem %in% c("pvalue", "qvalue", "mu", "sdev", "logD")]
    for (covar in col_rem) {
        S4Vectors::mcols(gi)[, covar] <- NULL
    }
    # set default covariates if need be: the minimum set of features available across chromosomes other than D and counts
    if (is.null(covariates)) {
        tally <- as.data.frame(base::table(colnames(S4Vectors::mcols(gi))), 
            stringsAsFactors = FALSE)
        covariates <- as.character(tally$Var1[tally$Freq == 1])
        covariates <- covariates[!covariates %in% c("counts", "D", "pvalue", "qvalue", "mu", "sdev", "logD")]
        if (length(covariates) == 0) 
            covariates <- NULL
    }
    if (!distance_type %in% c("spline","log")) stop("Valid options for distance_type are 'spline' or 'log'. See ?HiCDCPlus_chr for details.")
    if (!model_distribution %in% c("nb","nb_hurdle", "nb_vardisp")) stop("Valid options for model_distribution are 'nb','nb_hurdle', or 'nb_vardisp'. See ?HiCDCPlus_chr for details.")
    if (!splineknotting %in% c("uniform","count-based")) stop("Valid options for splineknotting are 'uniform' or 'count-based'. See ?HiCDCPlus_chr for details.")

        # filter into defined covariate rows
        if ('len'%in%(colnames(S4Vectors::mcols(gi)))){
            gi<-gi[S4Vectors::mcols(gi)$len!=0]
        }
        if ('width'%in%(colnames(S4Vectors::mcols(gi)))){
            gi<-gi[S4Vectors::mcols(gi)$width!=0]
        }
        # get distance eligible row indices
        D.eligible <- S4Vectors::mcols(gi)$D >= Dmin & S4Vectors::mcols(gi)$D <= Dmax
        if (distance_type == "spline") {
            bdpts <- range(S4Vectors::mcols(gi)$D)
            if (splineknotting == "uniform"){
                knots <- NULL
            } else {
                countsums<-as.data.frame(S4Vectors::mcols(gi)[D.eligible,])%>%dplyr::group_by(.data$D)%>%dplyr::summarize(counts=sum(.data$counts))%>%dplyr::filter(.data$counts>0)
                knots<-.weighted_quantile(x=countsums$D,w=countsums$counts,num=df-3)
                knots<-knots[!knots%in%bdpts]
            }
        } else {
            #new.x$logD <- log2(new.x$D + 1)
            S4Vectors::mcols(gi)$logD <- log2(S4Vectors::mcols(gi)$D + 1)
        }
        if (model_distribution == "nb_vardisp") {
            # get thetas for each D.range
            dat <- as.data.frame(S4Vectors::mcols(gi))%>%dplyr::select(.data$D,.data$counts)
            dat <- dat %>% dplyr::mutate(D.range = findInterval(.data$D, unique(c(seq(Dmin, min(Dmax, 1e+06), by = 50000), 
                seq(min(Dmax, 1e+06), min(Dmax, 2e+06), by = 1e+05))), rightmost.closed = TRUE))
            dispersion_DF <- data.frame(stringsAsFactors = FALSE)
            for (x in sort(unique(dat$D.range))) {
                zeroPairedBins_ix <- seq(nrow(dat))[dat$D.range == x & dat$counts == 0]
                countPairedBins_ix <- seq(nrow(dat))[dat$D.range == x & dat$counts != 0]
                dat_x_ix <- c(countPairedBins_ix[sample(length(countPairedBins_ix), floor(length(countPairedBins_ix) * ssize), 
                replace = FALSE)], zeroPairedBins_ix[sample(length(zeroPairedBins_ix), floor(length(zeroPairedBins_ix) * ssize), replace = FALSE)])
                bdpts_x <- range(S4Vectors::mcols(gi)[dat_x_ix,]$D)
                if (splineknotting == "uniform"){
                    knots_x <- NULL
                } else {
                    countsums_x<-as.data.frame(S4Vectors::mcols(gi)[dat_x_ix,])%>%dplyr::group_by(.data$D)%>%dplyr::summarize(counts=sum(.data$counts))%>%dplyr::filter(.data$counts>0)
                    knots_x<-.weighted_quantile(x=countsums_x$D,w=countsums_x$counts,num=df-3)
                    knots_x<-knots_x[!knots_x%in%bdpts_x]
                }
                fit <- suppressWarnings(.GLM(data=as.data.frame(S4Vectors::mcols(gi)[dat_x_ix,]), df=df, knots=knots_x, bdpts=bdpts_x, covariates=covariates, distance_type = distance_type, model_distribution="nb", splineknotting=splineknotting))
                dt <- data.frame(D = mean(S4Vectors::mcols(gi)[dat_x_ix,]$D), theta = fit$theta)
                dispersion_DF <- dplyr::bind_rows(dispersion_DF, dt)
            }
            alpha <- .dispersionfunction(dispersion_DF)
            rm(dat, dat_x_ix, bdpts_x, knots_x, fit, dt, zeroPairedBins_ix, countPairedBins_ix)
        }
        # get a stratified sample for modeling
        dat <- as.data.frame(S4Vectors::mcols(gi))%>%dplyr::select(.data$D,.data$counts)
        if (binned) {
            binsize <- .get_range(gi)
            dat<-dat%>%dplyr::mutate(D.range = findInterval(.data$D, seq(Dmin, to = Dmax, by = binsize), rightmost.closed = TRUE),
                                     ix=seq(nrow(dat)))
            bins.counts <- split((dat%>%dplyr::filter(.data$counts>0))$ix, (dat%>%dplyr::filter(.data$counts>0))$D.range)
            bins.zeros <- split((dat%>%dplyr::filter(.data$counts==0))$ix, (dat%>%dplyr::filter(.data$counts==0))$D.range)
            rm(dat)
            idx.counts <- unlist(lapply(bins.counts, function(x) {
                sample(x, floor(length(x) * ssize), replace = FALSE)
            }))
            idx.zeros <- unlist(lapply(bins.zeros, function(x) {
                sample(x, floor(length(x) * ssize), replace = FALSE)
            }))
            dat <- dplyr::bind_rows(as.data.frame(S4Vectors::mcols(gi))[idx.counts, ], as.data.frame(S4Vectors::mcols(gi))[idx.zeros, ])
            rm(bins.counts, bins.zeros, idx.counts, idx.zeros)
        } else {
            
            dat <- dplyr::bind_rows(as.data.frame(S4Vectors::mcols(gi)) %>%dplyr::filter(.data$counts>0) %>% dplyr::sample_frac(size = ssize), as.data.frame(S4Vectors::mcols(gi)) %>%dplyr::filter(.data$counts==0) %>% dplyr::sample_frac(size = ssize))
        }
        gc()
        # fit the model on the sample
        if (model_distribution %in% c("nb","nb_hurdle")) {
            fit <- suppressWarnings(.GLM(data=dat, df=df, knots=knots, bdpts=bdpts, covariates=covariates, distance_type=distance_type, model_distribution=model_distribution, splineknotting=splineknotting))
        } else {
            fit <- suppressWarnings(.GLM(data=dat, df=df, knots=knots, bdpts=bdpts, covariates=covariates, distance_type=distance_type, model_distribution=model_distribution, splineknotting=splineknotting, dispersion_DF=dispersion_DF))
        }
        # add predictions
        S4Vectors::mcols(gi)$mu <- NA
        S4Vectors::mcols(gi)$sdev <- NA
        mu<-S4Vectors::mcols(gi)$mu
        sdev<-S4Vectors::mcols(gi)$sdev
        if (model_distribution == "nb_vardisp"){
            sizes<-numeric(length(D.eligible))
        }
        if (model_distribution == "nb_hurdle"){
            phat_count<-numeric(length(D.eligible))
            phat<-numeric(length(D.eligible))
        }
        chunkstep<-100e6
        for (chunkstart in seq(1,length(D.eligible),chunkstep)){
            chunkend<-min(chunkstart+chunkstep,length(D.eligible))
            D.eligible.chunk<-rep(FALSE,length(D.eligible))
            D.eligible.chunk[chunkstart:chunkend]<-D.eligible[chunkstart:chunkend]
            if (!any(D.eligible.chunk)) next
        if (model_distribution == "nb") {
            mu[D.eligible.chunk] <- suppressWarnings(stats::predict(fit, newdata = S4Vectors::mcols(gi)[D.eligible.chunk, ], dispersion = fit$theta^(-1), type = "response"))
            sdev[D.eligible.chunk] <- sqrt(mu[D.eligible.chunk] + mu[D.eligible.chunk]^2/fit$theta)
        } else if (model_distribution == "nb_vardisp") {
            logmu.coefs <- fit@coef[grep("logmu+", names(fit@coef))]
            if (distance_type == "spline") {
                if (!is.null(covariates)) {
                if (splineknotting == "uniform"){
                mu[D.eligible.chunk] <- exp(logmu.coefs[1] + cbind(vapply(covariates, function(x) S4Vectors::mcols(gi)[D.eligible.chunk, x], 
                    rep(8, nrow(S4Vectors::mcols(gi)[D.eligible.chunk, ]))), splines::bs(S4Vectors::mcols(gi)$D[D.eligible.chunk], 
                    df = df, Boundary.knots = bdpts)) %*% logmu.coefs[2:length(logmu.coefs)])
                } else {
                    mu[D.eligible.chunk] <- exp(logmu.coefs[1] + cbind(vapply(covariates, function(x) S4Vectors::mcols(gi)[D.eligible.chunk, x], 
                    rep(8, nrow(S4Vectors::mcols(gi)[D.eligible.chunk, ]))), splines::bs(S4Vectors::mcols(gi)$D[D.eligible.chunk], 
                    knots = knots, Boundary.knots = bdpts)) %*% logmu.coefs[2:length(logmu.coefs)])                    
                }
                } else {
                    if (splineknotting == "uniform"){
                mu[D.eligible.chunk] <- exp(logmu.coefs[1] + cbind(splines::bs(S4Vectors::mcols(gi)$D[D.eligible.chunk], df = df, Boundary.knots = bdpts)) %*% 
                    logmu.coefs[2:length(logmu.coefs)])
                    }else{
                mu[D.eligible.chunk] <- exp(logmu.coefs[1] + cbind(splines::bs(S4Vectors::mcols(gi)$D[D.eligible.chunk], knots = knots, Boundary.knots = bdpts)) %*% 
                    logmu.coefs[2:length(logmu.coefs)])                        
                    }
                }
            } else {
                if (!is.null(covariates)) {
                mu[D.eligible.chunk] <- exp(logmu.coefs[1] + cbind(vapply(covariates, function(x) S4Vectors::mcols(gi)[D.eligible.chunk, x], 
                    rep(8, nrow(S4Vectors::mcols(gi)[D.eligible.chunk,]))), S4Vectors::mcols(gi)$logD[D.eligible.chunk]) %*% logmu.coefs[2:length(logmu.coefs)])
                } else {
                mu[D.eligible.chunk] <- exp(logmu.coefs[1] + cbind(S4Vectors::mcols(gi)$logD[D.eligible.chunk]) %*% logmu.coefs[2:length(logmu.coefs)])
                }
            }
            sizes[D.eligible.chunk] <- 1/alpha(S4Vectors::mcols(gi)$D[D.eligible.chunk])
            sdev[D.eligible.chunk] <- sqrt(mu[D.eligible.chunk] + mu[D.eligible.chunk]^2/sizes)
        } else {
            mu[D.eligible.chunk] <- suppressWarnings(stats::predict(fit, newdata = S4Vectors::mcols(gi)[D.eligible.chunk, ], dispersion = fit$theta^(-1), type = "count"))
            sdev[D.eligible.chunk] <- sqrt(mu[D.eligible.chunk] + mu[D.eligible.chunk]^2/fit$theta)
            phat_count[D.eligible.chunk] <- stats::dnbinom(x = 0, size = fit$theta, mu = mu[D.eligible.chunk])
            phat[D.eligible.chunk] <- 1 - (1 - phat_count[D.eligible.chunk]) * suppressWarnings(stats::predict(fit, newdata = S4Vectors::mcols(gi)[D.eligible.chunk, ], dispersion = fit$theta^(-1), 
                type = "zero"))
        }
        }
        rm(D.eligible.chunk)
        S4Vectors::mcols(gi)$mu[D.eligible] <- mu[D.eligible]
        S4Vectors::mcols(gi)$sdev[D.eligible] <- sdev[D.eligible]
        S4Vectors::mcols(gi)$pvalue <- NA
        S4Vectors::mcols(gi)$qvalue <- NA
        D.fill<-D.eligible&S4Vectors::mcols(gi)$counts>0
        pvalues<-rep(1,length(D.eligible))
        if (any(D.fill)){
        if (model_distribution == "nb") {
            pvalues[D.fill] <- stats::pnbinom(q = S4Vectors::mcols(gi)$counts[D.fill] - 1, size = fit$theta, mu = mu[D.fill], lower.tail = FALSE)
        } else if (model_distribution == "nb_vardisp") {
            pvalues[D.fill] <- stats::pnbinom(q = S4Vectors::mcols(gi)$counts[D.fill] - 1, size = sizes[D.fill], mu = mu[D.fill], lower.tail = FALSE)
        } else {
            pvalues[D.fill] <- ifelse(S4Vectors::mcols(gi)$counts[D.fill] == 0, 1, (1 - phat[D.fill])/(1 - phat_count[D.fill]) * (stats::pnbinom(q = S4Vectors::mcols(gi)$counts[D.fill], 
                size = fit$theta, mu = mu[D.fill], lower.tail = FALSE) + stats::dnbinom(x = S4Vectors::mcols(gi)$counts[D.fill], 
                size = fit$theta, mu = mu[D.fill])))
        }
        }
        qvalues <- stats::p.adjust(pvalues[D.eligible], method = "fdr")
        S4Vectors::mcols(gi)$pvalue[D.eligible] <- pvalues[D.eligible]
        S4Vectors::mcols(gi)$qvalue[D.eligible] <- qvalues
        chrom=as.character(GenomicRanges::seqnames(gi@regions[1,]))
        if (!is.null(model_filepath)) {
            message("Exporting fit object to file")
            fitpath <- path.expand(paste0(gsub("\\.rds$", "", model_filepath), "_Model-Fit_on_", chrom, ".rds"))
            fitpathdir<-gsub("/[^/]+$", "",fitpath)
            if (fitpathdir==fitpath){
                fitpathdir<-gsub("\\[^\\]+$", "",fitpath)
            }
            if (fitpathdir==fitpath){
                fitpathdir<-gsub("\\\\[^\\\\]+$", "",fitpath)
            }
            if (!fitpathdir==fitpath&!dir.exists(fitpathdir)){
                dir.create(fitpathdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            }
            saveRDS(fit, file = fitpath)
        }
    msg<-paste0("Chromosome ",chrom," complete.")
    message(msg)
    return(gi)
}
