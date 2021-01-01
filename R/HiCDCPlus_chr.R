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
#' defaults to 0
#'@param Dmax maximum distance (included) to check for significant interactions,
#' defaults to 2e6 or maximum in the data; whichever is minimum.
#'@param ssize Distance stratified sampling size. Can decrease for
#' large chromosomes. Increase recommended if
#' model fails to converge. Defaults to 0.01.
#'@return A valid \code{gi_list} instance with additional \code{mcols(.)} for
#'each chromosome: pvalue': significance \emph{P}-value, 'qvalue': FDR 
#'corrected \emph{P}-value, mu': expected counts, 'sdev': modeled standard
#'deviation of expected counts.
#'@examples gi_list<-generate_binned_gi_list(50e3,chrs='chr22')
#'gi_list<-add_hic_counts(gi_list,
#'hic_path<-system.file("extdata", "GSE63525_HMEC_combined_example.hic",
#' package = "HiCDCPlus"))
#'gi<-HiCDCPlus_chr(gi_list[[1]])
#'@export

HiCDCPlus_chr <- function(gi, covariates = NULL, distance_type = "spline", model_distribution = "nb", binned = TRUE, 
    df = 6, Dmin = 0, Dmax = 2e+06, ssize = 0.01) {
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
    
    
    # validations done load helper functions outlier removal functions
    remove_outliers_hurdle <- function(dat, mod) {
        mu <- suppressWarnings(stats::predict(mod, newdata = dat, dispersion = mod$theta^(-1), type = "count"))
        dat <- dat %>% dplyr::mutate(phat_count = stats::dnbinom(x = 0, size = mod$theta, mu = mu))
        dat <- dat %>% dplyr::mutate(phat = 1 - (1 - .data$phat_count) * suppressWarnings(stats::predict(mod, newdata = dat, dispersion = mod$theta^(-1), 
            type = "zero")))
        dat <- dat %>% dplyr::mutate(pvals = ifelse(.data$counts == 0, 1, (1 - .data$phat)/(1 - .data$phat_count) * (stats::pnbinom(q = .data$counts, 
            size = mod$theta, mu = mu, lower.tail = FALSE) + stats::dnbinom(x = .data$counts, size = mod$theta, mu = mu))))
        new.dat <- dat %>% dplyr::filter(.data$pvals >= 0.025) %>% dplyr::select(-.data$phat_count, -.data$phat, -.data$pvals)
        return(new.dat)
    }
    
    remove_outliers_nb <- function(dat, mod) {
        mu <- suppressWarnings(stats::predict(mod, newdata = dat, type = "response"))
        dat <- dat %>% dplyr::mutate(q = stats::qnbinom(0.975, size = mod$theta, mu = mu))
        new.dat <- dat %>% dplyr::filter(.data$counts <= .data$q) %>% dplyr::select(-.data$q)
        return(new.dat)
    }
    
    remove_outliers_nb_vardisp <- function(dat, covariates, dispersion_DF, mod) {
        dispersionfunction <- function(DF) {
            fit <- stats::smooth.spline(DF$D, 1/(DF$alpha))
            return(function(x) base::pmax(suppressWarnings(stats::predict(fit, x)$y), 0.001))
        }
        dispersion <- dispersionfunction(dispersion_DF)
        sizes <- 1/dispersion(dat$D)
        logmu.coefs <- mod@coef[grep("logmu+", names(mod@coef))]
        if (distance_type == "spline") {
            if (!is.null(covariates)) {
                mus <- exp(logmu.coefs[1] + cbind(vapply(covariates, function(x) dat[, x], rep(8, nrow(dat))), splines::bs(dat$D, 
                df = df, Boundary.knots = bdpts)) %*% logmu.coefs[2:length(logmu.coefs)])
            } else {
                mus <- exp(logmu.coefs[1] + cbind(splines::bs(dat$D, df = df, Boundary.knots = bdpts)) %*% logmu.coefs[2:length(logmu.coefs)])
            }
        } else {
            if (!is.null(covariates)) {
                mus <- exp(logmu.coefs[1] + cbind(vapply(covariates, function(x) dat[, x], rep(8, nrow(dat))), dat$logD) %*% 
                logmu.coefs[2:length(logmu.coefs)])
            } else {
                mus <- exp(logmu.coefs[1] + cbind(dat$logD) %*% logmu.coefs[2:length(logmu.coefs)])
            }
        }
        q <- stats::qnbinom(0.975, size = sizes, mu = mus)
        counts <- dat$counts
        counts[counts > q] <- NA
        new.dat <- dat %>% dplyr::filter(!is.na(counts))
        return(new.dat)
    }
    
    # glm.nb wrapper with poisson fallback if dispersion undetectable
    glm.nb.trycatch <- function(model.formula, data) {
        model <- tryCatch({
            strip_glm(MASS::glm.nb(model.formula, data, model = FALSE, y = FALSE, x = FALSE))
        }, error = function(e) {
            temp.model <- strip_glm(stats::glm(model.formula, data, family = "poisson", model = FALSE, y = FALSE, x = FALSE))
            temp.model$theta = 1e+16
            return(temp.model)
        })
        return(model)
    }
    
    GLM_nb <- function(data, df, bdpts, covariates, distance_type) {
        if (distance_type == "spline") {
            model.formula <- stats::as.formula(paste0("counts ~", paste(covariates, collapse = " + "), " + splines::bs(D,df=", 
                df, ",Boundary.knots = c(", paste(bdpts, collapse = ","), "))"))
        } else {
            model.formula <- stats::as.formula(paste0("counts ~", paste(covariates, collapse = " + "), " + logD"))
        }
        # Remove outliers
        mod <- glm.nb.trycatch(model.formula, data)
        new.dat <- remove_outliers_nb(data, mod)
        # Refit the model
        mod <- glm.nb.trycatch(model.formula, new.dat)
        return(mod)
    }
    
    
    GLM_nb_vardisp <- function(data, df, bdpts, covariates, dispersion_DF, distance_type) {
        dispersionfunction <- function(DF) {
            fit <- stats::smooth.spline(DF$D, 1/(DF$alpha))
            return(function(x) base::pmax(suppressWarnings(stats::predict(fit, x)$y), 0.001))
        }
        dispersion <- dispersionfunction(dispersion_DF)
        data$sizes <- 1/dispersion(data$D)
        if (distance_type == "spline") {
            if (!is.null(covariates)) {
                model.formula <- stats::as.formula(paste0("counts ~", paste(covariates, collapse = " + "), " + splines::bs(D,df=", 
                df, ",Boundary.knots = c(", paste(bdpts, collapse = ","), "))"))
            } else {
                model.formula <- stats::as.formula(paste0("counts ~ splines::bs(D,df=", df, ",Boundary.knots = c(", paste(bdpts, 
                collapse = ","), "))"))
            }
            mod <- glm.nb.trycatch(model.formula, data)
            start <- list(logmu = mod$coefficients)
            mod <- bbmle::mle2(counts ~ dnbinom(mu = exp(logmu), size = sizes), data = data, start = start, parameters = list(stats::as.formula(paste0("logmu ~ ", 
                paste(covariates, collapse = " + "), " + splines::bs(D,df=", df, ",Boundary.knots = c(", paste(bdpts, collapse = ","), 
                "))"))), method = "BFGS")
        } else {
            if (!is.null(covariates)) {
                model.formula <- stats::as.formula(paste0("counts ~", paste(covariates, collapse = " + "), " + logD"))
            } else {
                model.formula <- stats::as.formula(paste0("counts ~ logD"))
            }
            mod <- glm.nb.trycatch(model.formula, data)
            start <- list(logmu = mod$coefficients)
            if (!is.null(covariates)) {
                mod <- bbmle::mle2(counts ~ dnbinom(mu = exp(logmu), size = sizes), data = data, start = start, parameters = list(stats::as.formula(paste0("logmu ~", 
                paste(covariates, collapse = " + "), " + logD"))), method = "BFGS")
            } else {
                mod <- bbmle::mle2(counts ~ dnbinom(mu = exp(logmu), size = sizes), data = data, start = start, parameters = list(stats::as.formula(paste0("logmu ~ logD"))), 
                method = "BFGS")
            }
        }
        # Remove outliers
        new.dat <- remove_outliers_nb_vardisp(data, covariates, dispersion_DF, mod)
        logmu.coefs <- mod@coef[grep("logmu+", names(mod@coef))]
        start <- list(logmu = logmu.coefs)
        new.dat$sizes <- 1/dispersion(new.dat$D)
        # Refit the model
        if (distance_type == "spline") {
            if (!is.null(covariates)) {
                mod <- bbmle::mle2(counts ~ dnbinom(mu = exp(logmu), size = sizes), data = new.dat, start = start, parameters = list(stats::as.formula(paste0("logmu ~ ", 
                paste(covariates, collapse = " + "), " + splines::bs(D,df=", df, ",Boundary.knots = c(", paste(bdpts, collapse = ","), 
                "))"))), method = "BFGS")
            } else {
                mod <- bbmle::mle2(counts ~ dnbinom(mu = exp(logmu), size = sizes), data = new.dat, start = start, parameters = list(stats::as.formula(paste0("logmu ~  splines::bs(D,df=", 
                df, ",Boundary.knots = c(", paste(bdpts, collapse = ","), "))"))), method = "BFGS")
            }
        } else {
            if (!is.null(covariates)) {
                mod <- bbmle::mle2(counts ~ dnbinom(mu = exp(logmu), size = sizes), data = new.dat, start = start, parameters = list(stats::as.formula(paste0("logmu ~", 
                paste(covariates, collapse = " + "), " + logD"))), method = "BFGS")
            } else {
                mod <- bbmle::mle2(counts ~ dnbinom(mu = exp(logmu), size = sizes), data = new.dat, start = start, parameters = list(stats::as.formula(paste0("logmu ~ logD"))), 
                method = "BFGS")
            }
        }
        return(mod)
    }
    
    GLM <- function(data, df, bdpts, covariates, distance_type) {
        if (distance_type == "spline") {
            if (!is.null(covariates)) {
                mod <- strip_glm(pscl::hurdle(stats::as.formula(paste0("counts ~", paste(covariates, collapse = " + "), " + splines::bs(D,df=", 
                df, ",Boundary.knots = c(", paste(bdpts, collapse = ","), ")) | ", paste(covariates, collapse = " + "), 
                " + splines::bs(D,df=", df, ",Boundary.knots =c(", paste(bdpts, collapse = ","), "))")), data = data, dist = "negbin", 
                model = FALSE, y = FALSE, x = FALSE))
            } else {
                mod <- strip_glm(pscl::hurdle(stats::as.formula(paste0("counts ~ splines::bs(D,df=", df, ",Boundary.knots = c(", 
                paste(bdpts, collapse = ","), ")) | splines::bs(D,df=", df, ",Boundary.knots =c(", paste(bdpts, collapse = ","), 
                "))")), data = data, dist = "negbin", model = FALSE, y = FALSE, x = FALSE))
            }
        } else {
            if (!is.null(covariates)) {
                mod <- strip_glm(pscl::hurdle(stats::as.formula(paste0("counts ~", paste(covariates, collapse = " + "), " + logD | ", 
                paste(covariates, collapse = " + "), " + logD")), data = data, dist = "negbin", model = FALSE, y = FALSE, x = FALSE))
            } else {
                mod <- strip_glm(pscl::hurdle(stats::as.formula(paste0("counts ~ logD | logD")), data = data, dist = "negbin", 
                model = FALSE, y = FALSE, x = FALSE))
            }
        }
        # Remove outliers
        new.dat <- remove_outliers_hurdle(data, mod)
        # Refit the model
        if (distance_type == "spline") {
            if (!is.null(covariates)) {
                mod <- strip_glm(pscl::hurdle(stats::as.formula(paste0("counts ~", paste(covariates, collapse = " + "), " + splines::bs(D,df=", 
                df, ",Boundary.knots = c(", paste(bdpts, collapse = ","), ")) | ", paste(covariates, collapse = " + "), 
                " + splines::bs(D,df=", df, ",Boundary.knots =c(", paste(bdpts, collapse = ","), "))")), data = new.dat, 
                dist = "negbin", model = FALSE, y = FALSE, x = FALSE))
            } else {
                mod <- strip_glm(pscl::hurdle(stats::as.formula(paste0("counts ~ splines::bs(D,df=", df, ",Boundary.knots = c(", 
                paste(bdpts, collapse = ","), ")) |  splines::bs(D,df=", df, ",Boundary.knots =c(", paste(bdpts, collapse = ","), 
                "))")), data = new.dat, dist = "negbin", model = FALSE, y = FALSE, x = FALSE))
            }
        } else {
            if (!is.null(covariates)) {
                mod <- strip_glm(pscl::hurdle(stats::as.formula(paste0("counts ~", paste(covariates, collapse = " + "), " + logD | ", 
                paste(covariates, collapse = " + "), " + logD")), data = new.dat, dist = "negbin", model = FALSE, y = FALSE, 
                x = FALSE))
            } else {
                mod <- strip_glm(pscl::hurdle(stats::as.formula(paste0("counts ~ logD | logD")), data = new.dat, dist = "negbin", 
                model = FALSE, y = FALSE, x = FALSE))
            }
        }
        return(mod)
    }
    
    strip_glm = function(m1) {
        m1$data <- NULL
        m1$y <- NULL
        m1$linear.predictors <- NULL
        m1$weights <- NULL
        m1$fitted.values <- NULL
        m1$model <- NULL
        m1$prior.weights <- NULL
        m1$residuals <- NULL
        m1$effects <- NULL
        m1$qr$qr <- NULL
        m1$optim <- NULL
        m1$start <- NULL
        return(m1)
    }
    
    gm_mean = function(x, na.rm = TRUE) {
        exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
    }
    
    dispersionfunction <- function(DF) {
        fit <- stats::smooth.spline(DF$D, 1/(DF$alpha))
        return(function(x) base::pmax(suppressWarnings(stats::predict(fit, x)$y), 0.001))
    }
    # main routine for the function

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
            for (x in unique(dat$D.range)) {
                zeroPairedBins_ix <- seq(nrow(dat))[dat$D.range == x & dat$counts == 0]
                countPairedBins_ix <- seq(nrow(dat))[dat$D.range == x & dat$counts != 0]
                dat_x_ix <- c(countPairedBins_ix[sample(length(countPairedBins_ix), floor(length(countPairedBins_ix) * ssize), 
                replace = FALSE)], zeroPairedBins_ix[sample(length(zeroPairedBins_ix), floor(length(zeroPairedBins_ix) * ssize), replace = FALSE)])
                bdpts_x <- range(S4Vectors::mcols(gi)[dat_x_ix,]$D)
                fit <- suppressWarnings(GLM_nb(as.data.frame(S4Vectors::mcols(gi)[dat_x_ix,]), df, bdpts_x, covariates, distance_type = distance_type))
                dt <- data.frame(D = mean(S4Vectors::mcols(gi)[dat_x_ix,]$D), alpha = fit$theta)
                dispersion_DF <- dplyr::bind_rows(dispersion_DF, dt)
            }
            dispersion <- dispersionfunction(dispersion_DF)
            rm(dat, dat_x_ix, bdpts_x, fit, dt, zeroPairedBins_ix, countPairedBins_ix)
        }
        # get a stratified sample for modeling
        dat <- as.data.frame(S4Vectors::mcols(gi))%>%dplyr::select(.data$D,.data$counts)
        if (binned) {
            get.range <- function(gi) {
                return(stats::quantile(GenomicRanges::end(InteractionSet::regions(gi)) - GenomicRanges::start(InteractionSet::regions(gi)), 
                                       probs = c(0.025, 0.975)))
            }
            binsize <- min(unique(get.range(gi)))
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
        if (model_distribution == "nb") {
            fit <- suppressWarnings(GLM_nb(dat, df, bdpts, covariates, distance_type))
        } else if (model_distribution == "nb_vardisp") {
            fit <- suppressWarnings(GLM_nb_vardisp(dat, df, bdpts, covariates, dispersion_DF, distance_type))
        } else {
            fit <- suppressWarnings(GLM(dat, df, bdpts, covariates, distance_type))
        }
        # add predictions
        S4Vectors::mcols(gi)$mu <- NA
        S4Vectors::mcols(gi)$sdev <- NA
        if (model_distribution == "nb") {
            mu <- suppressWarnings(stats::predict(fit, newdata = S4Vectors::mcols(gi)[D.eligible, ], dispersion = fit$theta^(-1), type = "response"))
            sdev <- sqrt(mu + mu^2/fit$theta)
        } else if (model_distribution == "nb_vardisp") {
            logmu.coefs <- fit@coef[grep("logmu+", names(fit@coef))]
            if (distance_type == "spline") {
                if (!is.null(covariates)) {
                mu <- exp(logmu.coefs[1] + cbind(vapply(covariates, function(x) S4Vectors::mcols(gi)[D.eligible, x], 
                    rep(8, nrow(S4Vectors::mcols(gi)[D.eligible, ]))), splines::bs(S4Vectors::mcols(gi)$D[D.eligible], 
                    df = df, Boundary.knots = bdpts)) %*% logmu.coefs[2:length(logmu.coefs)])
                } else {
                mu <- exp(logmu.coefs[1] + cbind(splines::bs(S4Vectors::mcols(gi)$D[D.eligible], df = df, Boundary.knots = bdpts)) %*% 
                    logmu.coefs[2:length(logmu.coefs)])
                }
            } else {
                if (!is.null(covariates)) {
                mu <- exp(logmu.coefs[1] + cbind(vapply(covariates, function(x) S4Vectors::mcols(gi)[D.eligible, x], 
                    rep(8, nrow(S4Vectors::mcols(gi)[D.eligible, ]))), S4Vectors::mcols(gi)$logD[D.eligible]) %*% logmu.coefs[2:length(logmu.coefs)])
                } else {
                mu <- exp(logmu.coefs[1] + cbind(S4Vectors::mcols(gi)$logD[D.eligible]) %*% logmu.coefs[2:length(logmu.coefs)])
                }
            }
            sizes <- 1/dispersion(S4Vectors::mcols(gi)$D[D.eligible])
            sdev <- sqrt(mu + mu^2/sizes)
        } else {
            mu <- suppressWarnings(stats::predict(fit, newdata = S4Vectors::mcols(gi)[D.eligible, ], dispersion = fit$theta^(-1), type = "count"))
            sdev <- sqrt(mu + mu^2/fit$theta)
            phat_count <- stats::dnbinom(x = 0, size = fit$theta, mu = mu)
            phat <- 1 - (1 - phat_count) * suppressWarnings(stats::predict(fit, newdata = S4Vectors::mcols(gi)[D.eligible, ], dispersion = fit$theta^(-1), 
                type = "zero"))
        }
        S4Vectors::mcols(gi)$mu[D.eligible] <- mu
        S4Vectors::mcols(gi)$sdev[D.eligible] <- sdev
        S4Vectors::mcols(gi)$pvalue <- NA
        S4Vectors::mcols(gi)$qvalue <- NA
        if (model_distribution == "nb") {
            pvalues <- stats::pnbinom(q = S4Vectors::mcols(gi)$counts[D.eligible] - 1, size = fit$theta, mu = mu, lower.tail = FALSE)
        } else if (model_distribution == "nb_vardisp") {
            pvalues <- stats::pnbinom(q = S4Vectors::mcols(gi)$counts[D.eligible] - 1, size = sizes, mu = mu, lower.tail = FALSE)
        } else {
            pvalues <- ifelse(S4Vectors::mcols(gi)$counts[D.eligible] == 0, 1, (1 - phat)/(1 - phat_count) * (stats::pnbinom(q = S4Vectors::mcols(gi)$counts[D.eligible], 
                size = fit$theta, mu = mu, lower.tail = FALSE) + stats::dnbinom(x = S4Vectors::mcols(gi)$counts[D.eligible], 
                size = fit$theta, mu = mu)))
        }
        qvalues <- stats::p.adjust(pvalues, method = "fdr")
        S4Vectors::mcols(gi)$pvalue[D.eligible] <- pvalues
        S4Vectors::mcols(gi)$qvalue[D.eligible] <- qvalues
    print(paste0("Chromosome ",as.character(GenomicRanges::seqnames(gi@regions[1,]))," complete."))
    return(gi)
}
