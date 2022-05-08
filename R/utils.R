#' .get_range
#'
#' Infer binsize of a gi instance
#' @param gi a valid gi_list element
#' @noRd
.get_range <- function(gi) {
    return(min(unique(stats::quantile(GenomicRanges::end(InteractionSet::regions(gi)) - GenomicRanges::start(InteractionSet::regions(gi)), 
        probs = c(0.025, 0.975)))))
}


#' .strip_glm
#'
#' Make GLMs occupy less space by removing memory consuming elements
#' That would not impede with predict
#' @param mod a stats::glm model object
#' @noRd
    .strip_glm = function(mod) {
        mod$data <- NULL
        mod$y <- NULL
        mod$linear.predictors <- NULL
        mod$weights <- NULL
        mod$fitted.values <- NULL
        mod$model <- NULL
        mod$prior.weights <- NULL
        mod$residuals <- NULL
        mod$effects <- NULL
        mod$qr$qr <- NULL
        mod$optim <- NULL
        mod$start <- NULL
        return(mod)
    }

#' .dispersionfunction
#'  
#' Third degree polynomial fits on an empirical distance versus dispersion relationship
#' @param DF dataframe containing the median distance and dispersion associated with each distance range bin
#' @noRd
    .dispersionfunction <- function(DF) {
        fit <- stats::lm(1/(theta)~stats::poly(D,degree = 3),data=DF)
        return(function(x) base::pmax(suppressWarnings(stats::predict(fit, newdata=data.frame(D=x))), 0.001))
    }

#' .weighted_quantile
#'
#' Generate quantiles from a vector of values and associated weights. Values should not be duplicated.
#' @param x vector of values, should not have duplicate elements
#' @param w vector of weights, cannot be negative
#' @param num number of equally spaced quantiles to be generated
#' @noRd
    .weighted_quantile <- function(x, w, num) {
        ord <- order(x)
        w <- w[ord]
        x <- x[ord]
        w.ord <- cumsum(w) / sum(w)
        index <- seq(1,length(x),1)
        quantiles<-numeric(num)
        for (ix in seq(1,num,1)){
            prob<-ix/(num+1)
            if(min(w.ord) > prob) {
                lower.k.quant <- 1
            } else {
                lower.k.quant <- max(index[w.ord <= prob])
            }
            upper.k.quant <- min(index[w.ord > prob])
            if(w.ord[lower.k.quant] < prob) {
                quantiles[ix]<-x[upper.k.quant]
                next
            } else {
                quantiles[ix]<-(w[lower.k.quant] * x[lower.k.quant] +
                                    w[upper.k.quant] * x[upper.k.quant]) /
                    (w[lower.k.quant] + w[upper.k.quant])
                next
            }
        }
        quantiles<-quantiles[!duplicated(quantiles)]
        return(quantiles)
    }

#' .remove_outliers
#'
#' Remove observations with low probability of fitting the background model to increase detection power
#' @param dat data frame to be filtered
#' @param mod model to assess likelihood of fitting background model
#' @param model_distribution 'nb' uses a Negative Binomial model, 'nb_vardisp' uses a Negative Binomial model with a distance specific dispersion parameter inferred from the data, 'nb_hurdle' uses the legacy HiCDC model.
#' @param ... additional arguments for the specific model distribution to ensure prediction
#' @noRd
    .remove_outliers <- function(dat, mod, model_distribution,...) {
        if (model_distribution=="nb_hurdle"){
        mu <- suppressWarnings(stats::predict(mod, newdata = dat, dispersion = mod$theta^(-1), type = "count"))
        dat <- dat %>% dplyr::mutate(phat_count = stats::dnbinom(x = 0, size = mod$theta, mu = mu))
        dat <- dat %>% dplyr::mutate(phat = 1 - (1 - .data$phat_count) * suppressWarnings(stats::predict(mod, newdata = dat, dispersion = mod$theta^(-1), 
            type = "zero")))
        dat <- dat %>% dplyr::mutate(pvals = ifelse(.data$counts == 0, 1, (1 - .data$phat)/(1 - .data$phat_count) * (stats::pnbinom(q = .data$counts, 
            size = mod$theta, mu = mu, lower.tail = FALSE) + stats::dnbinom(x = .data$counts, size = mod$theta, mu = mu))))
        new.dat <- dat %>% dplyr::filter(.data$pvals >= 0.025) %>% dplyr::select(-.data$phat_count, -.data$phat, -.data$pvals)
        return(new.dat)
        }
        if (model_distribution=="nb"){
            mu <- suppressWarnings(stats::predict(mod, newdata = dat, type = "response"))
            dat <- dat %>% dplyr::mutate(q = stats::qnbinom(0.975, size = mod$theta, mu = mu))
            new.dat <- dat %>% dplyr::filter(.data$counts <= .data$q) %>% dplyr::select(-.data$q)
            return(new.dat)            
        }
        if (model_distribution=="nb_vardisp"){
            args<-list(...)
            alpha <- .dispersionfunction(args$dispersion_DF)
            sizes <- 1/alpha(dat$D)
            logmu.coefs <- mod@coef[grep("logmu+", names(mod@coef))]
            if (args$distance_type == "spline") {
                if (!is.null(args$covariates)) {
                    if (args$splineknotting == "uniform"){
                        mus <- exp(logmu.coefs[1] + cbind(vapply(args$covariates, function(x) dat[, x], rep(8, nrow(dat))), splines::bs(dat$D, 
                                                                                                                                   df = args$df, Boundary.knots = args$bdpts)) %*% logmu.coefs[2:length(logmu.coefs)])}
                    else{
                        mus <- exp(logmu.coefs[1] + cbind(vapply(args$covariates, function(x) dat[, x], rep(8, nrow(dat))), splines::bs(dat$D, 
                                                                                                                                   knots = args$knots, Boundary.knots = args$bdpts)) %*% logmu.coefs[2:length(logmu.coefs)])                   
                    }
                } else {
                    if (args$splineknotting == "uniform"){
                        mus <- exp(logmu.coefs[1] + cbind(splines::bs(dat$D, df = args$df, Boundary.knots = args$bdpts)) %*% logmu.coefs[2:length(logmu.coefs)])
                    }else{
                        mus <- exp(logmu.coefs[1] + cbind(splines::bs(dat$D, knots = args$knots, Boundary.knots = args$bdpts)) %*% logmu.coefs[2:length(logmu.coefs)])
                    }
                }
            } else {
                if (!is.null(args$covariates)) {
                    mus <- exp(logmu.coefs[1] + cbind(vapply(args$covariates, function(x) dat[, x], rep(8, nrow(dat))), dat$logD) %*% 
                                   logmu.coefs[2:length(logmu.coefs)])
                } else {
                    mus <- exp(logmu.coefs[1] + cbind(dat$logD) %*% logmu.coefs[2:length(logmu.coefs)])
                }
            }
            dat <- dat %>% dplyr::mutate(q = stats::qnbinom(0.975, size = sizes, mu = mus))
            new.dat <- dat %>% dplyr::filter(.data$counts <= .data$q) %>% dplyr::select(-.data$q)
            return(new.dat)
        }            
    }

#' .glm_nb_trycatch
#'
#' Wraps a MASS::glm.nb model with a poisson fallback (stats::glm(...,family="poisson")) in case
#' @param model_formula data frame to be filtered
#' @param data data used for modeling
#' @noRd
    .glm_nb_trycatch <- function(model_formula, data) {
        model <- tryCatch({
            .strip_glm(MASS::glm.nb(model_formula, data, model = FALSE, y = FALSE, x = FALSE))
        }, error = function(e) {
            temp_model <- .strip_glm(stats::glm(model_formula, data, family = "poisson", model = FALSE, y = FALSE, x = FALSE))
            temp_model$theta = 1e+16
            return(temp_model)
        })
        return(model)
    }


#'.GLM
#'
#' GLM models for the underlying distributions following the description in citation(package = "HiCDCPlus")
#' @param data data used for modeling
#' @param df degrees of freedom for the genomic distance spline function if \code{distance_type='spline'}. Defaults to 6, which corresponds to a cubic spline as explained in Carty et al. (2017)
#' @param knots knotting points for splines, if different from default
#' @param bdpts boundary points for the spline
#' @param covariates covariates to be considered in addition to genomic distance D. Defaults to all covariates besides  'D','counts','mu','sdev',pvalue','qvalue' 
#' @param distance_type distance covariate form: 'spline' or 'log'. Defaults to 'spline'.
#' @param model_distribution 'nb' uses a Negative Binomial model, 'nb_vardisp' uses a Negative Binomial model with a distance specific dispersion parameter inferred from the data, 'nb_hurdle' uses the legacy HiCDC model.
#' @param splineknotting pline knotting strategy. Either "uniform", uniformly spaced in distance, or placed based on distance distribution of counts "count-based" (i.e., more closely spaced where counts are more dense).
#' @param ... additional arguments for the specific model distribution to ensure prediction
#' @noRd
    .GLM<-function(data, df, knots, bdpts, covariates, distance_type, model_distribution, splineknotting, ...){
        if (distance_type == "spline") {
            if (splineknotting == "uniform"){
                model.formula.RHS<-paste0(paste(covariates, collapse = " + "), ifelse(!is.null(covariates)," + ",""), "splines::bs(D,df=", 
                                          df, ",Boundary.knots = c(", paste(bdpts, collapse = ","), "))")
            }else{
                model.formula.RHS<-paste0(paste(covariates, collapse = " + "), ifelse(!is.null(covariates)," + ",""), "splines::bs(D,knots=c(", 
                                          paste(knots, collapse = ","), "), Boundary.knots = c(", paste(bdpts, collapse = ","), "))")
            }
        } else {
                model.formula.RHS <- paste0(paste(covariates, collapse = " + "), ifelse(!is.null(covariates)," + ",""),"logD")
        }
        # Run model
        if (model_distribution %in% c("nb_vardisp","nb")){
        model.formula <- stats::as.formula(paste0("counts ~", model.formula.RHS))
        mod <- .glm_nb_trycatch(model.formula, data)
        }
        #Retrain on variable dispersion after warm start with fixed dispersion
        if (model_distribution=="nb_vardisp"){

            varg<-list(...)
            alpha <- .dispersionfunction(varg$dispersion_DF)
            data$sizes <- 1/alpha(data$D)    
            start <- list(logmu = mod$coefficients)
            mod <- bbmle::mle2(counts ~ dnbinom(mu = exp(logmu), size = sizes), data = data, start = start, parameters = list(stats::as.formula(paste0("logmu ~ ", model.formula.RHS))), method = "BFGS")
        } 
        if (model_distribution=="nb_hurdle") {
        model.formula <- stats::as.formula(paste0("counts ~", model.formula.RHS," | ",model.formula.RHS))
        mod <- .strip_glm(pscl::hurdle(model.formula, data = data, dist = "negbin", 
                                           model = FALSE, y = FALSE, x = FALSE))    
        }
        #Remove outliers
        if (model_distribution%in%c("nb","nb_hurdle")){
            new.dat <- .remove_outliers(dat=data, mod=mod, model_distribution = model_distribution)
        } else {
            new.dat <- .remove_outliers(dat=data, mod=mod, model_distribution = model_distribution, covariates = covariates, dispersion_DF = varg$dispersion_DF, splineknotting = splineknotting, df= df, bdpts=bdpts, knots=knots, distance_type=distance_type)
        }
        #Refit the model and return
        if (model_distribution=="nb"){
            mod <- .glm_nb_trycatch(model.formula, new.dat)
            return(mod)
        }  
        if (model_distribution=="nb_vardisp"){
            logmu.coefs <- mod@coef[grep("logmu+", names(mod@coef))]
            new.dat$sizes <- 1/alpha(new.dat$D)
            start <- list(logmu = logmu.coefs)
            mod <- bbmle::mle2(counts ~ dnbinom(mu = exp(logmu), size = sizes), data = new.dat, start = start, parameters = list(stats::as.formula(paste0("logmu ~ ", model.formula.RHS))), method = "BFGS")
            return(mod)
        } 
        if (model_distribution=="nb_hurdle"){
            mod <- .strip_glm(pscl::hurdle(model.formula, data = new.dat, dist = "negbin", 
                                           model = FALSE, y = FALSE, x = FALSE))  
            return(mod)
        }
    }



#'.TopDom
#'
#'Adapted version of the stable legacy TopDom package 
#'version 0.0.2 (no longer on CRAN or Bioconductor) 
#'written by Hanjun Shin(shanjun "at" usc.edu), contributions by 
#'Harris Lazaris(Ph.D Stduent, NYU), Dr. Gangqing Hu(Staff Scientist, NIH). 
#'If you're using this function, please cite 
#'TopDom according to the documentation at 
#'https://github.com/HenrikBengtsson/TopDom/blob/0.0.2/docs/.
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param matrix.file string, file address,
#'Has a structure of N * (N + 3), where N is the number of bins, see Vignette
#'at https://github.com/HenrikBengtsson/TopDom/blob/0.0.2/docs
#'@param window.size number of bins to extend. Defaults to 5
#'@param outFile path prefix for TAD annotations to write named across
#'chromosomes. Defaults to NULL (no file output).
#'@param statFilter whether a Wilcoxon rank sum (unpaired) test based filtering
#'of TAD boundaries would be
#' done based on a 0.05 \emph{P}-value threshold. Defaults to TRUE.
#'@param verbose Choose TRUE if you'd like to troubleshoot TopDom process.
#'@return A list of TAD annotations per each chromosome
#'@noRd
.TopDom <- function(matrix.file, window.size=5, outFile=NULL, statFilter=TRUE, verbose=FALSE)
    {
    #function definitions
    Get.Diamond.Matrix <- function(mat.data, i, size)
    {
        n_bins = nrow( mat.data )
        if(i==n_bins) return(NA)
        
        lowerbound = max( 1, i-size+1 )
        upperbound = min( i+size, n_bins)
        
        return( mat.data[lowerbound:i, (i+1):upperbound] )
    }
    
    
    Which.process.region <- function(rmv.idx, n_bins, min.size=3)
    {
        gap.idx = rmv.idx
        
        proc.regions = data.frame(start=numeric(0), end=numeric(0))
        proc.set = setdiff(seq(1,n_bins,1), gap.idx)
        n_proc.set = length(proc.set)
        
        i=1
        while(i < n_proc.set )
        {
            start = proc.set[i]
            j = i+1
            
            while(j <= n_proc.set)
            {
                if( proc.set[j] - proc.set[j-1] <= 1) j = j + 1
                else {
                    proc.regions = rbind(proc.regions, c(start=start, end=proc.set[j-1]) )
                    i = j
                    break
                }
            }
            
            if(j >= n_proc.set ) {
                proc.regions = rbind(proc.regions, c(start=start, end=proc.set[j-1]) )
                break
            }
        }
        
        colnames(proc.regions) = c("start", "end")
        proc.regions <- proc.regions[ which( abs(proc.regions[,"end"] - proc.regions[, "start"]) >= min.size ), ]
        
        return(proc.regions)
    }
    
    
    Which.Gap.Region <- function(matrix.data, w)
    {
        n_bins = nrow(matrix.data)
        gap = rep(0, n_bins)
        
        for(i in seq(1,n_bins,1))
        {
            if( sum( matrix.data[i, max(1, i-w):min(i+w, n_bins)] ) == 0 ) gap[i]=-0.5
        }
        
        idx = which(gap == -0.5)
        return(idx)
    }
    
    Detect.Local.Extreme <- function(x)
    {
        n_bins = length(x)
        ret = rep(0, n_bins)
        x[is.na(x)]=0
        
        if(n_bins <= 3)
        {
            ret[which.min(x)]=-1
            ret[which.max(x)]=1
            
            return(ret)
        }
        # Norm##################################################3
        new.point = Data.Norm(x=seq(1,n_bins,1), y=x)
        x=new.point$y
        ##################################################
        cp = Change.Point(x=seq(1,n_bins,1), y=x)
        
        if( length(cp$cp) <= 2 ) return(ret)
        if( length(cp$cp) == n_bins) return(ret)
        for(i in 2:(length(cp$cp)-1))
        {
            if( x[cp$cp[i]] >= x[cp$cp[i]-1] && x[cp$cp[i]] >= x[cp$cp[i]+1] ) ret[cp$cp[i]] = 1
            else if(x[cp$cp[i]] < x[cp$cp[i]-1] && x[cp$cp[i]] < x[cp$cp[i]+1]) ret[cp$cp[i]] = -1
            
            min.val = min( x[ cp$cp[i-1] ], x[ cp$cp[i] ] )
            max.val = max( x[ cp$cp[i-1] ], x[ cp$cp[i] ] )
            
            if( min( x[cp$cp[i-1]:cp$cp[i]] ) < min.val ) ret[ cp$cp[i-1] - 1 + which.min( x[cp$cp[i-1]:cp$cp[i]] ) ] = -1
            if( max( x[cp$cp[i-1]:cp$cp[i]] ) > max.val ) ret[ cp$cp[i-1] - 1 + which.max( x[cp$cp[i-1]:cp$cp[i]] ) ] = 1
        }
        
        return(ret)
    }
    

    Data.Norm <- function(x, y)
    {
        ret.x = rep(0, length(x))
        ret.y = rep(0, length(y))
        
        ret.x[1] = x[1]
        ret.y[1] = y[1]
        
        diff.x = diff(x)
        diff.y = diff(y)
        
        scale.x = 1 / mean( abs(diff(x) ) )
        scale.y = 1 / mean( abs( diff(y) ) )
        
        for(i in 2:length(x))
        {
            ret.x[i] = ret.x[i-1] + (diff.x[i-1]*scale.x)
            ret.y[i] = ret.y[i-1] + (diff.y[i-1]*scale.y)
        }
        
        return(list(x=ret.x, y=ret.y))
    }
    
    
    Change.Point <- function( x, y )
    {
        if( length(x) != length(y)) 
        {
            stop("ERROR : The length of x and y should be the same")
        }
        
        n_bins <- length(x)
        Fv <- rep(NA, n_bins)
        Ev <- rep(NA, n_bins)
        cp <- 1
        
        i=1
        Fv[1]=0
        while( i < n_bins )
        {
            j=i+1
            Fv[j] = sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 )
            
            while(j<n_bins)
            {
                j=j+1
                k=(i+1):(j-1)
                Ev[j] = ( sum( abs( (y[j]-y[i] )*x[k] - (x[j] -x[i])*y[k] - (x[i]*y[j]) + (x[j]*y[i]) ) ) / sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 ) )
                Fv[j] = sqrt( (x[j]-x[i])^2 + (y[j] - y[i])^2 ) - ( sum( abs( (y[j]-y[i] )*x[k] - (x[j] -x[i])*y[k] - (x[i]*y[j]) + (x[j]*y[i]) ) ) / sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 ) )
                
                #################################################
                #Not Original Code
                if( is.na(Fv[j]) || is.na(Fv[j-1]) ) {
                    j = j-1
                    cp <- c(cp, j)
                    break
                }
                ####################################################3
                if(Fv[j] < Fv[j-1] ) {
                    j = j - 1
                    cp <- c(cp, j )
                    break
                }
            }
            i=j
        }
        
        cp <- c(cp, n_bins)
        
        return(list(cp=cp, objF=Fv, errF=Ev))
    }
    
    
    Get.Pvalue <- function( matrix.data, size, scale=1 )
    {
        n_bins = nrow(matrix.data)
        pvalue <- rep(1, n_bins)
        
        for( i in seq(1,(n_bins-1),1) )
        {
            dia = as.vector( Get.Diamond.Matrix2(matrix.data, i, size=size) )
            ups = as.vector( Get.Upstream.Triangle(matrix.data, i, size=size) )
            downs = as.vector( Get.Downstream.Triangle(matrix.data, i, size=size) )
            
            wil.test =  stats::wilcox.test(x=dia*scale, y=c(ups, downs), alternative="less", exact=FALSE)
            pvalue[i] = wil.test$p.value  
            
        }
        
        pvalue[ is.na(pvalue) ] = 1
        return(pvalue)
    }
    
    
    Get.Upstream.Triangle <- function(mat.data, i, size)
    {
        n_bins = nrow(mat.data)
        
        lower = max(1, i-size)
        tmp.mat = mat.data[lower:i, lower:i]
        return( tmp.mat[ upper.tri( tmp.mat, diag=FALSE ) ] )
    }
    
    
    Get.Downstream.Triangle <- function(mat.data, i, size)
    {
        n_bins = nrow(mat.data)
        if(i==n_bins) return(NA)
        
        upperbound = min(i+size, n_bins)
        tmp.mat = mat.data[(i+1):upperbound, (i+1):upperbound]
        return( tmp.mat[ upper.tri( tmp.mat, diag=FALSE ) ] )
    }
    
    Get.Diamond.Matrix2 <- function(mat.data, i, size)
    {
        n_bins = nrow(mat.data)
        new.mat = matrix(rep(NA, size*size), nrow=size, ncol=size)
        
        for(k in seq(1,size,1))
        {
            if(i-(k-1) >= 1 && i < n_bins)
            {
                lower = min(i+1, n_bins)
                upper = min(i+size, n_bins)
                
                new.mat[size-(k-1), seq(1,(upper-lower+1),1)] = mat.data[i-(k-1), lower:upper]
            }
        }
        
        return(new.mat)
    }
    
    
    Convert.Bin.To.Domain <- function(bins, signal.idx, gap.idx, pvalues=NULL, pvalue.cut=NULL)
    {
        n_bins = nrow(bins)
        ret = data.frame(chr=character(0), from.id=numeric(0), from.coord=numeric(0), to.id=numeric(0), to.coord=numeric(0), tag=character(0), size=numeric(0))
        levels( x=ret[, "tag"] ) = c("domain", "gap", "boundary")
        
        rmv.idx = setdiff(seq(1,n_bins,1), gap.idx)
        proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
        from.coord = bins[proc.region[, "start"], "from.coord"]
        n_procs = nrow(proc.region)
        if(n_procs>0){
            gap = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("gap", n_procs), size=rep(0, n_procs), stringsAsFactors=FALSE)
        }else{
            gap=data.frame(stringsAsFactors = FALSE)
        }
        rmv.idx = union(signal.idx, gap.idx)
        proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
        n_procs = nrow(proc.region)
        from.coord = bins[proc.region[, "start"], "from.coord"]
        domain = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("domain", n_procs), size=rep(0, n_procs), stringsAsFactors=FALSE)
        
        rmv.idx = setdiff(seq(1,n_bins,1), signal.idx)
        proc.region = as.data.frame( Which.process.region(rmv.idx, n_bins, min.size=1) )
        n_procs = nrow(proc.region)
        if(n_procs>0)
        {
            from.coord = bins[proc.region[, "start"]+1, "from.coord"]  
            boundary = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("boundary", n_procs), size=rep(0, n_procs), stringsAsFactors=FALSE)
            ret = rbind(ret, boundary)
        }
        
        ret = rbind(gap, domain)
        ret = ret[order(ret[,3]), ]
        
        ret[, "to.coord"] = c(ret[2:nrow(ret), "from.coord"], bins[n_bins, "to.coord"])
        ret[, "from.id"] = match( ret[, "from.coord"], bins[, "from.coord"] )
        ret[, "to.id"] = match(ret[, "to.coord"], bins[, "to.coord"])
        ret[, "size"] = ret[,"to.coord"]-ret[,"from.coord"]
        
        if(!is.null(pvalues) && !is.null(pvalue.cut))
        {
            for(i in seq(1,nrow(ret),1))
            {
                if(ret[i, "tag"]=="domain")
                {
                    if (is.na(ret[i,"to.id"])){
                        ret[i,"to.id"]=ret[i+1,"from.id"]-1 
                    }
                    if (is.na(ret[i,"from.id"])){
                        if (i==1){ ret[i,"from.id"]=1
                        }else{
                            ret[i,"from.id"]=ret[i-1,"to.id"]+1 
                        }
                    }
                    domain.bins.idx = ret[i, "from.id"]:ret[i, "to.id"]
                    p.value.constr = which( pvalues[ domain.bins.idx ] < pvalue.cut )
                    
                    if( length(domain.bins.idx) == length(p.value.constr)) ret[i, "tag"] = "boundary"
                }
            }
        }
        
        return(ret)
    }
    #main body
        if (inherits(matrix.file, "TopDomData")) {
            bins <- matrix.file$bins
            matrix.data <- matrix.file$counts
            n_bins <- nrow(bins)
            mean.cf <- rep(0, times = n_bins)
            pvalue <- rep(1.0, times = n_bins)
            local.ext <- rep(-0.5, times = n_bins)
        } else {
            if(verbose) message("Step 0 : File Read ")
            window.size = as.numeric(window.size)
            matdf <- as.data.frame(data.table::fread(matrix.file, header=FALSE))
            
            if( ncol(matdf) - nrow(matdf) == 3) {
                colnames(matdf) <- c("chr", "from.coord", "to.coord")
            } else if( ncol(matdf) - nrow(matdf) ==4 ) {
                colnames(matdf) <- c("id", "chr", "from.coord", "to.coord")
            } else {
                stop("Unknown Type of matrix file")
            }
            n_bins = nrow(matdf)
            mean.cf <- rep(0, n_bins)
            pvalue <- rep(1, n_bins)
            
            local.ext = rep(-0.5, n_bins)
            
            bins <- data.frame(id=seq(1,n_bins,1), 
                               chr=matdf[, "chr"], 
                               from.coord=matdf[, "from.coord"], 
                               to.coord=matdf[, "to.coord"] )
            matrix.data <- as.matrix( matdf[, (ncol(matdf) - nrow(matdf)+1 ):ncol(matdf)] )
            
            if(verbose) message("Step 0 : Done !!")
        }
        
        
        if(verbose){
        message("Step 1 : Generating binSignals by computing bin-level contact frequencies")
        ptm <- proc.time()
        }
        for(i in seq(1,n_bins,1))
        {
            diamond = Get.Diamond.Matrix(mat.data=matrix.data, i=i, size=window.size)
            mean.cf[i] = mean(diamond)
        }
        
    if(verbose){
        eltm = proc.time() - ptm
        msg<-paste("Step 1 Running Time : ", eltm[3])
        message(msg)
        message("Step 1 : Done !!")
        message("Step 2 : Detect TD boundaries based on binSignals")
        ptm = proc.time()
    }

        gap.idx = Which.Gap.Region(matrix.data=matrix.data, window.size)
        
        proc.regions = Which.process.region(rmv.idx=gap.idx, n_bins=n_bins, min.size=3)
        
        
        for( i in seq(1,nrow(proc.regions),1))
        {
            start = proc.regions[i, "start"]
            end = proc.regions[i, "end"]
            
            if (verbose){
                msg<-paste("Process Regions from ", start, "to", end)
                message(msg)
            }
            
            local.ext[start:end] = Detect.Local.Extreme(x=mean.cf[start:end])
        }
        
        if(verbose){
        eltm = proc.time() - ptm
        msg<-paste("Step 2 Running Time : ", eltm[3])
        message(msg)
        message("Step 2 : Done !!")
        }
        
        if(statFilter)
        {
            if(verbose){
            message("Step 3 : Statistical Filtering of false positive TD boundaries")
            ptm = proc.time()  
            message("-- Matrix Scaling....")
            }
            scale.matrix.data = matrix.data
            for( i in seq(1,(2*window.size),1) )
            {
                scale.matrix.data[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ] = scale( matrix.data[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ] )
            }
            
            if(verbose) message("-- Compute p-values by Wilcoxon Rank Sum Test")
            for( i in seq(1,nrow(proc.regions),1))
            {
                start = proc.regions[i, "start"]
                end = proc.regions[i, "end"]
                
                if(verbose) {
                    msg<-paste("Process Regions from ", start, "to", end)
                    message(msg)
                }
                
                pvalue[start:end] <- Get.Pvalue(matrix.data=scale.matrix.data[start:end, start:end], size=window.size, scale=1)
            }
            if(verbose){
            message("-- Done!")
            message("-- Filtering False Positives")
            }
            local.ext[intersect( union(which( local.ext==-1), which(local.ext==-1)), which(pvalue<0.05))] = -2
            local.ext[which(local.ext==-1)] = 0
            local.ext[which(local.ext==-2)] = -1
            if(verbose){
            message("-- Done!")
            eltm = proc.time() - ptm
            msg<-paste("Step 3 Running Time : ", eltm[3])
            message(msg)
            message("Step 3 : Done!")
            }
        } else pvalue = 0
        
        domains = Convert.Bin.To.Domain(bins=bins, 
                                            signal.idx=which(local.ext==-1), 
                                            gap.idx=which(local.ext==-0.5), 
                                            pvalues=pvalue, 
                                            pvalue.cut=0.05)
        
        bins = cbind(bins, 
                     local.ext = local.ext,
                     mean.cf = mean.cf, 
                     pvalue = pvalue)
        
        bedform = domains[, c("chr", "from.coord", "to.coord", "tag")]
        colnames(bedform) = c("chrom", "chromStart", "chromEnd", "name")
        
        if( !is.null(outFile) ) {
            outBinSignal =  path.expand(paste(outFile, ".binSignal", sep=""))
            if(verbose){
            message("Writing Files")
            msg<-paste("binSignal File :", outBinSignal)
            message(msg)
            }
            outBinSignaldir<-gsub("/[^/]+$", "",outBinSignal)
            if (outBinSignaldir==outBinSignal){
                outBinSignaldir<-gsub("\\[^\\]+$", "",outBinSignal)
            }
            if (outBinSignaldir==outBinSignal){
                outBinSignaldir<-gsub("\\\\[^\\\\]+$", "",outBinSignal)
            }
            if (!outBinSignaldir==outBinSignal&!dir.exists(outBinSignaldir)){
                dir.create(outBinSignaldir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            }
            data.table::fwrite(bins, file=outBinSignal, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")  
            
            
            outDomain = path.expand(paste(outFile, ".domain", sep=""))
            if(verbose) {
                msg<-paste("Domain File :", outDomain)
                message(msg)
            }
            outDomaindir<-gsub("/[^/]+$", "",outDomain)
            if (outDomaindir==outDomain){
                outDomaindir<-gsub("\\[^\\]+$", "",outDomain)
            }
            if (outDomaindir==outDomain){
                outDomaindir<-gsub("\\\\[^\\\\]+$", "",outDomain)
            }
            if (!outDomaindir==outDomain&!dir.exists(outDomaindir)){
                dir.create(outDomaindir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            }
            data.table::fwrite( domains, file=outDomain, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
            
            outBed = path.expand(paste(outFile, ".bed", sep=""))
            if(verbose) {
                msg<-paste("Bed File : ", outBed)
                message(msg)
            }
            outBeddir<-gsub("/[^/]+$", "",outBed)
            if (outBeddir==outBed){
                outBeddir<-gsub("\\[^\\]+$", "",outBed)
            }
            if (outBeddir==outBed){
                outBeddir<-gsub("\\\\[^\\\\]+$", "",outBed)
            }
            if (!outBeddir==outBed&!dir.exists(outBeddir)){
                dir.create(outBeddir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            }
            data.table::fwrite( bedform, file=outBed, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
        }
        if(verbose){
        message("Job Complete !")
        }
        return(list(binSignal=bins, domain=domains, bed=bedform))
    }


#'.readDiffInputFiles
#'
#'This function prepares DESeq2 normalization factors and counts across conditions and replicates
#'@param conds conditions set, (names of input paths)
#'@param input_paths degrees of freedom for the genomic distance spline function if \code{distance_type='spline'}. Defaults to 6, which corresponds to a cubic spline as explained in Carty et al. (2017)
#'@param sigs filter_file data on hicdcdiff in data.table format
#'@param chrom name of the chromosome
#'@param Dmin minimum distance (included) to check for significant interactions,
#'defaults to 0. Put Dmin=1 to ignore D=0 bins in calculating
#'normalization factors.
#'@param Dmax maximum distance (included) to check for significant interactions,
#'defaults to 2e6 or maximum in the data; whichever is minimum.
#'@param dband set of distance bins in which normalization factors are calculated
#'@param bin_type 'Bins-uniform' if uniformly binned by binsize in
#'bp, or 'Bins-RE-sites' if binned by number of
#'restriction enzyme fragment cutsites!
#'@param binsize binsize in bp if bin_type='Bins-uniform' (or number of
#'RE fragments if bin_type='Bins-RE-sites'), e.g., default 5000
#'@noRd
    .readDiffInputFiles <- function(conds, input_paths, sigs, chrom, Dmin, Dmax, dband, bin_type, binsize) {
        retlist <- list()
        normfac <- NULL
        for (cond in conds) {
            numrep <- length(input_paths[[cond]])
            for (i in seq(1, numrep)) {
                prefix <- input_paths[[cond]][i]
                if (!grepl("\\.hic$", prefix, ignore.case = TRUE)){
                if (is.list(prefix)&methods::is(prefix[[1]], "GInteractions")){
                    gi_list_validate(prefix)
                    normfac_add <- prefix[[chrom]]
                    rm(prefix)
                } else {
                gi_list <- gi_list_read(path.expand(prefix))
                gi_list_validate(gi_list)
                normfac_add <- gi_list[[chrom]]
                rm(gi_list)
                }
                normfac_add <- data.frame(chr = chrom, startI = BiocGenerics::start(InteractionSet::anchors(normfac_add)$first), 
                    startJ = BiocGenerics::start(InteractionSet::anchors(normfac_add)$second), counts = mcols(normfac_add)$counts, 
                    D = mcols(normfac_add)$D, stringsAsFactors = FALSE)
                normfac_add <- normfac_add %>% dplyr::filter(.data$counts > 0)
                } else if (grepl("\\.hic$", prefix, ignore.case = TRUE)) {
                    if (bin_type == "Bins-uniform") {
                    if (!.Platform$OS.type=="windows"){
                    normfac_add <- tryCatch(
                            straw(norm = "NONE", fn = path.expand(prefix), bs = binsize, ch1 = gsub("chr", "", chrom), ch2 = gsub("chr", "", chrom), 
                                  u = "BP"),
                            error=function(e){
                                tryCatch(straw(norm = "NONE", fn = path.expand(prefix), bs = binsize, ch1 = chrom, ch2 = chrom, 
                                               u = "BP"),
                                         error=function(e){
                                             straw_dump(norm = "NONE",fn=path.expand(prefix),bs=binsize,ch1=gsub("chr", "", chrom),ch2=gsub("chr", "", chrom),u="BP")   
                                         })
                            })
                    }else{
                        normfac_add<-straw_dump(norm = "NONE",fn=path.expand(prefix),bs=binsize,ch1=gsub("chr", "", chrom),ch2=gsub("chr", "", chrom),u="BP")   
                        gc(reset=TRUE,full=TRUE)
                    }
                    } else {
                    if (!.Platform$OS.type=="windows"){
                    normfac_add <- tryCatch(
                        straw(norm = "NONE", fn = path.expand(prefix), bs = binsize, ch1 = gsub("chr", "", chrom), ch2 = gsub("chr", "", chrom), 
                             u = "FRAG"),
                        error=function(e){
                            tryCatch(straw(norm = "NONE", fn = path.expand(prefix), bs = binsize, ch1 = chrom, ch2 = chrom, 
                                    u = "FRAG"),
                                    error=function(e){
                                         straw_dump(norm = "NONE",fn=path.expand(prefix),bs=binsize,ch1=gsub("chr", "", chrom),ch2=gsub("chr", "", chrom),u="FRAG")   
                                     })
                        })
                    }else{
                        normfac_add<-straw_dump(norm = "NONE",fn=path.expand(prefix),bs=binsize,ch1=gsub("chr", "", chrom),ch2=gsub("chr", "", chrom),u="FRAG")   
                        gc(reset=TRUE,full=TRUE)
                    }
                    }
                    colnames(normfac_add) <- c("startI", "startJ", "counts")
                    normfac_add <- normfac_add %>% dplyr::mutate(chr = chrom, D = abs(.data$startI - .data$startJ))
                } else {
                    stop(paste0("File not found relating to ", prefix))
                }
                normfac_add <- normfac_add %>% dplyr::filter(.data$D >= Dmin & .data$D <= Dmax) %>% dplyr::mutate(Dband = as.numeric(cut(.data$D, 
                breaks = dband, labels = seq(1, (length(dband) - 1), 1), include.lowest = TRUE))) %>% dplyr::select(.data$Dband, 
                .data$chr, .data$startI, .data$startJ, .data$counts)
                # rename counts column
                colnames(normfac_add)[colnames(normfac_add) %in% "counts"] <- paste0(cond, ".", i)
                if (is.null(normfac)) {
                normfac <- normfac_add
                } else {
                # join counts across conditions and replicates
                normfac <- dplyr::inner_join(normfac, normfac_add)
                }
            }
        }
        # calculate geometric mean of all count columns
        countcols <- colnames(normfac)[!colnames(normfac) %in% c("chr", "startI", "startJ", "Dband")]
        normfac$geomean <- exp(rowMeans(log(normfac[countcols] + 0.5)))
        # calculate factors: count/geometric mean
        for (countcol in countcols) {
            normfac[, paste0(countcol, ".fac")] <- (normfac[, countcol] + 0.5)/normfac$geomean
        }
        # calculate median factor by distance band and store it as normfac.final
        countcols <- paste0(countcols, ".fac")
        normfac.final <- normfac %>% dplyr::group_by(.data$Dband) %>% dplyr::summarise_at(.vars = countcols, .funs = stats::median)
        normfac.final$geomean <- exp(rowMeans(log(normfac.final[countcols])))
        # join back with normfac with normfac.final columns .fac replaced with .norm and geomean removed
        colnames(normfac.final) <- gsub(".fac", ".norm", colnames(normfac.final))
        normfac <- dplyr::left_join(normfac, normfac.final %>% dplyr::select(-.data$geomean))
        # significant bins to be filtered
        sigbins <- sigs %>% dplyr::filter(.data$chr == chrom)
        sigbins <- paste0(sigbins$startI, ":", sigbins$startJ)
        retlist[["normfac"]] <- normfac
        retlist[["normfac.final"]] <- normfac.final
        retlist[["sigbins"]] <- sigbins
        return(retlist)
    }

#'.plotNormalizationFactors
#'
#'This function plots normalization factors
#'@param normfac.final normalization factors from retlist returned from .readDiffInputFiles
#'@param input_paths degrees of freedom for the genomic distance spline function if \code{distance_type='spline'}. Defaults to 6, which corresponds to a cubic spline as explained in Carty et al. (2017)
#'@param binsize binsize in basepair
#'@param chr name of the chromosome
#'@param dband set of distance bins in which normalization factors are calculated
#'@param conds conditions set, (names of input paths)
#'@param output_path the path to the folder and name prefix you want to
#'place DESeq-processed matrices (in a .txt file), plots
#'(if \code{diagnostics=TRUE}) and DESeq2 objects (if \code{DESeq.save=TRUE}).
#'Files will be generated for each chromosome.
#'@noRd
    .plotNormalizationFactors <- function(normfac.final, binsize, chr, dband, conds, output_path) {
        plotpath <- path.expand(paste0(output_path, "sizefactors_", chr, ".pdf"))
        plotpaths <- plotpath
        plotpathdir<-gsub("/[^/]+$", "",plotpath)
        if (plotpathdir==plotpath){
            plotpathdir<-gsub("\\[^\\]+$", "",plotpath)
        }
        if (plotpathdir==plotpath){
            plotpathdir<-gsub("\\\\[^\\\\]+$", "",plotpath)
        }
        if (!plotpathdir==plotpath&!dir.exists(plotpathdir)){
            dir.create(plotpathdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
        }
        grDevices::pdf(plotpath)
        graphics::par(mfrow = c(1, length(conds)))
        for (cond in conds) {
            numrep <- length(grep(paste0("^", cond, ".*norm$"), colnames(normfac.final)))
            condn <- paste0(cond, ".1.norm")
            num <- grep(condn, colnames(normfac.final))
            graphics::plot(dband[normfac.final$Dband], normfac.final[[condn]], xlab = "distance", ylab = "norm factors", ylim = range(normfac.final[, 
                num:(num + (numrep - 1))]), main = paste0(cond))
            if (numrep > 1) {
                for (i in seq(2, numrep, 1)) {
                condn <- paste0(cond, ".", i, ".norm")
                graphics::points(dband[normfac.final$Dband], normfac.final[[condn]], col = "red")
                }
            }
        }
        grDevices::dev.off()
        plotpath <- path.expand(paste0(output_path, "geomean_sizefactors_", chr, ".pdf"))
        plotpaths <- c(plotpaths, plotpath)
        plotpathdir<-gsub("/[^/]+$", "",plotpath)
        if (plotpathdir==plotpath){
            plotpathdir<-gsub("\\[^\\]+$", "",plotpath)
        }
        if (plotpathdir==plotpath){
            plotpathdir<-gsub("\\\\[^\\\\]+$", "",plotpath)
        }
        if (!plotpathdir==plotpath&!dir.exists(plotpathdir)){
            dir.create(plotpathdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
        }
        grDevices::pdf(plotpath)
        graphics::par(mfrow = c(1, 1))
        graphics::plot(dband[normfac.final$Dband], normfac.final$geomean, xlab = "distance", ylab = "geometric mean")
        grDevices::dev.off()
        return(plotpaths)
    }

#'.get_cache
#'
#'This function retrieves the BioCFileCache for this package 
#'@noRd  
    .get_cache <-
        function()
        {
            cache <- rappdirs::user_cache_dir(appname="HiCDCPlus")
            BiocFileCache::BiocFileCache(cache, ask = FALSE)
        }

#'.get_cache
#'
#'This function downloads juicer_tools as necessary
#'@noRd
.download_juicer <-
    function()
{
    fileURL <- "https://github.com/aidenlab/Juicebox/releases/download/v.2.13.07/juicer_tools.jar"

    bfc <- .get_cache()
    rid <- BiocFileCache::bfcquery(bfc, "juicer_tools", "rname")$rid
    if (!length(rid)) {
     message( "Downloading Juicer Tools" )
     rid <- names(BiocFileCache::bfcadd(bfc, "juicer_tools", fileURL ))
    }
    if (!isFALSE(BiocFileCache::bfcneedsupdate(bfc, rid)))
    BiocFileCache::bfcdownload(bfc, rid)

    BiocFileCache::bfcrpath(bfc, rids = rid)
    }

#'.hic2htcexp
#'
#'This function converts a .hic file into an HTCexp object
#'@param chrom1 first chromosome name
#'@param chrom2 second chromosome name
#'@param binsize binsize in basepair
#'@param hic_path path to the hic file
#'@param gen name of the species: e.g., default \code{'Hsapiens'}
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}
#'@noRd
.hic2htcexp<-function(chrom1,chrom2,binsize,hic_path,gen,gen_ver){
  chr_select1 <- gsub("chr", "", chrom1)
  chr_select2 <- gsub("chr", "", chrom2)
  if (!.Platform$OS.type=="windows"){
    count_matrix <- tryCatch(
      straw(norm = "NONE", fn = path.expand(hic_path), bs = binsize, ch1 = chr_select1, ch2 = chr_select2, 
            u = "BP"),
      error=function(e){
        tryCatch(straw(norm = "NONE", fn = path.expand(hic_path), bs = binsize, ch1 = chrom1, ch2 = chrom2, 
                       u = "BP"),
                 error=function(e){
                   straw_dump(norm = "NONE",fn=path.expand(hic_path),bs=binsize,ch1=chr_select1,ch2=chr_select2,u="BP")   
                 })
      })
  }else{
    count_matrix<-straw_dump(norm = "NONE",fn=path.expand(hic_path),bs=binsize,ch1=chr_select1,ch2=chr_select2,u="BP")   
  }
  xgi<-generate_binned_gi_list(binsize,chrs=chrom1,gen=gen,gen_ver=gen_ver)[[chrom1]]@regions
  ygi<-generate_binned_gi_list(binsize,chrs=chrom2,gen=gen,gen_ver=gen_ver)[[chrom2]]@regions
  names(xgi)<-paste0(GenomicRanges::seqnames(xgi),':',
                     GenomicRanges::start(xgi),'-',
                     GenomicRanges::end(xgi))
  names(ygi)<-paste0(GenomicRanges::seqnames(ygi),':',
                     GenomicRanges::start(ygi),'-',
                     GenomicRanges::end(ygi))
  minrow1<-min(GenomicRanges::start(xgi))
  minrow2<-min(GenomicRanges::start(ygi))
  maxrow1<-max(GenomicRanges::start(xgi))
  maxrow2<-max(GenomicRanges::start(ygi))
  count_matrix<-count_matrix%>%dplyr::filter(.data$x<=maxrow1&.data$y<=maxrow2)
  if (sum(count_matrix$x==maxrow1&count_matrix$y==maxrow2)<=0){
    count_matrix<-dplyr::bind_rows(count_matrix,data.frame(x=maxrow1,y=maxrow2,
                                                           counts=0))
  }
  if (sum(count_matrix$x==minrow1&count_matrix$y==minrow2)<=0){
    count_matrix<-dplyr::bind_rows(count_matrix,data.frame(x=minrow1,y=minrow2,
                                                           counts=0))
  }
  count_matrix$x<-subjectHits(GenomicRanges::findOverlaps(
    GenomicRanges::GRanges(chrom1,
                           IRanges::IRanges(count_matrix$x,count_matrix$x+1)),xgi,minoverlap=2))
  count_matrix$y<-subjectHits(GenomicRanges::findOverlaps(
    GenomicRanges::GRanges(chrom2,
                           IRanges::IRanges(count_matrix$y,count_matrix$y+1)),ygi,minoverlap=2))
  count_matrix<-count_matrix%>%dplyr::arrange(.data$x,.data$y)
  htc<-HiTC::normICE(HiTC::forceSymmetric(HiTC::HTCexp(
    intdata=Matrix::sparseMatrix(
      j=count_matrix$x,
      i=count_matrix$y,
      x=count_matrix$counts,symmetric = TRUE),
    ygi=ygi,
    xgi=xgi)),eps=1e-2, sparse.filter=0.025, max_iter=10)
  return(htc)
}
