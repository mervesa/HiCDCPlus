#' gi.list.read.R
#'
#'Reads a written gi_list instance using \code{gi.list.write} into a valid
#'gi_list instance.
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param fname path to the file to read from (can end with .txt,
#'.rds, or .txt.gz).
#'@param chrs select a subset of chromosomes' e.g.,
#'c('chr21','chr22'). Defaults to all chromosomes contained in the \code{fname}.
#'@param Dthreshold maximum distance (included) to check for significant
#'interactions, defaults to the maximum in the data.
#'@param features Select the subset of features (1-D or 2-D) to be added to the
#'gi_list instance (without the trailing I or J),
#'defaults to all features (score column gets ingested as 'score'). 
#'@param gen name of the species: e.g., default \code{'Hsapiens'}
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}
#'@return A valid gi_list instance with 1D features stored in regions metadata
#'handle of each list element (e.g., 
#'\code{gi_list[[1]]@regions@elementMetadata}) in the instance and with 2D 
#'features stored in metadata handle
#'(i.e., \code{mcols(gi)}).
#'@examples 
#'outputdir<-paste0(tempdir(check=TRUE),'/')
#'gi_list<-generate.binned.gi.list(1e6,chrs='chr22')
#'gi.list.write(gi_list,paste0(outputdir,'testgiread.txt'))
#'gi_list2<-gi.list.read(paste0(outputdir,'testgiread.txt'))
#'@export

gi.list.read <- function(fname, chrs = NULL, Dthreshold=NULL, features=NULL,
    gen="Hsapiens",gen_ver="hg19") {
    input.file.read <- function(filepath) {
        # reads files of RDS, .txt, or .txt.gz filepath: valid path ending in
        #.txt, .txt.gz, or .rds
        if (grepl("\\.txt.gz$", filepath) | grepl("\\.txt$", filepath)) {
            return(data.table::fread(filepath))
        } else if (grepl("\\.rds$", filepath)) {
            return(readRDS(filepath))
        } else {
            stop("Can only read in paths ending with .txt,.txt.gz, or .rds")
        }
    }
    data<-input.file.read(fname)
    if(is.list(data)&methods::is(data[[1]], "GInteractions")){
    gi.list.validate(data)
    return(data)
    }
    if(!all(c('chrI','startI','chrJ','startJ','D')%in%colnames(data))){
        stop("Some of columns 'chrI','startI','chrJ','startJ','D' do not
            exist in file")
    }
    if(is.null(chrs)) chrs<-unique(c(data$chrI,data$chrJ))
    if(is.null(Dthreshold)) Dthreshold<-max(data$D)
    if(is.null(features)){
        features<-colnames(data)[!colnames(data)%in%c('chrI','startI',
            'chrJ','startJ','D','strI','strJ','fragI','fragJ','endI','endJ')]
        features<-unique(gsub("J$","",gsub("I$","",features)))
    }
    data<-data%>%dplyr::filter(.data$chrI%in%chrs&.data$chrJ%in%chrs&
            .data$D<=Dthreshold)
    if(!"endI"%in%colnames(data)){
    df<-dplyr::bind_rows(data%>%dplyr::select(.data$chrI,.data$startI)%>%
            dplyr::rename('chr'='chrI','start'='startI'),
        data%>%dplyr::select(.data$chrJ,.data$startJ)%>%
            dplyr::rename('chr'='chrJ','start'='startJ'))%>%dplyr::distinct()
    }else{
    df<-dplyr::bind_rows(data%>%dplyr::select(.data$chrI,.data$startI,.data$endI)%>%
            dplyr::rename('chr'='chrI','start'='startI','end'='endI'),
        data%>%dplyr::select(.data$chrJ,.data$startJ,.data$endJ)%>%
            dplyr::rename('chr'='chrJ','start'='startJ','end'='endJ'))%>%dplyr::distinct()        
    }
    gi_list<-generate.df.gi.list(df,chrs=chrs,Dthreshold=Dthreshold,
        gen=gen,gen_ver=gen_ver)
    #detect and add 1D features
    features.1D<-colnames(data)[grepl("I$",colnames(data))|grepl("J$",colnames(data))]
    features.1D<-unique(gsub("J$","",gsub("I$","",features.1D)))
    features.1D<-features.1D[features.1D%in%features]
    if(length(features.1D)>0){
    df<-as.data.frame(data)[c("chrI","startI",paste0(features.1D,"I"))]%>%
            dplyr::rename_all(function(x) gsub("I", "", x))%>%dplyr::distinct()
    gi_list<-add.1D.features(gi_list,df,chrs=chrs,features=features.1D)
    }
    #detect and add 2D features
    features.2D<-colnames(data)[!(grepl("I$",colnames(data))|grepl("J$",colnames(data)))]
    features.2D<-features.2D[features.2D%in%features]
    features.2D<-features.2D[!features.2D%in%features.1D]
    if(length(features.2D)>0){
    data<-as.data.frame(data)[c("chrI","startI","startJ",features.2D)]%>%
        dplyr::rename("chr"="chrI")
    for (chrom in chrs){
        data_chr<-data%>%dplyr::filter(.data$chr==chrom)
        gi_list[[chrom]]<-add.2D.features(gi_list[[chrom]],data_chr,features=features.2D)
    }
    }
    return(gi_list)
}
