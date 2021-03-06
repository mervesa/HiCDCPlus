% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_bintolen_gi_list.R
\name{generate_bintolen_gi_list}
\alias{generate_bintolen_gi_list}
\title{generate_bintolen_gi_list}
\usage{
generate_bintolen_gi_list(
  bintolen_path,
  chrs = NULL,
  Dthreshold = 2e+06,
  binned = TRUE,
  binsize = NULL,
  gen = "Hsapiens",
  gen_ver = "hg19"
)
}
\arguments{
\item{bintolen_path}{path to the flat file containing columns named bins and
features}

\item{chrs}{select a subset of chromosomes' e.g.,
c('chr21','chr22'). Defaults to all chromosomes
specified in the bintolen file.}

\item{Dthreshold}{maximum distance (included) to check for significant
interactions, defaults to 2e6 or maximum in the data; whichever is smaller.}

\item{binned}{TRUE if the bintolen file is uniformly binned. Defaults to TRUE.}

\item{binsize}{bin size in bp to be generated for the object. Defaults to the
binsize in the bintolen file, if exists.}

\item{gen}{name of the species: e.g., default \code{'Hsapiens'}}

\item{gen_ver}{genomic assembly version: e.g., default \code{'hg19'}}
}
\value{
a valid gi_list instance with genomic features derived from specified
restriction enzyme cut patterns when generating the bintolen file using
\code{construct_features} (see \code{?construct_features} for help).  
Genomic 1D features are stored in the regions metadata handle
of each list element (e.g., \code{gi_list[[1]]@regions@elementMetadata}).
}
\description{
Generates a gi_list instance from a bintolen file generated by
\code{generate.features} (see \code{?generate.features}) for details).
}
\examples{
chrs<-'chr22'
bintolen_path<-system.file("extdata", "test_bintolen.txt.gz",
package = "HiCDCPlus")
gi_list<-generate_bintolen_gi_list(bintolen_path,chrs)
}
