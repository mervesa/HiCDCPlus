% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_hic_counts.R
\name{add_hic_counts}
\alias{add_hic_counts}
\title{add_hic_counts}
\usage{
add_hic_counts(gi_list, hic_path, chrs = NULL, add_inter = FALSE)
}
\arguments{
\item{gi_list}{valid, uniformly binned gi_list instance. 
See \code{?gi_list_validate} and \code{gi_list_binsize_detect} for details.}

\item{hic_path}{path to the .hic file}

\item{chrs}{a subset of chromosomes' e.g., c('chr21','chr22'). Defaults
to all chromosomes in the \code{gi_list} instance.}

\item{add_inter}{Interchromosomal interaction counts added as a 1D feature
named 'inter' on regions metadata handle of each gi_list element (e.g., 
\code{gi_list[[1]]@regions@elementMetadata} or not;
default FALSE}
}
\value{
\code{gi_list} instance with counts on the metadata (e.g., 
\code{mcols(gi_list[[1]])} handle on each list element, and 'inter' on 
regions metadata handle of each element if \code{add_inter=TRUE}.
}
\description{
This function adds counts from a .hic file into a valid, binned, gi_list
instance.
}
\examples{
gi_list<-generate_binned_gi_list(50e3,chrs='chr22')
gi_list<-add_hic_counts(gi_list,
hic_path=system.file("extdata", "GSE63525_HMEC_combined_example.hic",
package = "HiCDCPlus"))
}
