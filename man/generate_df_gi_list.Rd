% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_df_gi_list.R
\name{generate_df_gi_list}
\alias{generate_df_gi_list}
\title{generate_df_gi_list}
\usage{
generate_df_gi_list(
  df,
  chrs = NULL,
  Dthreshold = 2e+06,
  gen = "Hsapiens",
  gen_ver = "hg19"
)
}
\arguments{
\item{df}{DataFrame with columns named 'chr', 'start', (and optionally
'end', if the regions have gaps) and 1D features with their respective
column names.}

\item{chrs}{select a subset of chromosomes' e.g.,
c('chr21','chr22'). Defaults to all chromosomes
specified in \code{df}.}

\item{Dthreshold}{maximum distance (included) to check for significant
interactions, defaults to 2e6 or maximum in the data, whichever is smaller.}

\item{gen}{name of the species: e.g., default \code{'Hsapiens'}}

\item{gen_ver}{genomic assembly version: e.g., default \code{'hg19'}}
}
\value{
a valid gi_list instance with genomic features supplied from 
\code{df}. Genomic 1D features are stored in the regions metadata handle
of each list element (e.g., \code{gi_list[[1]]@regions@elementMetadata}).
}
\description{
Generates a gi_list instance from a data frame object describing the regions.
}
\examples{
df<-data.frame(chr='chr9',start=seq(1e6,10e6,1e6))
gi_list<-generate_df_gi_list(df)
}
