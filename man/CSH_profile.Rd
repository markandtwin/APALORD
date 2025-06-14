% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/APA_within_PAS.R
\name{CSH_profile}
\alias{CSH_profile}
\title{Check the 3'end distribution within the same PolyA site.}
\usage{
CSH_profile(
  reads,
  PAU_data,
  control,
  experimental,
  cores = 1,
  min_reads = 20,
  direct_RNA = F
)
}
\arguments{
\item{reads}{information from RNAseq samples}

\item{PAU_data}{PAU information from previous analysis}

\item{control}{which group in the data is used as the control group}

\item{experimental}{which group in the data is used as the experimental group}

\item{cores}{number of threads used for the computation}

\item{min_reads}{minium reads count required at PAS for heterogenity analysis}

\item{direct_RNA}{whether or not the data is direct RNAseq}
}
\value{
a table showing the uasage of top PASs for each single gene with enough depth in the dataset of each sample
}
\description{
This function checks the 3'end distribution within individual PAS to detect cleavage site shift.
}
