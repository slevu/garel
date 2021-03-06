% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnabin_functions.R
\name{prep_blast_align}
\alias{prep_blast_align}
\title{Process MSA and its corresponding BLAST match file for alignment}
\usage{
prep_blast_align(infile = "data/MSM_CRF02AG_splits/MSM_CRF02AG_nodrm.rds",
  subtype = FALSE, verbose = TRUE)
}
\arguments{
\item{infile}{Path to MSA}

\item{subtype}{Describe subtype distribution of hits}

\item{verbose}{explain stuff}
}
\value{
List of path to fasta files
}
\description{
create two fasta files
}
