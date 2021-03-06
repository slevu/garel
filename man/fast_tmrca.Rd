% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnabin_functions.R
\name{fast_tmrca}
\alias{fast_tmrca}
\title{Calculate TMRCA for one sequence with few same-subtype external references and other-subtype outgroup}
\usage{
fast_tmrca(pathlst = "lst.rds", outdir = ".", plots = FALSE)
}
\arguments{
\item{pathlst}{Path to list as .rds containing [seq = full MSA, suk = label of sequence to test, sla = labels of reference sequences, og = labels of outgroup sequences)}

\item{outdir}{Directory for output}

\item{plots}{Logical to draw plots}
}
\description{
Used in \code{\link{do_tmrca_outliers_splits}}
}
