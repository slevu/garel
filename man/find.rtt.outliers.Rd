% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnabin_functions.R
\name{find.rtt.outliers}
\alias{find.rtt.outliers}
\title{Root-to-tip regression outliers}
\usage{
find.rtt.outliers(tree, dates, mult = 2, plots = TRUE, title = NULL)
}
\arguments{
\item{tree}{phylo tree}

\item{dates}{with tip names}

\item{mult}{multiplier of sd to define outliers}

\item{plots}{plot}

\item{title}{title of plot}
}
\value{
output list of outliers from mult*SD and from grubbs
}
\description{
Superseded by \code{find.rtt.outliers2}, with cpp distance calculation
}
