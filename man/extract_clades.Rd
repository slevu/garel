% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_trees.R
\name{extract_clades}
\alias{extract_clades}
\title{Extract clades from list of LSD bootstrap trees}
\usage{
extract_clades(lsd_results, dirclades = "data/CLADES",
  max.clade.size = 1000, verbose = TRUE)
}
\arguments{
\item{lsd_results}{List of LSD results pathname (\code{*.rds})}

\item{dirclades}{Where to save clades}

\item{max.clade.size}{For \code{\link{find_clades_slice}}, maximum clade size}

\item{verbose}{Details}
}
\value{
List of clades filenames
}
\description{
Extract clades from list of LSD bootstrap trees
}
