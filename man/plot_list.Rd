% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plot_list}
\alias{plot_list}
\title{Make a list of plots and arrange in a grid}
\usage{
plot_list(lst, funplot, ncols = 2, ...)
}
\arguments{
\item{lst}{Named list of dataframes. Names are passed to \code{title =} of \code{funplot}}

\item{funplot}{Plotting function}

\item{ncols}{Number of columns if multiple plots in grid [default = 2]}

\item{...}{further options to pass to \code{funplot}}
}
\description{
Make a list of plots and arrange in a grid
}
