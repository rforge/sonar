% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Funktionen_Tracking.R
\name{track.eval.total}
\alias{track.eval.total}
\title{Evaluate tracking procedure}
\usage{
track.eval.total(track.evals, min.length, save.plot = FALSE, pfad.mult,
  fisch = "Fischart", a.2)
}
\arguments{
\item{track.evals}{Tracking output of \code{\link{preprocess.data}}}

\item{min.length}{Minimal length of tracks}

\item{save.plot}{\code{TRUE} (default): Plots are saved}

\item{pfad.mult}{Folder for resulting plots}

\item{fisch}{Type of fish (for main of the plots)}

\item{a.2}{Threshold for uncentered data (only for file name)}
}
\value{
Number of tracks with the minimal number of hotspots and numbers of 
 identified problems
}
\description{
This function evaluates the tracking procedure
}

