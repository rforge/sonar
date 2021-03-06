% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Funktionen_Preprocessing.R
\name{preprocess.data}
\alias{preprocess.data}
\title{Preprocess sonar video}
\usage{
preprocess.data(my.file, win.start, win.length, vers, y.lim, a.1, a.2, cut,
  c = 0, m.d.cs, pfad.mult, n.angle, signal.neg, do.plot, do.plot.vec,
  n.cores = 1, which.packs = NULL, frames.pack = 100, df.t = NULL,
  maxdist = 1, fastclust = FALSE, speeditup = FALSE, floodclust = FALSE,
  floodCpp = TRUE, plot.hots.mult = FALSE, only.watch = FALSE,
  schwarm.find = TRUE, regel.schwarm, which.regel.schwarm = "qda",
  save.preprocess = TRUE, save.jpg = FALSE, zeitstempel = "zeitstempel",
  save.movie = FALSE, delete.jpg = FALSE, each = 1, wait.jpg = FALSE,
  win)
}
\arguments{
\item{my.file}{Filename of .ddf file to be analyzed}

\item{win.start}{Start of the sonar window, possibly extracted before via 
\code{\link{get.version}}}

\item{win.length}{Length of the sonar window, possibly extracted before via 
\code{\link{get.version}}}

\item{vers}{Version of the .ddf file, possibly extracted before via 
\code{\link{get.version}}}

\item{y.lim}{Number of pixels to analyze in each beam, counted from the 
camera}

\item{a.1}{Threshold for centered data}

\item{a.2}{Threshold for uncentered data, default \code{0} means that no thresholding is done}

\item{cut}{Minimal size of cluster in number of pixels}

\item{c}{Number of white pixels as frames for the hotspot. Default \code{0} is a 
good choice}

\item{m.d.cs}{Maximal distance for which two hotspots can be assigned the 
same tracking number}

\item{pfad.mult}{Folder for resulting plots and output}

\item{n.angle}{Number of watch hands per quarter, i.e., total number of watch hands is \code{4 x n.angle}}

\item{signal.neg}{Is the signal negativ? Default is \code{FALSE}}

\item{do.plot}{\code{TRUE}: Plots are saved, default is \code{FALSE}}

\item{do.plot.vec}{\code{TRUE}: Plots are saved, default is \code{FALSE}}

\item{n.cores}{Number of cores if parallelization is needed, default is \code{1}}

\item{which.packs}{Default \code{NULL} results in analysis of the entire 
video. Alternatively, the indices of packages (per default 100 frames) to
be analyzed}

\item{frames.pack}{Number of frames to be analyzed in one package}

\item{df.t}{Number of degrees of freedom for the splines. Default \code{NULL} 
results in \code{c(100, 25, round(n.frames/2))} where \code{n.frames} is the 
total number of frames of the analyzed video}

\item{maxdist}{Maximal distance of two pixels to be assigned to the same cluster, 
not relevant for cluster methods floodclust and floodCpp}

\item{fastclust}{\code{TRUE}: fastclust method is used for clustering of 
pixels and definition of hotspots, default is \code{FALSE}}

\item{speeditup}{\code{TRUE}: fast version of fastclust method is used for clustering of 
pixels and definition of hotspots, default is \code{FALSE}}

\item{floodclust}{\code{TRUE}: floodclust method is used for clustering of 
pixels and definition of hotspots, default is \code{FALSE}}

\item{floodCpp}{\code{TRUE}: floodCpp method is used for clustering of 
pixels and definition of hotspots, default is \code{TRUE} - recommended method}

\item{plot.hots.mult}{\code{TRUE}: Plots are saved, default is \code{FALSE}}

\item{only.watch}{\code{TRUE}: Use watch hands instead of trapezoid variables - 
strongly recommended, default is \code{TRUE}}

\item{schwarm.find}{\code{TRUE}: A detection of shoal of fish is carried out, 
default is \code{FALSE}}

\item{regel.schwarm}{File name of classification rule for shoal of fish}

\item{which.regel.schwarm}{Type of classification rule for shoal of fish, i.e., lda, qda...}

\item{save.preprocess}{\code{TRUE}: Results of preprocessing are saved, default is \code{FALSE}}

\item{save.jpg}{\code{TRUE}: jpgs for the creating of the video are saved, 
default is \code{TRUE}}

\item{zeitstempel}{time stamp of .ddf file}

\item{save.movie}{\code{TRUE}: Video is created, 
default is \code{TRUE}}

\item{delete.jpg}{\code{TRUE}: jpgs for the creating of the video are deleted 
after creating the video, default is \code{TRUE}}

\item{each}{Interval of images for the video, for example default \code{5} 
means, that each 5th image is saved in the video}

\item{wait.jpg}{\code{TRUE}: Code pauses until creation of video is finished, 
default is \code{FALSE}}

\item{win}{Specify \code{TRUE} if you are running under windows, then, 
parallelization is not possible}
}
\value{
List with various elements necessary for further analysis steps
}
\description{
This function carries out the entire preprocessing of a sonar video, i.e.,
localizes hotspots and tracks objects over time
}

