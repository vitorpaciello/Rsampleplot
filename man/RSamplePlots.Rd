% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RSamplePlots.R
\name{RSamplePlots}
\alias{RSamplePlots}
\title{Random sample elements in a cartesian arena}
\usage{
RSamplePlots(plotdata, nsamples, dimx = NULL, dimy = NULL, xname = "gx",
  yname = "gy", shape = c("rectangle", "circle"), size = c(5, 5),
  distance = 0, showplots = FALSE, maxiter = 10^5)
}
\arguments{
\item{plotdata}{a data frame with a forest plot census data.}

\item{nsamples}{numeric value. The number of sample units.}

\item{dimx}{numeric vector with two positions minimum and  maximum value of the x coordinate.}

\item{dimy}{numeric vector with two positions minimum maximum value of the y coordinate.}

\item{xname}{character string with the name of the x coordinate variable.}

\item{yname}{character string with the name of the y coordinate variable.}

\item{shape}{a character string. Two shapes are possible: "rectangle" and "circle"}

\item{size}{numeric vector. For ''rectangle'' the dimension of x and y  sidesand for ''circle'' the radius size.}

\item{distance}{numeric value. The minimum distance among plots}

\item{showplots}{Boolean. True if plot should be plotted in a graph.}

\item{maxiter}{numeric value. indicating the maximum number of trails to sample units without overlaping.}
}
\value{
a data frame with the same columns as plotdata plus a 'sample' representing the sample unit number that the element are in.
}
\description{
Draw a sample from data with elements in a cartesian coordinate system
}
\section{Details}{


Draw sample from elements in a cartesian coordinate system and return the subset of data in the sample.
}

\section{References}{

}

\examples{

\dontrun{
data(bci7)
plotData <- RSamplePlots(bci7, nsamples = 100, shape = "circle", size = 10)
}
}
