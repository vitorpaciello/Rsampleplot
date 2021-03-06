% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulatePlots.R
\name{SimulatePlots}
\alias{SimulatePlots}
\title{Draw n samples from data with elements in a cartesian coordinate system}
\usage{
SimulatePlots(plotdata, nsamples, dimx = NULL, dimy = NULL, nsims,
  xname = "gx", yname = "gy", species = "species", family = "family",
  biomass = "biomass", shape = c("rectangle", "circle"), size = c(5, 5),
  distance = 0, maxiter = 10^5)
}
\arguments{
\item{plotdata}{a data frame with a forest plot census data.}

\item{nsamples}{numeric value. The number of sample units.}

\item{dimx}{numeric vector with two positions minimum and  maximum value of the x coordinate.}

\item{dimy}{numeric vector with two positions minimum maximum value of the y coordinate.}

\item{nsims}{numeric value. The number of simulations.}

\item{xname}{character string with the name of the x coordinate variable.}

\item{yname}{character string with the name of the y coordinate variable.}

\item{species}{character string with the name of the species variable.}

\item{family}{character string with the name of the family variable.}

\item{biomass}{character string with the name of the biomass variable.}

\item{shape}{a character string. Two shapes are possible: "rectangle" and "circle"}

\item{size}{numeric vector. For ''rectangle'' the dimension of x and y  sidesand for ''circle'' the radius size.}

\item{distance}{numeric value. The minimum distance among plots}

\item{maxiter}{numeric value. indicating the maximum number of trails to sample units without overlaping.}
}
\value{
a data frame with the same columns as plotdata plus a 'sample' representing the sample unit number that the element are in.
}
\description{
Draw n samples from data with elements in a cartesian coordinate system
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
