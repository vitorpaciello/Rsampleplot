#'  
#' 
#' Draw n samples from data with elements in a cartesian coordinate system  
#' 
#'
#' @param plotdata a data frame with a forest plot census data. 
#' @param nsamples numeric value. The number of sample units.
#' @param dimx numeric vector with two positions minimum and  maximum value of the x coordinate.
#' @param dimy  numeric vector with two positions minimum maximum value of the y coordinate.
#' @param nsims numeric value. The number of simulations.
#' @param xname character string with the name of the x coordinate variable.
#' @param yname character string with the name of the y coordinate variable.
#' @param species character string with the name of the species variable.
#' @param family character string with the name of the family variable.
#' @param biomass character string with the name of the biomass variable.
#' @param shape a character string. Two shapes are possible: "rectangle" and "circle"
#' @param size numeric vector. For ''rectangle'' the dimension of x and y  sidesand for ''circle'' the radius size.
#' @param distance numeric value. The minimum distance among plots
#' @param maxiter numeric value. indicating the maximum number of trails to sample units without overlaping.
#' 
#' @section Details:
#' 
#' Draw sample from elements in a cartesian coordinate system and return the subset of data in the sample.  
#' @return a data frame with the same columns as plotdata plus a 'sample' representing the sample unit number that the element are in. 
#' 
#' @section References:
#'
#' @examples
#' 
#' \dontrun{
#' data(bci7)
#' plotData <- RSamplePlots(bci7, nsamples = 100, shape = "circle", size = 10)
#' }
#' @export RSamplePlots
#' @import plotCircle


SimulatePlots <- function (plotdata, nsamples, dimx = NULL, dimy = NULL,
                           nsims, xname = "gx", yname = "gy", 
                           species = "species", family = "family", 
                           biomass = "biomass",
                           shape = c( "rectangle", "circle"), 
                           size = c(5,5), distance = 0,
                           maxiter = 10^5){
  
  shape <- match.arg(shape)
  
  if (is.null(dimx) | is.null(dimy)){
    dimx <- min(gx)
    dimx[2] <- max(gx)
    dimy <- min(gy)
    dimy[2] <- max(gy)
    cat("Arena dimensions not given. Using min and max coordinates of data.\n")
  }
  
  results <- data.frame (matrix(, nrow= 0 , ncol = 8))
  colnames(results) <- c("NPlots", "SampledArea", "NIndividuals", "Density", "Richness", "Shannon", "Simpson", "Biomass")
  
  for (j in 1:nsims){
    cat("\rstep", j, "of", nsims)
    for (i in nsamples:1){
      sim <- RSamplePlot(plotdata, dimx = dimx, dimy = dimy, xname = xname,
                         yname = yname, shape = shape, size = size, 
                         distance = distance, nsamples = i, maxiter = maxiter)
      tempframe <- ExtractPlotIndices(sim$treeData, sampled_area = sim$SampledArea, 
                                      family = family, species = species, biomass = biomass)
      results <- rbind(results, cbind(NPlots = i, tempframe))
    }
  }
  cat("\rSimulation completed\n")
  return(results)
}
