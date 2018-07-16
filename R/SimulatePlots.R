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
#' @param size numeric vector. For ''rectangle'' the dimension of x and y sides and for ''circle'' the radius size.
#' @param distance numeric value. The minimum distance among plots
#' @param maxiter numeric value. indicating the maximum number of trails to sample units without overlaping.
#' 
#' @section Details:
#' 
#' Draw sample from elements in a cartesian coordinate system and return the subset of data in the sample.  
#' 
#' @return a data frame with the same columns as plotdata plus a 'sample' representing the sample unit number that the element are in. 
#' 
#' @section References:
#'
#' @examples
#' 
#' \dontrun{
#' data(bci7)
#' plotData <- RSamplePlots(bci7, nsamples = 100, shape = "circle", size = 10, distance =5)
#' }
#' 
#' @export RSamplePlots
#' @import CheckPlotsOverlap
#' @import ExtractSampleIndices


SimulatePlots <- function (plotdata, nsamples, dimx = NULL, dimy = NULL,
                           nsims, xname = "gx", yname = "gy", 
                           species = "species", family = "family", 
                           biomass = "biomass",
                           shape = c( "rectangle", "circle"), 
                           size = NULL, distance = 0,
                           maxiter = 10^5){
  
  shape <- match.arg(shape)
  
  if (is.null(dimx) | is.null(dimy)){
    dimx <- min(plotdata[, xname])
    dimx[2] <- max(plotdata[, xname])
    dimy <- min(plotdata[, yname])
    dimy[2] <- max(plotdata[, yname])
    cat("Arena dimensions not given. Using min and max coordinates of data.\n")
  }
  if (shape == "rectangle"){
    if (length(size) != 2 | size[1] > dimx[2] | size[2] > dimy[2]){
      stop("Retangular sample needs two dimensions smallers then plot size.")
    }
    halfs <- size/2
  }
  if (shape == "circle"){
    if (length(size) != 1){
      cat("Circular sample needs only the radius dimensions in size.
          \rUsing only first dimension of size argument as radius \n")
    }
    if (size[1] >= (dimx[2]/2) | size[1] > (dimy[2]/2)){
      stop("Sample units should be smaller than plot size")
    }
    size[2] <- size[1] 
    halfs <- size
    }
  
  results <- data.frame (matrix(, nrow= 0 , ncol = 8))
  colnames(results) <- c("NPlots", "SampledArea", "NIndividuals", "Density", "Richness", "Shannon", "Simpson", "Biomass")
  
  for (j in 1:nsims){
    cat("\rSimulation", j, "of", nsims)
    for (i in nsamples:1){
      
      plotdata$sample <- NA
      xs <-  sample(seq(dimx[1] + halfs[1], dimx[2] - halfs[1]) ,i)
      ys <-  sample(seq(dimy[1] + halfs[2], dimy[2] - halfs[2]) ,i)
      iover <- CheckPlotsOverlap(xs,ys, shape = shape, size = size, 
                                 distance = distance, childfunc = TRUE)
      nover <- sapply(iover, length)
      sover <- which(nover > 0)
      nsover <- length(sover)
      niter <- 1  
      
      while (nsover > 0){
        newx <- sample(seq(dimx[1] + halfs[1], dimx[2] - halfs[1]) ,nsover)
        newy <- sample(seq(dimy[1] + halfs[2], dimy[2] - halfs[2]) ,nsover)
        xs[sover] <- newx
        ys[sover] <- newy
        iover <- CheckPlotsOverlap(xs,ys, shape = shape, size = size, 
                                   distance = distance, childfunc = TRUE)
        iover <- lapply(seq_along(iover), function(i){ iover[[i]]<i})
        nover <- sapply(iover, sum)
        sover <- which(nover > 0)
        nsover <- length(sover)
        niter <- niter + 1
        if (niter > maxiter){
          stop ("\rIt is very hard to place random samples without overlap. 
                \rPlease try smaller sample or increase maximum iteration (maxiter) number")
        }
        }
      
      sampTree <- rep(NA, nrow(plotdata))
      if (shape=="rectangle"){
        for (k in 1:i){
          sampTree[plotdata[, xname] >= xs[k] - halfs[1] & plotdata[, xname] < xs[k] + halfs[1] &
                     plotdata[, yname] >= ys[k] - halfs[2] & plotdata[, yname] < ys[k] + halfs[2]] <- i
        }
        sampledarea <- size[1]*size[2]*i
      }
      if (shape=="circle"){
        for (k in 1:i){
          sampTree[((plotdata[, xname] - xs[k])^2 + (plotdata[, yname] - ys[k])^2) < size[1]^2] <- i
        }
        sampledarea <- size[1]*size[1]*pi*i
      }
      
      plotdata$sample <- sampTree
      newdata <- plotdata[!is.na(sampTree),]
      
      if(nrow(newdata) > 0 ){
        tempframe <-  ExtractSampleIndices(newdata, sampled_area = sampledarea, 
                                           family = family, species = species, biomass = biomass)
        results <- rbind(results, cbind(NPlots = i, tempframe))
      }
      }
  }
  cat("\rSimulation completed                    \n")
  return(results)
  }