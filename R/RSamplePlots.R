#' Random sample elements in a cartesian arena
#'
#' Draw a sample from data with elements in a cartesian coordinate system  
#' 
#'
#' @param plotdata a data frame with a forest plot census data. 
#' @param nsamples numeric value. The number of sample units.
#' @param dimx numeric vector with two positions minimum and  maximum value of the x coordinate.
#' @param dimy  numeric vector with two positions minimum maximum value of the y coordinate.
#' @param xname character string with the name of the x coordinate variable.
#' @param yname character string with the name of the y coordinate variable.
#' @param shape a character string. Two shapes are possible: "rectangle" and "circle"
#' @param size numeric vector. For ''rectangle'' the dimension of x and y sides and for ''circle'' the radius size.
#' @param distance numeric value. The minimum distance among plots
#' @param showplots Boolean. True if plot should be plotted in a graph.
#' @param maxiter numeric value. indicating the maximum number of trails to sample units without overlaping.
#' @param childfunc 
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
#' @import PlotSampleCircle
#' @import CheckPlotsOverlap

RSamplePlots <- function(plotdata, nsamples, 
                        dimx = NULL, dimy = NULL,
                        xname = "gx", yname = "gy",  
                        shape = c("rectangle", "circle"),  
                        size = NULL,  distance = 0, 
                        showplots = FALSE, maxiter = 10^5){
  
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
      size <- size[1]    
    }
    if (size[1] >= (dimx[2]/2) | size[1] > (dimy[2]/2)){
      stop("Sample units should be smaller than plot size")
    }
    size[2] <- size[1] 
    halfs <- size
    
  }
  
  xs <-  sample(seq(dimx[1] + halfs[1], dimx[2] - halfs[1]) ,nsamples)
  ys <-  sample(seq(dimy[1] + halfs[2], dimy[2] - halfs[2]) ,nsamples)
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
  
  samPlots <- data.frame(x = xs, y = ys,
                         name = paste(ifelse(shape == "rectangle", "rect_", "circ_"), 1:nsamples, sep=""),
                         size = ifelse(shape == "rectangle", paste(size[1], size[2], sep="x"), size[1]))
  sampTree <- rep(NA, nrow(plotdata))
  if (shape=="rectangle"){
    for (i in 1:nsamples){
      sampTree[plotdata[, xname] >= xs[i] - halfs[1] & plotdata[, xname] < xs[i] + halfs[1] &
                 plotdata[, yname] >= ys[i] - halfs[2] & plotdata[, yname] < ys[i] + halfs[2]] <- paste("rect_", i, sep="")
    }
    sampledarea <- size[1]*size[2]*nsamples
  }
  if (shape=="circle"){
    for (i in 1:nsamples){
      sampTree[((plotdata[, xname] - xs[i])^2 + (plotdata[, yname] - ys[i])^2) < size[1]^2] <- paste("circle_", i, sep="")
    }
    sampledarea <- size[1]*size[1]*pi*nsamples
  }
  plotdata$sample <- sampTree
  newdata <- plotdata[!is.na(sampTree),]
  
  if(showplots){
    plot(NULL, xlim = dimx, ylim = dimy, xlab="x", ylab= "y")
    if(shape == "rectangle") {
      rect(xleft = xs - halfs[1], ybottom = ys - halfs[2], xright = xs + halfs[1],
           ytop = ys + halfs[2], col = rgb(0,1,0,0.4))
    } else {
      PlotSampleCircle(x = xs, y = ys, r = size[1], bg = rgb(0,1,0,0.4))
    }
  }
  
  return(list(treeData = newdata, sampleData = samPlots, SampledArea = sampledarea))
  }
