#' Random sample elements in a cartesian arena
#'
#' Draw a sample from data with elements in a cartesian coordinate system  
#' 
#'
#' @param plotdata a data frame with a forest plot census data. 
#' @param dimx numeric vector with two positions minimum and  maximum value of the x coordinate.
#' @param dimy  numeric vector with two positions minimum maximum value of the y coordinate.
#' @param xname character string with the name of the x coordinate variable.
#' @param yname character string with the name of the y coordinate variable.
#' @param nsample numeric value. The number of sample units.
#' @param shape a character string. Two shapes are possible: "rectangle" and "circle"
#' @param size numeric vector. For ''rectangle'' the dimension of x and y  sidesand for ''circle'' the radius size.
#' @param maxiter numeric value indicating the maximum number of trails to sample units without overlaping.
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
#' plotData <- singleStem(bci7, nsample = 100, shape = "circle", size = 10)
#' }
#' @export rsample
#' @import plotCircle
rsample <- function(plotdata, dimx = c(0,1000), dimy = c(0, 500), xname = "gx", yname = "gy", nsample = 100, shape = c( "rectangle", "circle"),  size = c(10,10), maxiter = 10^5){
    shape <- match.arg(shape)    
    if(shape == "rectangle")
    {
        if(length(size) != 2 | size[1] > dimx[2] | size[2] > dimy[2])
        {
            stop("Retangular sample needs two dimensions smallers then plot size")
        }
        mindist <- sqrt(sum((size/2)^2)) # this is not the minimum,but is a garantee of no overlap
    }
    if(shape == "circle")
    {
        if(length(size) != 1)
        {
            cat("Circular sample needs only the radius dimensions in size.\n
 Using only first dimension of size argument as radius \n")
            size <- size[1]    
        }
        if(size[1] >= (dimx[2]/2) | size[1] > (dimy[2]/2))
        {
            stop("Sample units should be smaller than plot size")
        }
        size <- rep(size[1]*2, 2) # from now size is diameter
        mindist <- size[1] # diameter (minimum distance)
        }
    halfs = size/2 # radius form circles and half dimensions for a rectangle
    xs <-  sample(seq(dimx[1] + halfs[1], dimx[2] - halfs[1]) ,nsample)
    ys <-  sample(seq(dimy[1] + halfs[2], dimy[2] - halfs[2]) ,nsample)
    ## above there is bias for the trees near border
    ## soluction below has a bias too and is more complicated and is not finished
    #########################################################
    ## verify if sample  is entire inside the plot dimension
    #########################################################
    xmin <- xs - halfs[1]
    ## xmin[xmin < dimx[1]] <- dimx[1]
    xmax <- xmin + size[1]
    ## outx <- which(xmax > dimx[2])
    ## soutx <- xmax[outx] - dimx[2]
    ## xmax[outx] <- xmax[outx] - soutx
    ## xmin[outx] <- xmin[outx] - soutx
    ## ## same for y axis
    ymin <- ys - halfs[2]
    ## ymin[ymin < dimy[1]] <- dimy[1]
    ymax <- ymin + size[2]
    ## outy <- which(ymax > dimy[2])
    ## souty <- ymax[outy] - dimy[2]
    ## ymax[outy] <- ymax[outy] - souty
    ## ymin[outy] <- ymin[outy] - souty
    ## check for overlap between samples
    sdist <- as.matrix( dist(cbind(xs,ys), upper=FALSE))
    sdist[upper.tri(sdist)] <- NA
    diag(sdist) <- NA
    iover <- apply(sdist, 1, function(x) {which(x < mindist,arr.ind=TRUE)})
    nover <- sapply(iover, length)
    sover <- which(nover > 0)
    niter = 1
    nsover <- length(sover)
    while(nsover > 0)
    {
        oldpar <- par(mar = c(4,4,2,2))
        plot(x= xmax, y= ymax, type= "n", xlim = dimx, ylim = dimy, xlab="x", ylab= "y")
       if(shape == "circle")
        {
            newx <- sample(seq(dimx[1] + halfs[1], dimx[2] - halfs[1]) ,nsover)
            newy <- sample(seq(dimy[1] + halfs[2], dimy[2] - halfs[2]) ,nsover)
            xs[sover] <- newx
            ys[sover] <- newy
            ## xmin[sover] <- newx - halfs[1]
            ## xmax[sover] <- newx + halfs[1]
            ## ymin[sover] <- newy - halfs[2]
            ## ymax[sover] <- newy + halfs[2]
            plotCircle(x = xs, y = ys, r = halfs[1], bg = rgb(0,1,0,0.3))
            plotCircle(x = xs[sover], y = ys[sover], r = halfs[1], bg = rgb(1,0,0,0.5))
        }
        if(shape == "rectangle")
        {
            newx <- sample(seq(dimx[1] + (size[1]/2), dimx[2] - (size[1]/2)), nsover)
            newy <- sample(seq(dimy[1] + (size[2]/2), dimy[2] - (size[2]/2)), nsover)
            xs[sover] <- newx
            ys[sover] <- newy
            xmin[sover] <- newx - halfs[1]
            xmax[sover] <- newx + halfs[1]
            ymin[sover] <- newy - halfs[2]
            ymax[sover] <- newy + halfs[2]
            rect(xleft = xmin, ybottom = ymin, xright = xmax, ytop = ymax, col = rgb(0,1,0,0.3))
            rect(xleft = xmin[sover], ybottom = ymin[sover], xright = xmax[sover], ytop = ymax[sover], col = rgb(1,0,0,0.5))
        }
        sdist <- as.matrix( dist(cbind(xs,ys), upper=FALSE))
        sdist[upper.tri(sdist)] <- NA
        diag(sdist) <- NA
        iover <- apply(sdist, 1, function(x) {which(x < mindist,arr.ind=TRUE)})
        nover <- sapply(iover, length)
        sover <- which(nover > 0)
        nsover <- length(sover)
        niter <- niter + 1
        if(niter > maxiter)
        {
            stop("It is very hard to place random samples without overlap. Please try smaller sample or increase maximum iteration (maxiter) number")
        }
        samPlots <- data.frame(x = xs, y = ys, name =  paste(ifelse(shape == "rectangle", "rect_", "circ_"), 1:nsample, sep=""), size = ifelse(shape == "rectangle",paste(size[1], size[2],sep="x"), size[1]/2))
    }
    ## select data inside sample
    sampTree <- rep(NA, nrow(plotdata))
    if(shape=="rectangle")
    {
        for(i in 1:nsample)
        {
            sampTree[plotdata[ ,xname] >= xmin[i] & plotdata[ ,xname] < xmax[i] &
                     plotdata[ ,yname] >= ymin[i] & plotdata[ ,yname] < ymax[i]] <- paste("rect_", i, sep="")
        }
    }
    if(shape=="circle")
    {
        for(i in 1:nsample)
        {
            sampTree[((plotdata[ ,xname] - xs[i])^2 +  (plotdata[ ,yname] - ys[i])^2)  < size[1]^2] <- paste("circle_", i, sep="")
        }
    }
    plotdata$sample <- sampTree
    newdata <- plotdata[!is.na(sampTree),]
    plot(x= xmax, y= ymax, type= "n", xlim = dimx, ylim = dimy, xlab="x", ylab= "y")
   #rect(xleft = xmin, ybottom = ymin, xright = xmax, ytop = ymax, col = rgb(0,1,0,0.4))
    if(shape == "rectangle") {
        rect(xleft = xmin, ybottom = ymin, xright = xmax, ytop = ymax, col = rgb(0,1,0,0.4))
    }else{
        plotCircle(x = xs, y = ys, r = size[1]/2, bg = rgb(0,1,0,0.4))
    }
    par(oldpar)
    return(list(treeData = newdata, sampleData = samPlots ))
}
