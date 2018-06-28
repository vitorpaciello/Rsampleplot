#' @title Checks if plots in a cartesian arena overlap
#' 
#' @param x_coords vector of x coordinates of the plots center.
#' @param y_coords vector of x coordinates of the plots center.
#' @param angles vector of the orientation angles of the plots center.
#' @param shape a character string. Two shapes are possible: "rectangle" and "circle"
#' @param size numeric vector. For ''rectangle'' the dimension of x and y sides and for ''circle'' the radius size.
#' @param distance numeric value. The minimum distance among plots
#'
#'
#' @details 
#' 
#' @return A list of the plots with their overlapping plots
#' 
#' @examples 
#' \dontrun{
#' 
#' x <- sample(seq(0, 1000), 50)
#' y <- sample(seq(0, 1000) , 50)
#' radius <- 10
#' overlap <- CheckPlotsOverlap(x, y, shape = "circle", size = radius, distance = 10)
#' }
#' 
#' @export CheckPlotsOverlap

CheckPlotsOverlap <- function(x_coords, y_coords, angles = NULL, 
                              shape = c("circle", "rectangle"),
                              size = c(10,10), distance = 0,
                              childfunc = FALSE){
  
  
  if (!childfunc){
    shape <- match.arg(shape)
    if (length(x_coords) != length(y_coords)){
      stop("Coordinate vectors should have same length.")
    }
    if (!is.null(angles)){
      cat("Angles functionality under development. All plots are assumed to have the same orientation.\n")
    }
    if (shape == "circle"){
      if (length(size) != 1){
        cat("Circular plots need only the radius dimensions in size.
            \rUsing only first dimension of size argument as radius \n")
      }
    }
    if (shape == "rectangle"){
      if (length(size) != 2 | size[1] < 0 | size[2] < 0){
        stop("Retangular plots need two positive dimensions.")
      }
    }
  }
  if (shape == "circle"){
    mindist <- 2*size[1] + distance
    sdist <- as.matrix( dist(cbind(x_coords,y_coords), upper=FALSE))
    sdist[upper.tri(sdist)] <- t(sdist)[upper.tri(sdist)]
    diag(sdist) <- mindist +1
    output <- apply(sdist, 1, function(x){
                    which(x < mindist,arr.ind=TRUE)})
    return(output)
  }
  if (shape == "rectangle"){
    mindist_x <- size[1] + distance
    mindist_y <- size[2] + distance
    mxcoords <- matrix(x_coords, nrow = length(x_coords), ncol = length(x_coords))
    mycoords <- matrix(y_coords, nrow = length(y_coords), ncol = length(y_coords))
    mdistx <- abs(mxcoords-t(mxcoords))-mindist_x
    mdisty <- abs(mycoords-t(mycoords))-mindist_y
    comparedist <- mdistx * mdisty
    comparedist[comparedist <= 0] <- 0
    diag(comparedist) <- 0
    comparedist[mdistx > 0] <- 0
    output <- apply(comparedist, 1, function(x){
                    which(x > 0,arr.ind=TRUE)})
  }
  return (output)
}
