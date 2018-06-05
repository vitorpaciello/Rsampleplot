#' Internal function for samplePlotR package
#'
#' Functions to calculate the parameters from the
#' sample.
#'
#' @aliases dbh2barea barea2dbh plotCircle
#' @inheritParams basefunctions
#' @param x a numeric vector. For plot should be a x cartesian data.
#' @param y a numeric vector from a cartesian y data
#' @param r a numeric vetor representing a radius of a circle
#' @param col a R color indicator
#' @param bg a R color indicator for the graphic background.
#' @param mscale a character indicating the mesured scale (mm or cm). 
#' @param statistics a function that calculates the statistics of interest from the dataframe.
#' The first argument should be the dataframe with the data and preferably should
#' @section Details:
#' 
#' Basic function to calculate descriptive statistics from a sample of a forest inventory. Basal area is calculate in 'cm' and the parameter \code{mscale} is the scale of the imput data.   
#' @return a vector of the \code{statistics}
#' 
#' @section References:
#'
#' @examples
#' 
#' \dontrun{
#' dbh2barea(x=seq(10:100, by=10), mscale= "mm")
#' }
#' @export dbh2barea barea2dbh
#' @import 
dbh2barea <- function (x, mscale = c("mm", "cm")){
    mscale <- match.arg(mscale)
    if(mscale == "mm"){
        x <- x*0.1
    }
    pi*(x/2)^2
}

barea2dbh <- function (x){    
        round((4* x/pi)^(1/2), 2)
}
############################
## Function to draw a circle
PlotSampleCircle <- Vectorize (function (x, y, r, col = 1, bg = NA) {
    ang <- seq(0,2*pi,length.out=360)
    lines(r*cos(ang)+x,r*sin(ang)+y, col = 1)
    polygon(r*cos(ang)+x, r*sin(ang)+y, col = bg)
})
