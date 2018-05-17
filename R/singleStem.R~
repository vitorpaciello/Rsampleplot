#' Calculate dbh for a tree with multiple stems
#'
#' Functions to estimate a single dbh from a multiple stem tree or shrub.
#' 
#'
#' @inheritParams basefunctions
#' @param  plotdata a data frame with a forest plot census data. 
#' @param mscale a character indicating the mesured scale (mm or cm).
#' @param tag a character string with the name of the variable unique for each tree.
#' @param dbh a character string with the name of the variable with the diameter of the trees
#' @section Details:
#' 
#' Calculate the sum of basal areas for multiple stem tree and them return the dbh associated with the sum of basala area.   
#' @return a data frame with the same columns as plotdata plus a 'dbh.cm' column. 
#' 
#' @section References:
#'
#' @examples
#' 
#' \dontrun{
#' data(bci7)
#' plotData <- singleStem(bci7, mscale= "mm")
#' }
#' @export singleStem
#' @import dbh2barea barea2dbh
singleStem <-
function(plotdata, tag = "tag", dbh = "dbh",  mscale = "mm"){
    #remove record with NA in dbh field 
    pd <- plotdata[!(is.na(plotdata[,dbh]) | is.na(plotdata[,tag])),]
    #subset the data (usually status or census field)
    ## check if there is multiple stem trees
    if(sum(duplicated(pd[,tag])) == 0)
    {
        cat("The data frame has no duplicated tags. Go ahead!\n")
        stop("There is no multiple stem plants. You do not need to run this function")
    }
    new_namedbh <- paste(dbh,".cm", sep="")
    basalarea <- dbh2barea(pd[,dbh], mscale)
    bareaTree <- aggregate(basalarea, list(tag = pd[,tag]), FUN = sum)
    data.cm <- pd[!duplicated(pd[, tag]), ]
    mtag <- match(data.cm[,tag], bareaTree[,tag])
    data.cm[,new_namedbh] <- bareaTree[mtag, "x"]
    data.cm <- data.cm[, !(names(data.cm) %in% c(dbh, "fuste", "stem", "stemID")) ]
    return(data.cm)
}
