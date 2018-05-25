#' @title Above ground biomas calculation
#' @description
#' There several alometric equations to use. Possibilities are:
#' \code{tiepolo}
#' Tiepolo et al. (2002) 
#' $$ AGB = 21.297 - 6.953 dbh +  0.740 (dbh)^2 $$
#' where:
#' \code{dbh}  diameter at breast high (cm)
#'
#' \code{scatena}
#' Scatena et al. 1993
#' DBH only:
#' $$ AGB = exp(-2.399 + 2.475 x ln(dbh)) $$
#' DBH and Plant Height
#' $$ AGB = exp(-3.228 + 0.95 x ln(dbh^2 x H)) $$
#' DBH, Height and WD
#' $$ AGB = exp(- 3.098 + 0.994 x ln(dbh^2 x H x wd)) $$
#' Chave 2014
#' \code{chave}
#' In this pan-tropical allometric model, the specific wood density is required.
#' $$ AGB_{est} = 0.0673 x (wd x dbh^2 x H)^0.976 $$
#'
#' \code{burgdeli}
#' Burger e Delitti 2008
#' $$ AGB_{est} = exp(-3.068 + 2.522 * log(dbh)) $$
#' $$ AGB_{est} = exp(-3.676  + 0.951 x ln(dbh^2 x H)) $$
#' 
#' @aliases chave scatena tiepolo burgdeli
#' @param dbh numeric vector with diameter at breast height in cm
#' @param H numeric vector with plants height
#' @param wd numeric vector with wood density
#' @return the above ground biomass
#' @author Alexandre Adalardo de Oliveira and Paulo Inacio Prado
#' \email{aleadalardo@@gmail.com}
#'
#' @references
#' SCATENA, F.N., SILVER, W., SICCAMA, T., JOHNSON, A. & SÁNCHEZ, M.J.  1993.
#' Biomass  and  nutrient  content  of  the  Bisley  Experimental Watersheds,
#' Luquillo Experimental Forest, Puerto Rico, before and after Hurricane Hugo, 1989.
#' Biotropica 25(1):15-27.
#' 
#' TIEPOLO,  G.,  CALMON,  M.  &  FERETTI,  A.R.  2002.  Measuring  and Monitoring
#' Carbon  Stocks  at  the  Guaraqueçaba  Climate  Action Project, Paraná, Brazil.
#' In: International Symposium on Forest Carbon Sequestration and Monitoring.
#' Extension Serie Taiwan Forestry Research Institute 153:98-115.
#'
#' BURGER   D.M.;   DELITTI,   W.B.C.   allometric   models   for   estimating   the
#' phytomass of a secondary  Atlantic  Forest  area  of  southeastern  Brazil. 
#' Biota  Neotrop.8 (4)
#'
#' @seealso \code{\link{bootree}}
#' 
#' @examples
#' 
#' \dontrun{
#' bio = data.frame(dbh = rnorm(10, 10, 2), H= runif(10, 2, 8),
#'                  wd = rep(0.52,10), plots = rep(c(1,2),5))
#' wich(bio, scatena(dbh, H, wd))
#' wich(bio, chave(dbh, H, wd))
#' wich(bio, burgdeli(dbh, H, wd))
#' }
#' @export chave burgdeli scatena tiopolo
#'
#############################################
tiepolo <- function(dbh, H = NULL, wd = NULL)
{
    21.297 - 6.953*dbh +  0.74*(dbh^2)
}
scatena <- function(dbh, H = NULL , wd = NULL)
{
    if(is.null(wd)  & !is.null(H))
    {
       agb <- exp(-3.282 + 0.95 * log(dbh^2 * H))
    }
    if(is.null(wd) & is.null(H))
    {
        agb <- exp(-2.399 + 2.475 * log(dbh))
    }
    if(!is.null(dbh)  & !is.null(H)  & !is.null(wd))
    {
        agb <- exp(- 3.098 + 0.994 * log(dbh^2 * H * wd))
    }
    agb
}
#############################
chave <- function(dbh, H, wd)
{
    0.0673 * (wd * dbh^2 * H)^0.976
}

burgdeli <- function(dbh, H= NULL, wd = NULL)
{
    if(is.null(H))
    {
        exp(-3.068 + 2.522 * log(dbh))
    }else{

        exp(-3.676 + 0.951 * log(dbh^2 * H))
    }
}
##################################################################
#' @title Bootstrap Confidence Intervals for tree biomass samples
#' @description
#' \code{bootree} return confidence interval for biomass (AGB) from samples of trees in a plot
#' 
#' @param plotdata a data frame with a forest sample plots data. 
#' @param samplename character vector with the name of the column in plotdata that indicate the plots codes.
#' @param area numeric, plot size in square meters. 
#' @param varnames vector of names variables variable diameter at breast hight, plant height and wood density, note that they should be  named as 'dbh', 'H' and 'wd' respectivelly. Put only those varibles that are in plotdata and needed for alometric calculation.
#' @param alom alometric function to be used, possibilities are \strong{chave} \strong{burgdeli}, \strong{scatena} and \strong{tiopolo}. Those are functions so should be write without quotes.
#' @param quantil Bootstrap quantiles to return, default is 95\% limite confidence interval (0.25 and 0.975 quatiles);
#' @param nsim numeric number of resamples from data to calculate bootstrap.
#' @param convtoMgha logical, if true convert AGB to Mg (tons) per hectar, otherwise the variable is in Kg per sample size.
#' @author Alexandre Adalardo de Oliveira and Paulo Inacio Prado
#' \email{aleadalardo@@gmail.com}
#' @return biomass calculation and bootstrap confidence interval, as well as standard deviation for the simulated bootstrap distribution 
#'
#' @details
#' Bootstrap resamples above ground biomass of each sample, reordering the sample tag.
#' 
#' @references
#' SCATENA, F.N., SILVER, W., SICCAMA, T., JOHNSON, A. & SÁNCHEZ, M.J.  1993.
#' Biomass  and  nutrient  content  of  the  Bisley  Experimental Watersheds,
#' Luquillo Experimental Forest, Puerto Rico, before and after Hurricane Hugo, 1989.
#' Biotropica 25(1):15-27.
#' 
#' TIEPOLO,  G.,  CALMON,  M.  &  FERETTI,  A.R.  2002.  Measuring  and Monitoring
#' Carbon  Stocks  at  the  Guaraqueçaba  Climate  Action Project, Paraná, Brazil.
#' In: International Symposium on Forest Carbon Sequestration and Monitoring.
#' Extension Serie Taiwan Forestry Research Institute 153:98-115.
#'
#' BURGER   D.M.;   DELITTI,   W.B.C.   allometric   models   for   estimating   the
#' phytomass of a secondary  Atlantic  Forest  area  of  southeastern  Brazil. 
#' Biota  Neotrop.8 (4)
#'
#' @seealso \code{\link{chave}} \code{\link{scatena}} \code{\link{tiepolo}} \code{\link{burgdeli}}
#' 
#' @examples
#' 
#' \dontrun{
#' bio = data.frame(dbh = rnorm(10, 10, 2), H= runif(10, 2, 8),
#'                  wd = rep(0.52,10), plots = rep(c(1,2),5))
#' bootree(bio, samplename = "plot", varnames =c(dbh="dbh", H = "H", wd = "wd"), nsim = 10)
#' metaEr(cl=10, rw=10, f0=0.2, pi=0.2, ce=0.15, tmax=100)
#' metaCiEr(cl=10, rw=10, f0=0.2, ci=0.2, ce=0.15, tmax=100)
#' }
#' @export bootree
#' @import chave burgdeli scatena tiopolo
#' @export bootree
bootree <- function(plotdata, samplename = "fragpar", area = 100,  varnames = c(dbh = "diameter.breast.height", H = "plant.height", wd = "wood.density"), alom = chave, quantil = c(0.025, 0.975), nsim=1000, convtoMgha = TRUE)
{
    if(sum(!(varnames %in% names(plotdata)))>0)
    {
        stop("Variables (varnames) are not in the data frame")
    }
    bio <- as.data.frame(plotdata[, varnames])
    names(bio) <- names(varnames)
    agb <- with(bio, alom(dbh, H, wd))
    agbaggr <- tapply(agb, INDEX = plotdata[,samplename], sum)
    codpar <- names(agbaggr)
    bootagb <- matrix(NA, nrow = length(agbaggr), ncol = nsim)
    bootagb[,1] <- agbaggr
    for(i in 2:nsim) 
    {
        sagbplot <- sapply(tapply(agb, INDEX = plotdata[,samplename], FUN = sample, replace=TRUE), sum)
        bootagb[,i] <- sagbplot[codpar]
    }
    agbci <- t(apply(bootagb, 1, quantile, probs= quantil))
    agbsd <- apply(bootagb, 1, sd)
    agball <- cbind(agbaggr, agbsd, agbci)
    if(convtoMgha)
    {
        convha <- 10/area
        agball <-  agball* convha
    }
    return(agball)
}
