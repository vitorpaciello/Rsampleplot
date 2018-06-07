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
#' @param plotdata a data frame with a forest plot census data.
#' @return the above ground biomass of each individual
#' 
#' 
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
#' 
#' @examples
#' 
#' \dontrun{
#' bio = data.frame(dbh = rnorm(10, 10, 2), H= runif(10, 2, 8),
#'                  wd = rep(0.52,10))
#' bio = CalculateBiomass(bio, method = "chave", wd = "wd", H = "H")
#' 
#' }
#' @export CalculateBiomass
#'
CalculateBiomass <- function(plotdata, method = c("chave", "tiepolo", "burgdeli", "scatena"), 
                             wd = NULL, H = NULL, dbh = "dbh"){
  method <- match.arg(method)
  switch(method,
         "chave" = {
           plotdata$biomass <- 0.0673 * (plotdata[, wd] * plotdata[, dbh]^2 * plotdata[, H])^0.976
         },
         "tiepolo" = {
           plotdata$biomass <- 21.297 - 6.953*plotdata[, dbh] +  0.74*(plotdata[, dbh]^2)
         },
         "burgdeli" = {
           if (is.null(H)){
             plotdata$biomass <- exp(-3.068 + 2.522 * log(plotdata[, dbh]))
             } else {
               plotdata$biomass <- exp(-3.676 + 0.951 * log(plotdata[, dbh]^2 * plotdata[, H]))
               }
         },
         "scatena" = {
           if (is.null(wd)  & !is.null(H)){
              plotdata$biomass  <- exp(-3.282 + 0.95 * log(plotdata[, dbh]^2 * plotdata[,H]))
           }
           if (is.null(wd) & is.null(H)){
              plotdata$biomass <- exp(-2.399 + 2.475 * log(plotdata[, dbh]))
           }
           if(!is.null(dbh)  & !is.null(H)  & !is.null(wd)){
              plotdata$biomass <- exp(- 3.098 + 0.994 * log(plotdata[, dbh]^2 * plotdata[, H] * plotdata[, wd]))
           }
         },
         cat("Method  '", method, "'  not found\n"))
  plotdata[!is.na(plotdata$biomass), "biomass"] <- mean(plotdata[!is.na(plotdata$biomass), "biomass"])
  return(plotdata)
}