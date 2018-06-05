#' @title Extract biodiversity indices from sampled plots
#' 
#' @param data dataframe with a forest plot census data.
#' @param sampled_area numeric value. Total sampled area.
#' @param species character string with the name of the species variable.
#' @param family character string with the name of the family variable.
#' @param biomass character string with the name of the biomass variable.
#' @param sample character string with the name of the sample variable got from the rsample function.
#'
#'
#' @details 
#' 
#' @return A data frame with the calculated biodiversity indices
#' 
#' @examples 
#' \dontrun{
#' 
#' 
#' 
#' 
#' 
#' }
#' 
#' @export ExtractSampleIndices


ExtractSampleIndices<- function(data, sampled_area = NA,
                              species = "species", family  = "family",
                              biomass = "biomass", sample = "sample"){
  
  results<- data.frame(sampled_area = numeric(0))
  results[1,1]<- sampled_area
  results$n_individuals<- nrow(data)
  results$density <- results[1,2]/sampled_area
  results$richness <- length(unique(data[, species][!is.na(data[, species])]))
  species_proportions <- aggregate((count = data[, species]), list(value = data[, species]), length)[2]
  species_proportions <- species_proportions/ sum(species_proportions)
  results$shannon <- (-1) * sum(species_proportions * log(species_proportions))
  results$simpson <- 1 - sum(species_proportions^2)
  if (biomass %in% colnames(data)){
    results$biomass <- sum(data[, biomass])
  } else {
    results$biomass <- NA
  }
  return (results)
}
