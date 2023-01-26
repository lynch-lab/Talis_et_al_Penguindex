library(mapppdr) # brings in sites df, among others
# see: https://github.com/CCheCastaldo/mapppdr

# script to load in MAPPPD z states from the .rda files from Chris
#  & format them into better csv file

# also, write csvs for each species with y/n observation data for each time series (to run stats on data dens)

# function to convert MAPPPD abundance_lambdas to better dataframe format
convert_abundance_lambda <- function(species_id, abundance_lambda){
  years <- c(1970:2019)
  
  if (species_id == "ADPE"){
    binomial = "Pygoscelis_adeliae"
  } else if (species_id == "CHPE"){
    binomial = "Pygoscelis_antarcticus"
  } else if (species_id == "GEPE"){
    binomial = "Pygoscelis_papua"
  }
  
  x <- c("ID", "Binomial",
         "Species_id", "Site_ID", "Site_name",
         "Latitude", "Longitude",
         "Region", "Ccamlr_id",
         sapply(years, toString, how = "replace"))
  species_modeled_df <- data.frame(matrix(ncol = length(x), nrow = 0))
  colnames(species_modeled_df) <- x
  
  species_obs_df <- data.frame(matrix(ncol = length(x), nrow = 0))
  colnames(species_obs_df) <- x
  
  site_list <- unique(abundance_lambda$site_id)
  
  for (i in 1:length(site_list)){
    site_id <- site_list[i]
    inds <- which(abundance_lambda$site_id == site_id)
    site_name <- abundance_lambda$site_name[inds[1]]
    lat <- abundance_lambda$latitude[inds[1]]
    long <- abundance_lambda$longitude[inds[1]]
    region <- sites$region[which(sites$site_id == site_id)]
    if (length(region) == 0) {
      if (site_id == "AMBU"){
        region <- "Northeast Antarctic Peninsula"
      } else if (site_id == "PTWD" || site_id == "CRBA"){
        region <- "Elephant Island"
      } else if (site_id == "SAND"){
        region <- "South Orkney Islands"
      } else{
        print(paste(site_id, site_name))
        region <- NA
      }
    }
    ccamlr_id <- abundance_lambda$ccamlr_id[inds[1]]
    
    modeled_abundance <- c()
    obs <- c()
    for (yr in 1:length(years)){
      obs[yr] <- abundance_lambda$obs_1[(abundance_lambda$site_id == site_id) & (abundance_lambda$season == years[yr])]
      if (obs[yr] == 0){
        obs[yr] <- NA
      }
      modeled_abundance_raw <- abundance_lambda$lza_mean[(abundance_lambda$site_id == site_id) & (abundance_lambda$season == years[yr])]
      if (is.na(modeled_abundance_raw)){
        modeled_abundance[yr] <- 0
      } else{
        modeled_abundance[yr] <- exp(modeled_abundance_raw)
      }
    }
    
    species_modeled_df[nrow(species_modeled_df) + 1,] = c(i, binomial, species_id, site_id, site_name,
                                                        lat, long, region, ccamlr_id,
                                                        modeled_abundance)
    species_obs_df[nrow(species_obs_df) + 1,] = c(i, binomial, species_id, site_id, site_name,
                                                          lat, long, region, ccamlr_id,
                                                          obs)
    
  }
  return(list(species_modeled_df, species_obs_df))
}


# adelies
load("ADPE_abundance_lambda.rda")
convert_abundance_lambda_output <- convert_abundance_lambda("ADPE", ADPE_abundance_lambda)
adelie_modeled_df <- convert_abundance_lambda_output[[1]]
adelie_obs_df <- convert_abundance_lambda_output[[2]]

# chinstraps
load("CHPE_abundance_lambda.rda")
convert_abundance_lambda_output <- convert_abundance_lambda("CHPE", CHPE_abundance_lambda)
chinstrap_modeled_df <- convert_abundance_lambda_output[[1]]
chinstrap_obs_df <- convert_abundance_lambda_output[[2]]

# gentoos
load("GEPE_abundance_lambda.rda")
convert_abundance_lambda_output <- convert_abundance_lambda("GEPE", GEPE_abundance_lambda)
gentoo_modeled_df <- convert_abundance_lambda_output[[1]]
gentoo_obs_df <- convert_abundance_lambda_output[[2]]


# bind all 3 species together
all_species_modeled_df <- rbind(adelie_modeled_df, chinstrap_modeled_df, gentoo_modeled_df)


# write csvs for each species and for all 3
write.csv(adelie_modeled_df,"prepared_ADPE.csv", row.names = FALSE)
write.csv(chinstrap_modeled_df,"prepared_CHPE.csv", row.names = FALSE)
write.csv(gentoo_modeled_df,"prepared_GEPE.csv", row.names = FALSE)
write.csv(all_species_modeled_df,"prepared_all_species.csv", row.names = FALSE)



# bind all 3 species together (obs)
all_species_obs_df <- rbind(adelie_obs_df, chinstrap_obs_df, gentoo_obs_df)

# write csvs for each species' observation data
write.csv(adelie_obs_df,"prepared_ADPE_obs.csv", row.names = FALSE)
write.csv(chinstrap_obs_df,"prepared_CHPE_obs.csv", row.names = FALSE)
write.csv(gentoo_obs_df,"prepared_GEPE_obs.csv", row.names = FALSE)
write.csv(all_species_obs_df,"prepared_all_species_obs.csv", row.names = FALSE)



