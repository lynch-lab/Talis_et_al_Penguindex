library(dplyr)

#
obs <- read.csv("../2 Filtering time series/prepared_all_species_obs_atleast2_1980.csv")
gepe <- read.csv("../2 Filtering time series/prepared_all_species_atleast2_1980.csv")
gepe <- gepe[which(lpi_data_all$Species_id == "GEPE"),]
years <- c(1980:2019)
site_list <- gepe$Site_ID
percent_change <- c()
not_increasing <- c()
percent_change_gent <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(percent_change_gent) <- c("Site_ID", "Population_name", "Start_year", "End_year", "Avg_p_chg")
for (i in 1:length(site_list)){
  site_id <- gepe[i, "Site_ID"]
  site_name <- gepe[i, "Site_name"]
  N <- as.numeric(gepe[i,] %>% dplyr::select(starts_with("X")))
  obs_N <- as.numeric(obs[i,] %>% dplyr::select(starts_with("X")))
  N <- N[length(N)-length(years)+1:length(years)]
  obs_N <- obs_N[length(N)-length(years)+1:length(years)]
  obs_inds <- which(!is.na(obs_N))
  N <- replace(N, which(is.na(obs_N)), NA)
  start_year <- which.min(N) + 1980 - 1
  end_year <- max(obs_inds) + 1980 - 1
  if (start_year == end_year){
    not_increasing <- c(not_increasing, site_id)
    next
  }
  if (N[max(obs_inds)] == 0){
    not_increasing <- c(not_increasing, site_id)
    next
  }
  N <- N[obs_inds]
  percent_change[i] <- (N[length(N)] - min(N))/min(N)*100
  if (min(N) == 0){
    percent_change[i] <- (N[length(N)] - min(N))*100
  }
  percent_change_gent[nrow(percent_change_gent) + 1,] <- c(site_id, site_name, start_year, end_year, percent_change[i])

}
write.csv(percent_change_gent, file = "gentoo_pop_percent_change.csv")


total_avg_percent_change_gent <- mean(percent_change, na.rm = TRUE)

save(total_avg_percent_change_gent, file = "total_avg_percent_change_gent.RData")
