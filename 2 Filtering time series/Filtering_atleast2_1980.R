library(dplyr)

# load in prepared MAPPPD observation 0/1 csvs
# filter time series to include only those with 2 or more obs
# then load in actual MAPPPD z state csvs and filter those according to those filtered obss

lpi_data <- read.csv("prepared_all_species_obs_atleast2.csv", na.strings = "NULL")

years <- c(1980:2019)


site_list <- lpi_data$Site_ID
print(length(site_list))

atleast2 <- c()
not <- c()
for (i in 1:length(site_list)){
  N <- as.numeric(lpi_data[i,] %>% dplyr::select(starts_with("X")))[1:length(years)]
  obs_inds <- sum(N, na.rm = TRUE)
  if (obs_inds >= 2){
    atleast2 <- c(atleast2, i)
  } else{
    not <- c(not, i)
    #print(N)
  }
}

site_list <- site_list[atleast2]
print(length(site_list))


obs_df <- lpi_data[is.element(lpi_data$Site_ID, site_list),]


lpi_data <- read.csv("prepared_all_species_atleast2.csv", na.strings = "NULL")
modeled_df <- lpi_data[is.element(lpi_data$Site_ID, site_list),]



# write csv (obs)
write.csv(obs_df,"prepared_all_species_obs_atleast2_1980.csv", row.names = FALSE)

# write csv (actual z states)
write.csv(modeled_df,"prepared_all_species_atleast2_1980.csv", row.names = FALSE)





##### region stats
lpi_data_all <- read.csv("prepared_all_species_atleast2_1980.csv", na.strings = "NULL")

lpi_data_all <- lpi_data_all %>%
  mutate(Region_id = case_when(
    Region == "Central-west Antarctic Peninsula" ~ 1,
    Region == "Northwest Antarctic Peninsula" ~ 1,
    Region == "Southwest Antarctic Peninsula" ~ 2,
    Region == "Palmer Land" ~ 2,
    Region == "Elephant Island" ~ 3,
    Region == "South Orkney Islands" ~ 3,
    Region == "South Shetland Islands" ~ 3,
    Region == "Northeast Antarctic Peninsula" ~ 4,
    Ccamlr_id == "88.1" ~ 5,
    Ccamlr_id == "88.2" ~ 5,
    Ccamlr_id == "88.3" ~ 6,
    Ccamlr_id == "58.4.1" ~ 7,
    Ccamlr_id == "58.4.2" ~ 8), .after = Region) %>%   
  mutate(Ccamlr = case_when(
    Ccamlr_id == "48.1" ~ 1,
    Ccamlr_id == "48.2" ~ 1,
    Ccamlr_id == "88.1" ~ 2,
    Ccamlr_id == "88.2" ~ 2,
    Ccamlr_id == "88.3" ~ 2,
    Ccamlr_id == "58.4.1" ~ 2,
    Ccamlr_id == "58.4.2" ~ 2), .after = Ccamlr_id)

lpi_data_all_region_stats <- lpi_data_all %>%
  group_by(Region_id, Species_id) %>%
  summarise(count=n())

write.csv(lpi_data_all_region_stats, "species_region_stats_atleast2_1980.csv", row.names = FALSE)

