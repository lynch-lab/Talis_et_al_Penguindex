library(dplyr)


# load in prepared MAPPPD observation 0/1 csvs
# filter time series to include only those with 2 or more obs
# then load in actual MAPPPD z state csvs and filter those according to those filtered obss


# adelies
lpi_data_ad <- read.csv("../1 BSSM outputs/prepared_ADPE_obs.csv", na.strings = "NULL")

site_list_ad <- lpi_data_ad$Site_ID
print(length(site_list_ad))

atleast2 <- c()
#not <- c()
for (i in 1:length(site_list_ad)){
  N <- as.numeric(lpi_data_ad[i,] %>% dplyr::select(starts_with("X")))
  obs_inds <- sum(N, na.rm = TRUE)
  if (obs_inds >= 2){
    atleast2 <- c(atleast2, i)
  } #else{
    #not <- c(not, i)
    #print(N)
  #}
}
site_list_ad <- site_list_ad[atleast2]
print(length(site_list_ad))

adelie_obs_df <- lpi_data_ad[is.element(lpi_data_ad$Site_ID, site_list_ad),]


lpi_data_ad <- read.csv("../1 BSSM outputs/prepared_ADPE.csv", na.strings = "NULL")
adelie_modeled_df <- lpi_data_ad[is.element(lpi_data_ad$Site_ID, site_list_ad),]


# chinstraps
lpi_data_ch <- read.csv("../1 BSSM outputs/prepared_CHPE_obs.csv", na.strings = "NULL")

site_list_ch <- lpi_data_ch$Site_ID
print(length(site_list_ch))

atleast2 <- c()
for (i in 1:length(site_list_ch)){
  N <- as.numeric(lpi_data_ch[i,] %>% dplyr::select(starts_with("X")))
  obs_inds <- sum(N, na.rm = TRUE)
  if (obs_inds >= 2){
    atleast2 <- c(atleast2, i)
  }
}
site_list_ch <- site_list_ch[atleast2]
print(length(site_list_ch))

chinstrap_obs_df <- lpi_data_ch[is.element(lpi_data_ch$Site_ID, site_list_ch),]


lpi_data_ch <- read.csv("../1 BSSM outputs/prepared_CHPE.csv", na.strings = "NULL")
chinstrap_modeled_df <- lpi_data_ch[is.element(lpi_data_ch$Site_ID, site_list_ch),]




# gentoos
lpi_data_ge <- read.csv("../1 BSSM outputs/prepared_GEPE_obs.csv", na.strings = "NULL")

site_list_ge <- lpi_data_ge$Site_ID
print(length(site_list_ge))

atleast2 <- c()
for (i in 1:length(site_list_ge)){
  N <- as.numeric(lpi_data_ge[i,] %>% dplyr::select(starts_with("X")))
  obs_inds <- sum(N, na.rm = TRUE)
  if (obs_inds >= 2){
    atleast2 <- c(atleast2, i)
  }
}
site_list_ge <- site_list_ge[atleast2]
print(length(site_list_ge))

gentoo_obs_df <- lpi_data_ge[is.element(lpi_data_ge$Site_ID, site_list_ge),]


lpi_data_ge <- read.csv("../1 BSSM outputs/prepared_GEPE.csv", na.strings = "NULL")
gentoo_modeled_df <- lpi_data_ge[is.element(lpi_data_ge$Site_ID, site_list_ge),]



# bind all 3 species together (obs)
all_species_obs_df <- rbind(adelie_obs_df, chinstrap_obs_df, gentoo_obs_df)
# write csvs (obs)
write.csv(adelie_obs_df,"prepared_ADPE_obs_atleast2.csv", row.names = FALSE)
write.csv(chinstrap_obs_df,"prepared_CHPE_obs_atleast2.csv", row.names = FALSE)
write.csv(gentoo_obs_df,"prepared_GEPE_obs_atleast2.csv", row.names = FALSE)
write.csv(all_species_obs_df,"prepared_all_species_obs_atleast2.csv", row.names = FALSE)




# bind all 3 species together (actual z states)
all_species_modeled_df <- rbind(adelie_modeled_df, chinstrap_modeled_df, gentoo_modeled_df)
# write csvs (actual z states)
write.csv(adelie_modeled_df,"prepared_ADPE_atleast2.csv", row.names = FALSE)
write.csv(chinstrap_modeled_df,"prepared_CHPE_atleast2.csv", row.names = FALSE)
write.csv(gentoo_modeled_df,"prepared_GEPE_atleast2.csv", row.names = FALSE)
write.csv(all_species_modeled_df,"prepared_all_species_atleast2.csv", row.names = FALSE)



