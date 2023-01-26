library(dplyr)
library(ggplot2)

# Global - species - region
# GLOBAL penguin LPI / Penguindex -- for each species, calc lambdas for each region (weight by # of pops), aggregate over spp

# input: MAPPPD output full **1980**-2019 TIME SERIES WITH AT LEAST 2 OBSERVATIONS
# 0 pop values (which are, according to MAPPPD model, observed no occupancy at end/beginning of time series)
# are changed to 1% of the mean pop count for that time series
# output: penguindex null models

start_time <- Sys.time()

lpi_data_all <- read.csv("../2 Filtering time series/prepared_all_species_atleast2_1980.csv", na.strings = "NULL")

species_list <- unique(lpi_data_all$Species_id) # "ADPE", "CHPE", "GEPE"

years <- c(1980:2019) ###############

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
region_list <- sort(unique(lpi_data_all$Region_id)) # 1-8
lpi_data_all_region_stats <- lpi_data_all %>%
  group_by(Species_id, Region_id) %>%
  summarise(count=n())

species_regions <- list()
species_regions_site_list_lengths <- list()
# species_regions is a list of vectors, one for each species, of the actual Region_ids that have that species
for (sp_i in 1:length(species_list)){
  species <- species_list[sp_i]
  vec <- lpi_data_all_region_stats$Region_id[which(lpi_data_all_region_stats$Species_id == species)]
  species_regions[[sp_i]] <- vec
  vec <- lpi_data_all_region_stats$count[which(lpi_data_all_region_stats$Species_id == species)]
  species_regions_site_list_lengths[[sp_i]] <- vec
}

species_total_num_pops <- c()
for (sp_i in 1:length(species_list)){
  species <- species_list[sp_i]
  species_total_num_pops[sp_i] <- length(lpi_data_all$Site_ID[which(lpi_data_all$Species_id == species)])
}


######## process error
# Adelie in CCAMLAR 48.1/48.2 --------- 0.262
# Adelie in CCAMLAR 88.1/88.2/88.3 ---- 0.169
# Adelie in CCAMLAR 58.4.1/58.4.2 ----- 0.169
# Chinstrap --------------------------- 0.329
# Gentoo ------------------------------ 0.293
sigma_proc_adelie <- c(rep(0.262, 4), rep(0.169, 4))
sigma_proc_chinstrap <- c(rep(0.329, 8))
sigma_proc_gentoo <- c(rep(0.293, 8))
#
sigma_proc_reg <- list(sigma_proc_adelie, sigma_proc_chinstrap, sigma_proc_gentoo)


############
#
# null model
#
#
#
##############
#
#
# load in zstates with posterior mean and sd
#
adelie_z_df <- read.csv("../1 BSSM outputs/ADPE_abundance_lambda.csv", na.strings = "NULL")
chinstrap_z_df <- read.csv("../1 BSSM outputs/CHPE_abundance_lambda.csv", na.strings = "NULL")
gentoo_z_df <- read.csv("../1 BSSM outputs/GEPE_abundance_lambda.csv", na.strings = "NULL")
all_species_z_df <- rbind(adelie_z_df, chinstrap_z_df, gentoo_z_df)
#
#
iter_num <- 1000 ######
nullmodel_d <- array(NA, dim = c(length(species_list), length(region_list), max(unlist(species_regions_site_list_lengths)), length(years)))
nullmodel_I_speciesregionlevel <- array(NA, dim = c(length(species_list), length(region_list), iter_num, length(years)))
nullmodel_I_specieslevel <- array(NA, dim = c(length(species_list), iter_num, length(years)))
nullmodel_I_global <- matrix(NA, nrow = iter_num, ncol = length(years))
for (iter in 1:iter_num){
  #
  nullmodel_d_bar_region <- array(NA, dim = c(length(species_list), length(region_list), length(years)))
  nullmodel_d_bar_species <- matrix(NA, nrow = length(species_list), ncol = length(years))
  nullmodel_I_global[iter,1] <- 1
  #################
  for (sp_i in 1:length(species_list)){
    species <- species_list[sp_i]
    lpi_data_species <- lpi_data_all[which(lpi_data_all$Species_id == species),]
    region_weights <- c()
    for (reg in species_regions[[sp_i]]){
      reg_i_rel <- which(species_regions[[sp_i]] == reg)
      lpi_data <- lpi_data_species[which(lpi_data_species$Region_id == reg),]
      site_list <- lpi_data$Site_ID # number of sites in region for that species
      # process error
      sigma <- sigma_proc_reg[[sp_i]][reg]
      for (i in 1:length(site_list)){
        site_id <- lpi_data[i, "Site_ID"]
        N <- as.numeric(lpi_data[i,] %>% dplyr::select(starts_with("X")))
        N <- N[length(N)-length(years)+1:length(years)]
        firstcol <- which(colnames(lpi_data) == "X1980") # first col
        N_sim <- rep(NA, times = length(N))
        # start with initial abundance estimate (1980)
        mean_i1 <- as.numeric(all_species_z_df$lza_mean[(all_species_z_df$species_id == species) & (all_species_z_df$site_id == site_id) & (all_species_z_df$season == years[1])])
        sd_i1 <- as.numeric(all_species_z_df$lza_sd[(all_species_z_df$species_id == species) & (all_species_z_df$site_id == site_id) & (all_species_z_df$season == years[1])])
        N_sim[1] <- rlnorm(1, meanlog = mean_i1, sdlog = sd_i1)
        for (t in 2:length(years)){
          N_sim[t] <- exp(rnorm(1, log(N_sim[t-1]), sigma[sp_i]))
          nullmodel_d[sp_i,reg,i,t] <- log10(N_sim[t]/N_sim[t-1])
        }
      }
      # now for each region, calculate d_bar(spp,reg,t) & region-level index
      nullmodel_I_speciesregionlevel[sp_i,reg,iter,1] <- 1
      for (t in 2:length(years)){
        nullmodel_d_bar_region[sp_i,reg,t] <- (1/species_regions_site_list_lengths[[sp_i]][reg_i_rel])*sum(nullmodel_d[sp_i,reg,,t], na.rm = TRUE) # over each site
        nullmodel_I_speciesregionlevel[sp_i,reg,iter,t] <- nullmodel_I_speciesregionlevel[sp_i,reg,iter,t-1]*10^(nullmodel_d_bar_region[sp_i,reg,t])
      }
      # region weights
      region_weights[reg] <- species_regions_site_list_lengths[[sp_i]][reg_i_rel]/species_total_num_pops[sp_i]
    }
    region_weights <- c(region_weights, rep(NA, times=length(region_list)-length(region_weights)))
    # now for each spp, calculate d_bar(spp,t) & species-level index
    nullmodel_I_specieslevel[sp_i,iter,1] <- 1
    for (t in 2:length(years)){
      nullmodel_d_bar_species[sp_i,t] <- (1/length(species_regions[[sp_i]]))*sum(nullmodel_d_bar_region[sp_i,,t]*region_weights, na.rm = TRUE)
      nullmodel_I_specieslevel[sp_i,iter,t] <- nullmodel_I_specieslevel[sp_i,iter,t-1]*10^(nullmodel_d_bar_species[sp_i,t])
    }
  }
  # calculate global index
  for (t in 2:length(years)){
    nullmodel_I_global[iter,t] <- nullmodel_I_global[iter,t-1]*10^(mean(nullmodel_d_bar_species[,t]))
  }
}

end_time <- Sys.time()

end_time-start_time


###############################################################
# SAVE
save(lpi_data_all, species_regions_site_list_lengths, nullmodel_d, nullmodel_d_bar_region, nullmodel_d_bar_species, nullmodel_I_speciesregionlevel, nullmodel_I_specieslevel, nullmodel_I_global, file = "nullmodel_data.RData")


###############################################################
# LOAD
load("nullmodel_data.RData")

###############################################################

# GLOBAL LEVEL
lower_ind <- (iter_num - iter_num*.95)/2
upper_ind <- iter_num - lower_ind
mean_nullmodel_I_global <- c()
lower_nullmodel_I_global <- c()
upper_nullmodel_I_global <- c()
for (t in 1:length(years)){
  mean_nullmodel_I_global[t] <- mean(nullmodel_I_global[,t])
  lower_nullmodel_I_global[t] <- sort(nullmodel_I_global[,t])[lower_ind]
  upper_nullmodel_I_global[t] <- sort(nullmodel_I_global[,t])[upper_ind]
}

# plot
pdf(file="Figures/global_penguindex_plot_nullmodel.pdf")
plot(years, nullmodel_I_global[1,], type = "l", col = alpha("gray", 0.075), ylim = c(0.95, 1.05), ylab = "Global Pygoscelis Penguindex", xlab = "Season", axes = FALSE)
axis(1)
axis(2)
for (iter in 2:iter_num){
  lines(years, nullmodel_I_global[iter,], type = "l", col = alpha("gray", 0.075))
}
lines(years, mean_nullmodel_I_global, type = "l", col = "black", lwd = 2)
lines(years, lower_nullmodel_I_global, type = "l", col = "white", lwd = 1.5)
lines(years, upper_nullmodel_I_global, type = "l", col = "white", lwd = 1.5)
dev.off()

LPI_df <- data.frame(LPI_final = mean_nullmodel_I_global, CI_low = lower_nullmodel_I_global, CI_high = upper_nullmodel_I_global, row.names = years)


###############################################################

# SPECIES LEVEL
species_list_longname <- c("Adelie", "Chinstrap", "Gentoo")
LPI_specieslevel_df <- data.frame()
colnames_LPI_specieslevel_df <- c()
for (sp_i in 1:length(species_list)){
  mean_nullmodel <- c()
  lower_nullmodel <- c()
  upper_nullmodel <- c()
  # pick central 9500 I values for each year
  for (t in 1:length(years)){
    mean_nullmodel[t] <- mean(nullmodel_I_specieslevel[sp_i,,t])
    lower_nullmodel[t] <- sort(nullmodel_I_specieslevel[sp_i,,t])[lower_ind]
    upper_nullmodel[t] <- sort(nullmodel_I_specieslevel[sp_i,,t])[upper_ind]
    LPI_specieslevel_df[t,3*(sp_i-1)+1] <- mean_nullmodel[t]
    LPI_specieslevel_df[t,3*(sp_i-1)+2] <- lower_nullmodel[t]
    LPI_specieslevel_df[t,3*(sp_i-1)+3] <- upper_nullmodel[t]
  }
  colnames_LPI_specieslevel_df <- c(colnames_LPI_specieslevel_df,paste("LPI_final",species_list[sp_i],sep="_"))
  colnames_LPI_specieslevel_df <- c(colnames_LPI_specieslevel_df,paste("CI_low",species_list[sp_i],sep="_"))
  colnames_LPI_specieslevel_df <- c(colnames_LPI_specieslevel_df,paste("CI_high",species_list[sp_i],sep="_"))
  
  # plot
  pdf(file=paste("Figures/",species_list[sp_i],"_penguindex_plot_nullmodel.pdf",sep=""))
  plot(years, nullmodel_I_specieslevel[sp_i,1,], type = "l", col = alpha("gray", 0.75), ylim = c((min(lower_nullmodel)*0.9),(max(upper_nullmodel)*1.1)), ylab = paste("Global",species_list_longname[sp_i],"Penguindex",sep=" "), xlab = "Season", axes = FALSE)
  axis(1)
  axis(2)
  for (iter in 2:iter_num){
    lines(years, nullmodel_I_specieslevel[sp_i,iter,], type = "l", col = alpha("gray", 0.75))
  }
  lines(years, mean_nullmodel, type = "l", col = "black", lwd = 2)
  lines(years, lower_nullmodel, type = "l", col = "white", lwd = 1.5)
  lines(years, upper_nullmodel, type = "l", col = "white", lwd = 1.5)
  dev.off()
}

colnames(LPI_specieslevel_df) <- colnames_LPI_specieslevel_df
rownames(LPI_specieslevel_df) <- years


###############################################################

# REGION LEVEL
LPI_speciesregionlevel_df <- list(data.frame(), data.frame(), data.frame())
for (sp_i in 1:length(species_list)){
  colnames_LPI_speciesregionlevel_df <- c()
  for (reg in species_regions[[sp_i]]){
    mean_nullmodel <- c()
    lower_nullmodel <- c()
    upper_nullmodel <- c()
    reg_i_rel <- which(species_regions[[sp_i]] == reg)
    for (t in 1:length(years)){
      mean_nullmodel[t] <- mean(nullmodel_I_speciesregionlevel[sp_i,reg,,t])
      lower_nullmodel[t] <- sort(nullmodel_I_speciesregionlevel[sp_i,reg,,t])[lower_ind]
      upper_nullmodel[t] <- sort(nullmodel_I_speciesregionlevel[sp_i,reg,,t])[upper_ind]
      LPI_speciesregionlevel_df[[sp_i]][t,3*(reg_i_rel-1)+1] <- mean_nullmodel[t]
      LPI_speciesregionlevel_df[[sp_i]][t,3*(reg_i_rel-1)+2] <- lower_nullmodel[t]
      LPI_speciesregionlevel_df[[sp_i]][t,3*(reg_i_rel-1)+3] <- upper_nullmodel[t]
    }
    colnames_LPI_speciesregionlevel_df <- c(colnames_LPI_speciesregionlevel_df,paste("LPI_final",reg,sep="_"))
    colnames_LPI_speciesregionlevel_df <- c(colnames_LPI_speciesregionlevel_df,paste("CI_low",reg,sep="_"))
    colnames_LPI_speciesregionlevel_df <- c(colnames_LPI_speciesregionlevel_df,paste("CI_high",reg,sep="_"))
    
    # plot
    pdf(file=paste("Figures/",species_list[sp_i],"_region",reg,"_penguindex_plot_nullmodel.pdf",sep=""))
    plot(years, nullmodel_I_speciesregionlevel[sp_i,reg,1,], type = "l", col = alpha("gray", 0.075), ylim = c(min(0.90,min(lower_nullmodel)*0.9),max(1.30,max(upper_nullmodel)*1.1)), ylab = paste(species_list_longname[sp_i],"Penguindex Region",reg,sep=" "), xlab = "Season", axes = FALSE)
    axis(1)
    axis(2)
    for (iter in 2:iter_num){
      lines(years, nullmodel_I_speciesregionlevel[sp_i,reg,iter,], type = "l", col = alpha("gray", 0.075))
    }
    lines(years, mean_nullmodel, type = "l", col = "black", lwd = 2)
    lines(years, lower_nullmodel, type = "l", col = "white", lwd = 1.5)
    lines(years, upper_nullmodel, type = "l", col = "white", lwd = 1.5)
    dev.off()
  }
  colnames(LPI_speciesregionlevel_df[[sp_i]]) <- colnames_LPI_speciesregionlevel_df
  rownames(LPI_speciesregionlevel_df[[sp_i]]) <- years
}



######## print & save
LPI_null_df <- LPI_df
LPI_specieslevel_null_df <- LPI_specieslevel_df
LPI_speciesregionlevel_null_df <- LPI_speciesregionlevel_df

print(paste("global", round(LPI_null_df$LPI_final[length(LPI_null_df$LPI_final)],4)))


for (sp_i in 1:length(species_list)){
  print(paste(species_list[sp_i], round(LPI_specieslevel_null_df[length(LPI_specieslevel_null_df[,1]),3*(sp_i-1)+1],4)))
  for (reg in species_regions[[sp_i]]){
    reg_i_rel <- which(species_regions[[sp_i]] == reg)
    print(paste(region_list[reg],round(tail(LPI_speciesregionlevel_null_df[[sp_i]][which(!is.na(LPI_speciesregionlevel_null_df[[sp_i]][,3*(reg_i_rel-1)+1])),3*(reg_i_rel-1)+1], n=1),4)))
  }
}


save(LPI_null_df, LPI_specieslevel_null_df, LPI_speciesregionlevel_null_df, file = "nullmodel_finisheddata.RData")


