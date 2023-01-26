library(dplyr)
library(ggplot2)

# Global - species - region
# GLOBAL penguin LPI / Penguindex -- for each species, calc lambdas for each region (weight by # of pops), aggregate over spp

# input: MAPPPD output full **1980**-2019 TIME SERIES WITH AT LEAST 2 OBSERVATIONS
# 0 pop values (which are, according to MAPPPD model, observed no occupancy at end/beginning of time series)
# are changed to 1% of the mean pop count for that time series
# output: all Penguindex indices

start_time <- Sys.time()

lpi_data_all <- read.csv("../Filtering time series/prepared_all_species_atleast2_1980.csv", na.strings = "NULL")

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


############
#
#
# confidence intervals
#
#
# sampling from the posteriors
#
#
#
##############
#
# load in zstates with posterior mean and sd
#
adelie_z_df <- read.csv("../BSSM outputs/ADPE_abundance_lambda.csv", na.strings = "NULL")
chinstrap_z_df <- read.csv("../BSSM outputs/CHPE_abundance_lambda.csv", na.strings = "NULL")
gentoo_z_df <- read.csv("../BSSM outputs/GEPE_abundance_lambda.csv", na.strings = "NULL")
all_species_z_df <- rbind(adelie_z_df, chinstrap_z_df, gentoo_z_df)
#
#
iter_num <- 1000 ######
samppost_d <- array(NA, dim = c(length(species_list), length(region_list), max(unlist(species_regions_site_list_lengths)), length(years)))
samppost_I_speciesregionlevel <- array(NA, dim = c(length(species_list), length(region_list), iter_num, length(years)))
samppost_I_specieslevel <- array(NA, dim = c(length(species_list), iter_num, length(years)))
samppost_I_global <- matrix(NA, nrow = iter_num, ncol = length(years))
for (iter in 1:iter_num){
  #
  # sample from posteriors
  #
  lpi_data_all_prime <- lpi_data_all
  samppost_d_bar_region <- array(NA, dim = c(length(species_list), length(region_list), length(years)))
  samppost_d_bar_species <- matrix(NA, nrow = length(species_list), ncol = length(years))
  samppost_I_global[iter,1] <- 1
  #################
  for (sp_i in 1:length(species_list)){
    species <- species_list[sp_i]
    lpi_data_species <- lpi_data_all[which(lpi_data_all$Species_id == species),]
    region_weights <- c()
    for (reg in species_regions[[sp_i]]){
      reg_i_rel <- which(species_regions[[sp_i]] == reg)
      lpi_data <- lpi_data_species[which(lpi_data_species$Region_id == reg),]
      site_list <- lpi_data$Site_ID # number of sites in region for that species
      for (i in 1:length(site_list)){
        site_id <- lpi_data[i, "Site_ID"]
        N <- as.numeric(lpi_data[i,] %>% dplyr::select(starts_with("X")))
        N <- N[length(N)-length(years)+1:length(years)]
        firstcol <- which(colnames(lpi_data) == "X1980") # first col
        N_prime <- rep(NA, times = length(N))
        for (t in 1:length(years)){
          ###
          mean_it <- as.numeric(all_species_z_df$lza_mean[(all_species_z_df$species_id == species) & (all_species_z_df$site_id == site_id) & (all_species_z_df$season == years[t])])
          sd_it <- as.numeric(all_species_z_df$lza_sd[(all_species_z_df$species_id == species) & (all_species_z_df$site_id == site_id) & (all_species_z_df$season == years[t])])
          if (!is.na(sd_it)){
            N_it_prime <- rlnorm(1, meanlog = mean_it, sdlog = sd_it)
          } else{
            N_it_prime <- 0.01*mean(N[which(N != 0)]) # onepercent
          }
          N_prime[t] <- N_it_prime
          if (t >= 2){
            samppost_d[sp_i,reg,i,t] <- log10(N_prime[t]/N_prime[t-1])
          }
        }
        lpi_data_all_prime[(lpi_data_all_prime$Species_id == species) & (lpi_data_all_prime$Site_ID == site_id),firstcol:ncol(lpi_data_all_prime)] <- N_prime
      }
      # now for each region, calculate d_bar(spp,reg,t) & region-level index
      samppost_I_speciesregionlevel[sp_i,reg,iter,1] <- 1
      for (t in 2:length(years)){
        samppost_d_bar_region[sp_i,reg,t] <- (1/species_regions_site_list_lengths[[sp_i]][reg_i_rel])*sum(samppost_d[sp_i,reg,,t], na.rm = TRUE) # over each site
        samppost_I_speciesregionlevel[sp_i,reg,iter,t] <- samppost_I_speciesregionlevel[sp_i,reg,iter,t-1]*10^(samppost_d_bar_region[sp_i,reg,t])
      }
      # region weights
      region_weights[reg] <- species_regions_site_list_lengths[[sp_i]][reg_i_rel]/species_total_num_pops[sp_i]
    }
    region_weights <- c(region_weights, rep(NA, times=length(region_list)-length(region_weights)))
    # now for each spp, calculate d_bar(spp,t) & species-level index
    samppost_I_specieslevel[sp_i,iter,1] <- 1
    for (t in 2:length(years)){
      samppost_d_bar_species[sp_i,t] <- (1/length(species_regions[[sp_i]]))*sum(samppost_d_bar_region[sp_i,,t]*region_weights, na.rm = TRUE)
      samppost_I_specieslevel[sp_i,iter,t] <- samppost_I_specieslevel[sp_i,iter,t-1]*10^(samppost_d_bar_species[sp_i,t])
    }
  }
  # calculate global index
  for (t in 2:length(years)){
    samppost_I_global[iter,t] <- samppost_I_global[iter,t-1]*10^(mean(samppost_d_bar_species[,t]))
  }
}

end_time <- Sys.time()

end_time-start_time


###############################################################
# SAVE
save(lpi_data_all, species_regions_site_list_lengths, samppost_d, samppost_d_bar_region, samppost_d_bar_species, samppost_I_speciesregionlevel, samppost_I_specieslevel, samppost_I_global, file = "penguindex_data.RData")


###############################################################
# LOAD
load("penguindex_data.RData")

# include all iterations?
iter_num <- 1000 
samppost_I_speciesregionlevel <- samppost_I_speciesregionlevel[,,1:iter_num,]
samppost_I_specieslevel <- samppost_I_specieslevel[,1:iter_num,]
samppost_I_global <- samppost_I_global[1:iter_num,]

###############################################################

# GLOBAL LEVEL
lower_ind <- (iter_num - iter_num*.95)/2
upper_ind <- iter_num - lower_ind
mean_samppost_I_global <- c()
lower_samppost_I_global <- c()
upper_samppost_I_global <- c()
for (t in 1:length(years)){
  mean_samppost_I_global[t] <- mean(samppost_I_global[,t])
  lower_samppost_I_global[t] <- sort(samppost_I_global[,t])[lower_ind]
  upper_samppost_I_global[t] <- sort(samppost_I_global[,t])[upper_ind]
}

# plot
pdf(file="Figures/global_penguindex_plot.pdf")
plot(years, samppost_I_global[1,], type = "l", col = alpha("gray", 0.075), ylim = c(0.90, 1.30), ylab = "Global Pygoscelis Penguindex", xlab = "Season", axes = FALSE)
axis(1)
axis(2)
for (iter in 2:iter_num){
  lines(years, samppost_I_global[iter,], type = "l", col = alpha("gray", 0.075))
}
lines(years, mean_samppost_I_global, type = "l", col = "black", lwd = 2)
lines(years, lower_samppost_I_global, type = "l", col = "white", lwd = 1.5)
lines(years, upper_samppost_I_global, type = "l", col = "white", lwd = 1.5)
dev.off()

LPI_df <- data.frame(LPI_final = mean_samppost_I_global, CI_low = lower_samppost_I_global, CI_high = upper_samppost_I_global, row.names = years)


###############################################################

# SPECIES LEVEL
species_list_longname <- c("Adelie", "Chinstrap", "Gentoo")
LPI_specieslevel_df <- data.frame()
colnames_LPI_specieslevel_df <- c()
for (sp_i in 1:length(species_list)){
  mean_samppost <- c()
  lower_samppost <- c()
  upper_samppost <- c()
  # pick central 9500 I values for each year
  for (t in 1:length(years)){
    mean_samppost[t] <- mean(samppost_I_specieslevel[sp_i,,t])
    lower_samppost[t] <- sort(samppost_I_specieslevel[sp_i,,t])[lower_ind]
    upper_samppost[t] <- sort(samppost_I_specieslevel[sp_i,,t])[upper_ind]
    LPI_specieslevel_df[t,3*(sp_i-1)+1] <- mean_samppost[t]
    LPI_specieslevel_df[t,3*(sp_i-1)+2] <- lower_samppost[t]
    LPI_specieslevel_df[t,3*(sp_i-1)+3] <- upper_samppost[t]
  }
  colnames_LPI_specieslevel_df <- c(colnames_LPI_specieslevel_df,paste("LPI_final",species_list[sp_i],sep="_"))
  colnames_LPI_specieslevel_df <- c(colnames_LPI_specieslevel_df,paste("CI_low",species_list[sp_i],sep="_"))
  colnames_LPI_specieslevel_df <- c(colnames_LPI_specieslevel_df,paste("CI_high",species_list[sp_i],sep="_"))
  
  # plot
  pdf(file=paste("Figures/",species_list[sp_i],"_penguindex_plot.pdf",sep=""))
  plot(years, samppost_I_specieslevel[sp_i,1,], type = "l", col = alpha("gray", 0.075), ylim = c(min(0.90,min(lower_samppost)*0.9),max(1.30,max(upper_samppost)*1.1)), ylab = paste("Global",species_list_longname[sp_i],"Penguindex",sep=" "), xlab = "Season", axes = FALSE)
  axis(1)
  axis(2)
  for (iter in 2:iter_num){
    lines(years, samppost_I_specieslevel[sp_i,iter,], type = "l", col = alpha("gray", 0.075))
  }
  lines(years, mean_samppost, type = "l", col = "black", lwd = 2)
  lines(years, lower_samppost, type = "l", col = "white", lwd = 1.5)
  lines(years, upper_samppost, type = "l", col = "white", lwd = 1.5)
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
    mean_samppost <- c()
    lower_samppost <- c()
    upper_samppost <- c()
    reg_i_rel <- which(species_regions[[sp_i]] == reg)
    for (t in 1:length(years)){
      mean_samppost[t] <- mean(samppost_I_speciesregionlevel[sp_i,reg,,t])
      lower_samppost[t] <- sort(samppost_I_speciesregionlevel[sp_i,reg,,t])[lower_ind]
      upper_samppost[t] <- sort(samppost_I_speciesregionlevel[sp_i,reg,,t])[upper_ind]
      LPI_speciesregionlevel_df[[sp_i]][t,3*(reg_i_rel-1)+1] <- mean_samppost[t]
      LPI_speciesregionlevel_df[[sp_i]][t,3*(reg_i_rel-1)+2] <- lower_samppost[t]
      LPI_speciesregionlevel_df[[sp_i]][t,3*(reg_i_rel-1)+3] <- upper_samppost[t]
    }
    colnames_LPI_speciesregionlevel_df <- c(colnames_LPI_speciesregionlevel_df,paste("LPI_final",reg,sep="_"))
    colnames_LPI_speciesregionlevel_df <- c(colnames_LPI_speciesregionlevel_df,paste("CI_low",reg,sep="_"))
    colnames_LPI_speciesregionlevel_df <- c(colnames_LPI_speciesregionlevel_df,paste("CI_high",reg,sep="_"))
    
    # plot
    pdf(file=paste("Figures/",species_list[sp_i],"_region",reg,"_penguindex_plot.pdf",sep=""))
    plot(years, samppost_I_speciesregionlevel[sp_i,reg,1,], type = "l", col = alpha("gray", 0.075), ylim = c(min(0.90,min(lower_samppost)*0.9),max(1.30,max(upper_samppost)*1.1)), ylab = paste(species_list_longname[sp_i],"Penguindex Region",reg,sep=" "), xlab = "Season", axes = FALSE)
    axis(1)
    axis(2)
    for (iter in 2:iter_num){
      lines(years, samppost_I_speciesregionlevel[sp_i,reg,iter,], type = "l", col = alpha("gray", 0.075))
    }
    lines(years, mean_samppost, type = "l", col = "black", lwd = 2)
    lines(years, lower_samppost, type = "l", col = "white", lwd = 1.5)
    lines(years, upper_samppost, type = "l", col = "white", lwd = 1.5)
    dev.off()
  }
  colnames(LPI_speciesregionlevel_df[[sp_i]]) <- colnames_LPI_speciesregionlevel_df
  rownames(LPI_speciesregionlevel_df[[sp_i]]) <- years
}


###############################################################

######## print & save
print(paste("global", round(LPI_df$LPI_final[length(LPI_df$LPI_final)],4)))

for (sp_i in 1:length(species_list)){
  print(paste(species_list[sp_i], round(LPI_specieslevel_df[length(LPI_specieslevel_df[,1]),3*(sp_i-1)+1],4)))
  for (reg in species_regions[[sp_i]]){
    reg_i_rel <- which(species_regions[[sp_i]] == reg)
    print(paste(region_list[reg],round(tail(LPI_speciesregionlevel_df[[sp_i]][which(!is.na(LPI_speciesregionlevel_df[[sp_i]][,3*(reg_i_rel-1)+1])),3*(reg_i_rel-1)+1], n=1),4)))
  }
}

save(LPI_df, LPI_specieslevel_df, LPI_speciesregionlevel_df, file = "penguindex_finisheddata.RData")


write.csv(LPI_df, file = "global_finisheddata.csv")
write.csv(LPI_specieslevel_df, file = "species_level_finisheddata.csv")
write.csv(LPI_speciesregionlevel_df, file = "region_level_finisheddata.csv")



###############################################################
###############################################################
#
#
# break point/change point analysis & final plots
#
#
###############################################################
###############################################################

library(segmented)
library(stats)

#######################
#
# load null model data
#
########################

load("nullmodel_finisheddata.RData")


########################################
# GLOBAL LEVEL
breakpoints_global <- list()
global <- data.frame(year = years, index = LPI_df$LPI_final)
model <- glm(index~year, data=global)
os <- selgmented(model, Kmax=10, type="bic")
breakpoints <- os$psi
breakpoints_global <- breakpoints[,2]

# plot
pdf(file="Figures/global_penguindex_plot_full.pdf")
plot(years, samppost_I_global[1,], type = "l", col = alpha("gray", 0.075), ylim = c(0.90, 1.30), ylab = "Global Pygoscelis Penguindex", xlab = "Season", axes = FALSE)
axis(1)
axis(2)
for (brkpt in 1:length(breakpoints[,2])){
  abline(v = as.numeric(breakpoints[,2][brkpt]), lwd = 2, col = alpha("palegreen4", 0.55), lty = "longdash")
}
lines(years, LPI_null_df$LPI_final, type = "l", lty = "solid", col = alpha("royalblue1",0.99), lwd = 2.25)
for (iter in 2:iter_num){
  lines(years, samppost_I_global[iter,], type = "l", col = alpha("gray", 0.075))
}
lines(years, mean_samppost_I_global, type = "l", col = "black", lwd = 2)
lines(years, lower_samppost_I_global, type = "l", col = "white", lwd = 1.5)
lines(years, upper_samppost_I_global, type = "l", col = "white", lwd = 1.5)
dev.off()

########################################
# SPECIES LEVEL
breakpoints_species <- list(list(), list(), list())
for (sp_i in 1:length(species_list)){
  spp <- data.frame(year = years, index = LPI_specieslevel_df[,(3*(sp_i-1)+1)])
  model <- glm(index~year, data=spp)
  os <- selgmented(model, Kmax=10, type="bic")
  breakpoints <- os$psi
  breakpoints_species[[sp_i]] <- breakpoints[,2]
  
  # plot
  pdf(file=paste("Figures/",species_list[sp_i],"_penguindex_plot_full.pdf",sep=""))
  plot(years, samppost_I_specieslevel[sp_i,1,], type = "l", col = alpha("gray", 0.075), ylim = c(min(0.90,min(LPI_specieslevel_df[,3*(sp_i-1)+2])*0.9),max(1.30,max(LPI_specieslevel_df[t,3*(sp_i-1)+3])*1.1)), ylab = paste("Global",species_list_longname[sp_i],"Penguindex",sep=" "), xlab = "Season", axes = FALSE)
  axis(1)
  axis(2)
  for (brkpt in 1:length(breakpoints[,2])){
    abline(v = as.numeric(breakpoints[,2][brkpt]), lwd = 2, col = alpha("palegreen4", 0.55), lty = "longdash")
  }
  lines(years, LPI_specieslevel_null_df[,3*(sp_i-1)+1], type = "l", lty = "solid", col = alpha("royalblue1",0.99), lwd = 2.25)
  for (iter in 2:iter_num){
    lines(years, samppost_I_specieslevel[sp_i,iter,], type = "l", col = alpha("gray", 0.075))
  }
  lines(years, LPI_specieslevel_df[,3*(sp_i-1)+1], type = "l", col = "black", lwd = 2)
  lines(years, LPI_specieslevel_df[,3*(sp_i-1)+2], type = "l", col = "white", lwd = 1.5)
  lines(years, LPI_specieslevel_df[,3*(sp_i-1)+3], type = "l", col = "white", lwd = 1.5)
  dev.off()
}

########################################
# REGION LEVEL
empty <- list(list(), list(), list(), list(), list(), list(), list(), list())
breakpoints_species_region <- list(empty, empty, empty)
for (sp_i in 1:length(species_list)){
  for (reg in species_regions[[sp_i]]){
    reg_i_rel <- which(species_regions[[sp_i]] == reg)
    regg <- data.frame(year = years, index = LPI_speciesregionlevel_df[[sp_i]][,3*(reg_i_rel-1)+1]) 
    model <- glm(index~year, data=regg)
    os <- selgmented(model, Kmax=10, type="bic")
    breakpoints <- os$psi
    breakpoints_species_region[[sp_i]][[reg]] <- breakpoints[,2]
    
    # plot
    pdf(file=paste("Figures/",species_list[sp_i],"_region",reg,"_penguindex_plot_full.pdf",sep=""))
    plot(years, samppost_I_speciesregionlevel[sp_i,reg,1,], type = "l", col = alpha("gray", 0.075), ylim = c(min(0.90,min(LPI_speciesregionlevel_df[[sp_i]][,3*(reg_i_rel-1)+2])*0.9),max(1.30,max(LPI_speciesregionlevel_df[[sp_i]][,3*(reg_i_rel-1)+3])*1.1)), ylab = paste(species_list_longname[sp_i],"Penguindex Region",reg,sep=" "), xlab = "Season", axes = FALSE)
    axis(1)
    axis(2)
    for (brkpt in 1:length(breakpoints[,2])){
      abline(v = as.numeric(breakpoints[,2][brkpt]), lwd = 2, col = alpha("palegreen4", 0.55), lty = "longdash")
    }
    lines(years, LPI_speciesregionlevel_null_df[[sp_i]][,3*(reg_i_rel-1)+1], type = "l", lty = "solid", col = alpha("royalblue1",0.99), lwd = 2.25)
    for (iter in 2:iter_num){
      lines(years, samppost_I_speciesregionlevel[sp_i,reg,iter,], type = "l", col = alpha("gray", 0.075))
    }
    lines(years, LPI_speciesregionlevel_df[[sp_i]][,3*(reg_i_rel-1)+1], type = "l", col = "black", lwd = 2)
    lines(years, LPI_speciesregionlevel_df[[sp_i]][,3*(reg_i_rel-1)+2], type = "l", col = "white", lwd = 1.5)
    lines(years, LPI_speciesregionlevel_df[[sp_i]][,3*(reg_i_rel-1)+3], type = "l", col = "white", lwd = 1.5)
    dev.off()
  }
}



save(breakpoints_global, breakpoints_species, breakpoints_species_region, file = "penguindex_breakpointdata.RData")

