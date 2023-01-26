library(dplyr)

####### download Living Planet Database here: https://www.livingplanetindex.org/data_portal
# ----> "LPD2022_public.csv"

#
lpd <- read.csv("LPD2022_public.csv")
pop_list <- lpd$ID
spp_list <- unique(lpd$Binomial)
avg_percent_change <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(avg_percent_change) <- c("Binomial", "Common_name", "Num_pops", "Start_year", "End_year", "Avg_p_chg")
for (i in 1:length(spp_list)){
  percent_change <- c()
  binomial <- spp_list[i]
  lpd_spp <- lpd[which(lpd$Binomial == binomial),]
  common_name <- unique(lpd_spp$Common_name)
  if (length(common_name) != 1) common_name <- common_name[which.max(nchar(common_name))] 
  num_pop <- length(lpd_spp$ID)
  start_year <- c()
  end_year <- c()
  
  for (j in 1:num_pop){
    N <- as.numeric(lpd_spp[j,] %>% dplyr::select(starts_with("X")))
    obs_inds <- which(!is.na(N))
    if (length(obs_inds) == 0) next
    start_year[j] <- which.min(N) + 1950 - 1
    end_year[j] <- max(obs_inds) + 1950 - 1
    if (start_year[j] == end_year[j]) next
    if (N[max(obs_inds)] == 0) next
    N <- N[obs_inds]
    if (length(N) == 0) next
    percent_change[j] <- (N[length(N)] - min(N))/min(N)*100
    if (min(N) == 0){
      percent_change[j] <- (N[length(N)] - min(N))*100
    }
  }
  
  if (length(start_year) == 0) next
  if (length(end_year) == 0) next
  if (length(percent_change) == 0) next
  
  num_pos_pop <- length(which(!is.na(percent_change)))
  
  avg_p_chg_spp <- round(mean(percent_change, na.rm = TRUE))
  avg_percent_change[nrow(avg_percent_change) + 1,] <- c(binomial, common_name, num_pos_pop, min(start_year), max(end_year), avg_p_chg_spp)
}

save(avg_percent_change, file = "lpd_percent_change.RData")
write.csv(avg_percent_change, file = "lpd_percent_change.csv")

write.csv(as.data.frame(quantile(as.numeric(avg_percent_change$Avg_p_chg), probs = seq(0,1,0.05))), file = "lpd_percent_change_quantile.csv")



####################
# determine quantile for average gentoo penguin growth
load("total_avg_percent_change_gent.RData")
ecdf(as.numeric(avg_percent_change$Avg_p_chg))(total_avg_percent_change_gent)

