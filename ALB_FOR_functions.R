
# DATA FUNCTIONS ----------------------------------------------------------

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# GPS functions -----------------------------------------------------------

calc_speed <- function(mydf) {
  
  n = nrow(mydf)
  
  mydf$dist_next <- c(gcd.hf(mydf$longitude[2:n],
                             mydf$latitude[2:n],
                             mydf$longitude[1:(n-1)],
                             mydf$latitude[1:(n-1)]),NA)
  
  dt = c(as.numeric(difftime(mydf$datetime[2:n], mydf$datetime[1:(n-1)], 
                             units = "secs")), NA)
  mydf$calc_speed <- mydf$dist_next*1000/dt
  
  return(mydf$calc_speed)
  
}




gcd.hf <- function(long1, lat1, long2, lat2) { 
  R <- 6371 # Earth mean radius [km]
  deg2rad <- function(deg) return(deg*pi/180)
  long1 <- deg2rad(long1)
  long2 <- deg2rad(long2)
  lat1 <- deg2rad(lat1)
  lat2 <- deg2rad(lat2)
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  
  coo2 <- function(x){2 * asin(min(1,sqrt(x)))}
  c <- lapply(a, coo2)
  c <- do.call("c", c)
  d = R * c
  return(d) # Distance in km
}




# STATS FUNCTIONS ---------------------------------------------------------

log_to_percent <- function(coef) {
  odds_ratio <- exp(coef)
  return((odds_ratio - 1) * 100)
}




process_interaction_estimates <- function(var_name, fixef_data, posterior_samples.df, colony_interaction) {
  
  # Get interaction SE
  main_effect_samples <- posterior_samples.df[[var_name]]
  interaction_samples <- posterior_samples.df[[colony_interaction]]
  
  ## Calculate covariance between the two estimates
  covariance <- cov(main_effect_samples, interaction_samples)
  
  ## Standard errors from fixef table
  se_main <- fixef_data$Est.Error[which(rownames(fixef_data) == gsub("b_", "", var_name))] 
  se_interaction <- fixef_data$Est.Error[which(rownames(fixef_data) == gsub("b_", "", colony_interaction))]
  
  ## Combine the SEs using the formula
  se_total <- sqrt(se_main^2 + se_interaction^2 + 2 * covariance)
  
  # Calculate estimates and confidence intervals
  est <- fixef_data$Estimate[which(rownames(fixef_data) == gsub("b_", "", var_name))] +
    fixef_data$Estimate[which(rownames(fixef_data) == gsub("b_", "", colony_interaction))]
  error <- se_total
  LCI <- fixef_data$`l-95% CI`[which(rownames(fixef_data) == gsub("b_", "", var_name))] +
    fixef_data$`l-95% CI`[which(rownames(fixef_data) == gsub("b_", "", colony_interaction))]
  UCI <- fixef_data$`u-95% CI`[which(rownames(fixef_data) == gsub("b_", "", var_name))] +
    fixef_data$`u-95% CI`[which(rownames(fixef_data) == gsub("b_", "", colony_interaction))]
  
  # Compute probabilities
  samples <- posterior_samples.df[[var_name]] + posterior_samples.df[[colony_interaction]]
  prop_rope <- round(prop_rope.func(samples), digits = 2)
   
  # Create a row to append
  result <- cbind(est, error, LCI, UCI, matrix(NA, 1, 3), exp(est), exp(LCI), exp(UCI), prop_rope)
  
  rownames(result) <- paste("col", var_name, sep = ".")
  
  return(result)
}



prop_rope.func <- function (variable) {  
  mean(variable > rope_lower & variable < rope_upper) }

