
# DATA FUNCTIONS ----------------------------------------------------------

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


save_dharma_diag <- function(model, response_var, covariate_vec, filename, title) {
  # save residual outputs
  
  preds <- posterior_predict(model)
  
  sim_res <- createDHARMa(
    simulatedResponse = t(preds), 
    observedResponse = model$data[[response_var]], 
    fittedPredictedResponse = apply(preds, 2, median),
    integerResponse = TRUE)
  
  png(filename, width = 1500, height = 500, res = 150)
  par(mfrow = c(1, 3), oma = c(0, 0, 2, 0), mar = c(5, 5, 4, 2)) 
  
  plotQQunif(sim_res)
  plotResiduals(sim_res, main = "Residuals vs. Predicted")
  plotResiduals(sim_res, covariate_vec, 
                quantreg = TRUE, 
                xlab = "Trip Duration")
  mtext(title, outer = TRUE, cex = 1.2, font = 2)
  
  dev.off() 
}


# STATS FUNCTIONS ---------------------------------------------------------

# model = bbal_incub.brms
# behaviour_data = bbal_incub.behaviour
# behaviour_var = "median_trip.days"
# values = c(2,5)
# colony = "birdisland"

calculate_prob_change <- function(model, 
                                  behaviour_data, 
                                  behaviour_var, 
                                  values, 
                                  colony) {
  
  # Create new dataframe with specified values
  base_row <- data.frame(
    colony = colony,
    median_trip.days = mean(behaviour_data$median_trip.days, na.rm = TRUE),
    diff_trip_abs.days = mean(behaviour_data$diff_trip_abs.days, na.rm = TRUE),
    diff_var_abs.days = mean(behaviour_data$diff_var_abs.days, na.rm = TRUE),
    season = NA  )
  
  new_df <- base_row[rep(1, length(values)), ]
  new_df[[behaviour_var]] <- values
  
  # Predict posterior probabilities
  pp <- posterior_epred(model, newdata = new_df, allow_new_levels = TRUE)
  pp_summary <- apply(pp, 2, function(x) c(mean = mean(x), quantile(x, probs = c(0.025, 0.975))))
  
  # Calculate differences
  prob_change <- (round(pp_summary[1, 2] - pp_summary[1, 1], 2)) * 100
  ci_change_low <- (pp_summary[2, 2] - pp_summary[2, 1]) * 100
  ci_change_high <- (pp_summary[3, 2] - pp_summary[3, 1]) * 100
  ci_change_vals <- sort(c(ci_change_low, ci_change_high))
  ci_change <- paste0("[", round(ci_change_vals[1], 2), "-", round(ci_change_vals[2], 2), "]")
  
  return(list(prob_change = prob_change, ci_change = ci_change))
}






# brms_fit = bbal_incub.brms_weak
# main_terms = predictors
# int_terms = interaction_terms
# colony_names = c("Bird Island", "Kerguelen")
# 
# extract_posteriors(brms_fit, main_effects, int_terms, colony_names)

extract_posteriors <- function(brms_fit, main_terms, int_terms, colony_names) {
  
  draws <- as_draws_df(brms_fit)
  out <- vector("list", length(main_terms))
  
  for (i in seq_along(main_terms)) {
    
    ref_est <- draws[[main_terms[i]]]
    int_est <- draws[[int_terms[i]]]
    
    out[[i]] <- data.frame(
      variable = main_terms[i],
      value = c(ref_est, ref_est + int_est),
      colony = rep(colony_names, each = length(ref_est)) )
  }
  
  bind_rows(out)
}




log_to_percent <- function(coef) {
  odds_ratio <- exp(coef)
  return((odds_ratio - 1) * 100)
}




# var_name = "b_debt.days"
# fixef_data = waal_brooding.fixef
# posterior_samples.df = posterior_samples.brooding_waal
# colony_interaction = "b_debt.days:colonycrozet"

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
  LCI <- NA 
  UCI <- NA
  
  # Compute probabilities
  samples <- posterior_samples.df[[var_name]] + posterior_samples.df[[colony_interaction]]
  prop_rope <-as.numeric(rope(samples, range = c(rope_lower, rope_upper), ci = 0.95))
  prop_rope <- round(prop_rope, digits = 2)
   
  # Create a row to append
  result <- cbind(est, error, LCI, UCI, matrix(NA, 1, 3), exp(est), exp(LCI), exp(UCI), prop_rope)
  
  rownames(result) <- paste("col", var_name, sep = ".")
  
  return(result)
}

