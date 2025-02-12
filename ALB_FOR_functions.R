
# LOG_TO_PERCENT ----------------------------------------------------------
# Convert logits to percentages
log_to_percent <- function(coef) {
  odds_ratio <- exp(coef)
  return((odds_ratio - 1) * 100)
}

# PROCESS_INTERACTION_ESTIMATES: CALCULATE ESTIMATES FOR EACH INTE --------
# var_name = "b_debt.days"
# fixef_data = bbal_incub.fixef
# posterior_samples.df = posterior_samples.incub_bba
# colony_interaction = "b_debt.days:colonykerguelen"

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


# CALCUATE ROPE AND MARGINAL ROPE -------------------

prop_rope.func <- function (variable) {  
  mean(variable > rope_lower & variable < rope_upper) }

