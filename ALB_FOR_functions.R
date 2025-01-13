
# LOG_TO_PERCENT ----------------------------------------------------------
# Convert logits to percentages
log_to_percent <- function(coef) {
  odds_ratio <- exp(coef)
  return((odds_ratio - 1) * 100)
}

# PROCESS_INTERACTION_ESTIMATES: CALCULATE ESTIMATES FOR EACH INTE --------

process_interaction_estimates <- function(var_name, fixef_data, posterior_samples.df, colony_interaction) {
  # Calculate estimates and confidence intervals
  est <- fixef_data$Estimate[which(rownames(fixef_data) == gsub("b_", "", var_name))] +
    fixef_data$Estimate[which(rownames(fixef_data) == gsub("b_", "", colony_interaction))]
  LCI <- fixef_data$`l-95% CI`[which(rownames(fixef_data) == gsub("b_", "", var_name))] +
    fixef_data$`l-95% CI`[which(rownames(fixef_data) == gsub("b_", "", colony_interaction))]
  UCI <- fixef_data$`u-95% CI`[which(rownames(fixef_data) == gsub("b_", "", var_name))] +
    fixef_data$`u-95% CI`[which(rownames(fixef_data) == gsub("b_", "", colony_interaction))]
  
  # Compute probabilities
  samples <- posterior_samples.df[[var_name]] + posterior_samples.df[[colony_interaction]]
  above_zero <- round(above_zero.func(samples), digits = 2)
  below_zero <- round(below_zero.func(samples), digits = 2)
  
  # Create a row to append
  result <- cbind(matrix(NA, 1, 12), est, LCI, UCI, below_zero, above_zero)
  rownames(result) <- paste("col", var_name, sep = ".")
  
  return(result)
}


# ABOVE AND BELOW ZERO: CALCULATE PROPORTION OF SAMPLES -------------------

above_zero.func <- function (variable) { mean(variable > 0) }
below_zero.func <- function (variable) { mean(variable < 0) }

