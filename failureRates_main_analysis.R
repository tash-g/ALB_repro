## ---------------------------
##
## Script name: failureRates_main_analysis.R
##
## Purpose of script: Calculate metrics of behaviour from trip data e.g. mismatch, anomalies
## and analyse how impact breeding success
##
## Author: Dr. Natasha Gillies
##
## Created: 2024-02-21
##
## Email: gilliesne@gmail.com
##
## ---------------------------


# Load functions, packages, & data ---------------------------------------------

# Functions for GPS processing and plotting
source("ALB_FOR_functions.R")

# Define the packages
packages <- c("dplyr", "magrittr", "ggplot2", "brms", "tidyverse", "tidybayes",
              "ggpubr", "extrafont", "cowplot")

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())
# 
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages], dependencies = TRUE)
# }

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Suppress dplyr summarise warning
options(dplyr.summarise.inform = FALSE)

# Make sure using dplyr select
select <- dplyr::select


# Define common parameters ------------------------------------------------

### Set colour parameters ------------------------------------------------------

ker_col <-  "#156082" 
ker_col.high <- "#155470"
ker_col.low <- "#96d0eb"

bi_col <- "#E97132" 
bi_col.low <- "#f59b6c"
bi_col.high <- "#E97132" 


### Set labels ------------------------------------

# my_labels <- c(
#   "Within-pair trip<br>variability difference, <b><i>V</i><sub>i</sub></b>",
#   "Within-pair trip duration<br>difference, <b><i>D</i><sub>i</sub></b>",
#   "Pair-level median<br>trip duration, <b><i>T</i><sub>i</sub></b>",
#   "Fasting debt, <b><i>F</i><sub>i</sub></b>"
# )

my_labels <- c(
  expression(italic("V"[i])),
  expression(italic("D"[i])),
  expression(italic("T"[i])),
  expression(italic("F"[i]))
)


# Set up variable labeller
variable_labels <- c(
  "debt.days:colony" =  expression("Fasting debt,"~italic("F")),
  "median_trip.days:colony" = expression("Pair-level median trip duration,"~italic("T")),
  "diff_trip.days:colony" = expression("Within-pair trip duration difference,"~italic("D")),
  "diff_var.days:colony" = expression("Within-pair trip variability difference,"~italic("V"))
)

variable_labels2 <- c(
  "debt.days" =  expression("Fasting debt,"~italic("F")),
  "median_trip.days" = expression("Pair-level median trip duration,"~italic("T")),
  "diff_trip.days" = expression("Within-pair trip duration difference,"~italic("D")),
  "diff_var.days" = expression("Within-pair trip variability difference,"~italic("V"))
)

### Define ROPE bounds ------------------------------------------------------

rope_lower = -0.1
rope_upper = 0.1


### Set species ----------------------------------------------------------------

species <- c("baBI", "baKer", "waBI", "waCro")


# ______________________________ ####
# FAILURE RATES SUMMARY --------------------------------------------------------

load("Data_inputs/BAS_demo_complete.RData")
bas_fail <- bas_rs.full %>% filter(rs != "UNKNOWN" & rs != "NON-BREEDER"& !is.na(species)) %>%
  select(species, colony, season, rs, pairID) %>% distinct()

load("Data_inputs/Chize_demo_complete.RData")
chize_fail <- chize_rs.full %>% filter(rs != "UNKNOWN" & rs != "NON-BREEDER" & !is.na(species)) %>%
  select(species, colony, season, rs, pairID) %>% distinct()

failure <- rbind(bas_fail, chize_fail)

failure_summary <- failure %>%
  group_by(species, colony) %>%
  summarise(
    sample = n(),
    failed_chick = sum(rs == "FAILED_CHICK")/sample*100,
    failed_egg = sum(rs == "FAILED_EGG")/sample*100,
    failed_unknown = sum(rs == "FAILED")/sample*100,
    successsful = sum(rs == "SUCCESSFUL")/sample*100
  )

t(failure_summary)



# ______________________________ ####
# ~ * INCUBATION * ~ ###########################################################


# ANALYSIS =====================================================================
    
load("Data_inputs/all_pair_behaviour_incub.RData")

# Datasets
behaviour.bba_incub <- pair_behaviour.incub %>% filter(species == "BBA") 
behaviour.waal_incub <- pair_behaviour.incub %>% filter(species == "WAAL")

#### Examine correlations ---------------------------------------------------------------

# BBA
GGally::ggpairs(
  data = behaviour.bba_incub %>%
    select(colony, debt.days, median_trip.days, diff_trip.days, diff_var.days),
  columns = 2:5,
  mapping = ggplot2::aes(colour = colony),
  title = "Behavioural correlations by colony") +
  scale_colour_manual(values = c(bi_col, ker_col)) +
  scale_fill_manual(values = c(bi_col, ker_col))

# WAAL
GGally::ggpairs(
  data = behaviour.waal_incub %>%
    select(colony, debt.days, median_trip.days, diff_trip.days, diff_var.days),
  columns = 2:5,
  mapping = ggplot2::aes(colour = colony),
  title = "Behavioural correlations by colony") +
  scale_colour_manual(values = c(bi_col, ker_col)) +
  scale_fill_manual(values = c(bi_col, ker_col))

#### Specify priors ---------------------------------------------------------------

# BBA 
priors.bba <- c(set_prior("normal(-0.5, 1)", class = "b", coef = "debt.days"),
                set_prior("normal(-0.25, 0.5)", class = "b", coef = "median_trip.days"),       
                set_prior("normal(-0.5, 0.5)", class = "b", coef = "diff_trip.days"),
                set_prior("normal(-0.12, 0.5)", class = "b", coef = "diff_var.days"))

# WAAL
priors.waal <- c( set_prior("normal(-0.5, 1)", class = "b", coef = "debt.days"),
                  set_prior("normal(-0.25, 0.5)", class = "b", coef = "median_trip.days"),       
                  set_prior("normal(-0.5, 0.5)", class = "b", coef = "diff_trip.days"),
                  set_prior("normal(-0.12, 0.5)", class = "b", coef = "diff_var.days"))

#### Specify models ----------------------------------------------------------------

# Model structure
bf.incub <- brms::bf(breeding_outcome.bin ~
                       debt.days * colony +
                       median_trip.days * colony +
                       diff_trip.days * colony +
                       diff_var.days * colony +
                       (1|season) + (1|pairID), 
                       family = "bernoulli")

## Fit models
incub_brms.bba <- brm(bf.incub,
              data = behaviour.bba_incub, 
              cores = 4, chains = 4, 
              iter = 20000, warmup = 10000, thin = 10,
              control = list(adapt_delta = 0.9999, max_treedepth = 14),
              prior = priors.bba)

summary(incub_brms.bba)
pp_check(incub_brms.bba, ndraw = 50)
save(incub_brms.bba, file = "Data_outputs/bba_incub_brms_model.RData")

incub_brms.waal <- brm(bf.incub,
                      data = behaviour.waal_incub, 
                      cores = 4, chains = 4, 
                      iter = 20000, warmup = 10000, thin = 10,
                      control = list(adapt_delta = 0.9999, max_treedepth = 14),
                      prior = priors.waal)

summary(incub_brms.waal)
pp_check(incub_brms.waal, ndraw = 50)
save(incub_brms.waal, file = "Data_outputs/waal_incub_brms_model.RData")


# ______________________________ ####
# VISUALISE --------------------------------------------------------------------

load("Data_inputs/all_pair_behaviour_incub.RData")

load("Data_outputs/bba_incub_brms_model.RData")
load("Data_outputs/waal_incub_brms_model.RData")

# Datasets
behaviour.bba_incub <- pair_behaviour.incub %>% filter(species == "BBA") 
behaviour.waal_incub <- pair_behaviour.incub %>% filter(species == "WAAL")

### Posterior estimates --------------------------------------------------------

#### BBAL  ---------------------------------------------------------------------
##### Estimates ----------------------------------------------------------------
posterior_samples.incub_bba <- posterior_samples(incub_brms.bba)

# Extract posterior samples by colony
predictors <- colnames(posterior_samples.incub_bba)[2:6]
predictors <- predictors[-2]
interaction_terms <- colnames(posterior_samples.incub_bba)[7:10]
posteriors.incub_bba <- data.frame()

for (i in 1:length(predictors)) {
  
  # Extract the main effect (Bird Island)
  bird_est <- posterior_samples.incub_bba[[predictors[i]]]
  
  # Extract the interaction effect and add it to the main effect for Kerguelen
  kerg_est <- bird_est + posterior_samples.incub_bba[[interaction_terms[i]]]
  
  temp_df <- data.frame(
    variable = predictors[i],
    value = c(bird_est, kerg_est),
    colony = rep(c("Bird Island", "Kerguelen"), each = length(bird_est))
  )
  
  posteriors.incub_bba <- rbind(posteriors.incub_bba, temp_df)

}


## Relevel factors
posteriors.incub_bba %<>%
  mutate(variable = fct_rev(fct_relevel(variable, 
                                        c("b_debt.days", 
                                          "b_median_trip.days", 
                                          "b_diff_trip.days", 
                                          "b_diff_var.days"))))


##### Plot ----------------------------------------------------------------

posteriors_plot.incub_bba.horizontal <- 
  ggplot() +
  stat_halfeye(data = subset(posteriors.incub_bba, colony == "Bird Island"), 
               side = "bottom", 
               aes(x = value, y = variable,
                   fill = ifelse(after_stat(x) < rope_lower, "bi_low", 
                                 ifelse(after_stat(x) > rope_upper, "bi_high", "grey80"))),
               height = 0.7) +
  stat_halfeye(data = subset(posteriors.incub_bba, colony == "Kerguelen"), 
               side = "top", 
               aes(x = value, y = variable,
                   fill = ifelse(after_stat(x) < rope_lower, "ker_low", 
                                 ifelse(after_stat(x) > rope_upper, "ker_high", "grey80"))),
               height = 0.7) +
    scale_fill_manual(values = c(
    "bi_low" = bi_col.low, "bi_high" = bi_col.high,
    "ker_low" = ker_col.low, "ker_high" = ker_col.high ) ) +
  scale_y_discrete(labels = my_labels) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, col = "grey20") +
   geom_rect(aes(xmin = rope_lower, xmax = rope_upper, ymin = -Inf, ymax = Inf),
              fill = "lightcyan4", alpha = 0.25) +
  labs(y = "", 
       x = expression(logit(italic(p[i])))) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none")


#### WAAL  ---------------------------------------------------------------------
##### Estimates ---------------------------------------------------------------------
posterior_samples.incub_waal <- posterior_samples(incub_brms.waal)

# Extract posterior samples by colony
predictors <- colnames(posterior_samples.incub_waal)[2:6]
predictors <- predictors[-2]
interaction_terms <- colnames(posterior_samples.incub_waal)[7:10]
posteriors.incub_waal <- data.frame()

for (i in 1:length(predictors)) {
  
  # Extract the main effect (Bird Island)
  bird_est <- posterior_samples.incub_waal[[predictors[i]]]
  
  # Extract the interaction effect and add it to the main effect for Kerguelen
  cro_est <- bird_est + posterior_samples.incub_waal[[interaction_terms[i]]]
  
  temp_df <- data.frame(
    variable = predictors[i],
    value = c(bird_est, cro_est),
    colony = rep(c("Bird Island", "Crozet"), each = length(bird_est))
  )
  
  posteriors.incub_waal <- rbind(posteriors.incub_waal, temp_df)
  
}


## Relevel factors
posteriors.incub_waal %<>%
  mutate(variable = fct_rev(fct_relevel(variable, 
                                        c("b_debt.days", 
                                          "b_median_trip.days", 
                                          "b_diff_trip.days", 
                                          "b_diff_var.days"))))

##### Plot ---------------------------------------------------------------

posteriors_plot.incub_waal.horizontal <- 
  ggplot() +
  stat_halfeye(data = subset(posteriors.incub_waal, colony == "Bird Island"), 
               side = "bottom", 
               aes(x = value, y = variable,
                   fill = ifelse(after_stat(x) < rope_lower, "bi_low", 
                                 ifelse(after_stat(x) > rope_upper, "bi_high", "grey80"))),
               height = 0.7) +
  stat_halfeye(data = subset(posteriors.incub_waal, colony == "Crozet"), 
               side = "top", 
               aes(x = value, y = variable,
                   fill = ifelse(after_stat(x) < rope_lower, "cro_low", 
                                 ifelse(after_stat(x) > rope_upper, "cro_high", "grey80"))),
               height = 0.7) +
  scale_fill_manual(values = c(
    "bi_low" = bi_col.low, "bi_high" = bi_col.high,
    "cro_low" = ker_col.low, "cro_high" = ker_col.high ) ) +
  scale_y_discrete(labels = my_labels) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, col = "grey20") +
  geom_rect(aes(xmin = rope_lower, xmax = rope_upper, ymin = -Inf, ymax = Inf),
            fill = "lightcyan4", alpha = 0.25) +
  labs(y = "", 
       x = expression(logit(italic(p[i])))) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none")


### Conditional effects  -------------------------------------------------------

effects <- 6:9

#### BBAL  ------------------------------------------------------------
pb <- txtProgressBar(min = 0, max = length(effects), style = 3)

for (i in 1:length(effects)) {
  
  setTxtProgressBar(pb, i)
  
  cond_term <- conditional_effects(incub_brms.bba)[effects[i]]
  var_name <- names(cond_term)[1]
  var_name.data <- gsub(":colony|colony:", "", var_name) 
  var_label <- variable_labels[var_name]
  cond_term.df <- as.data.frame(cond_term[[1]])
 
  cond_plot <- ggplot() +
    geom_ribbon(aes(x = effect1__, ymin = lower__, ymax = upper__, fill = colony), 
                data = cond_term.df, alpha = 0.5) +
    geom_line(aes(x = effect1__, y = estimate__, col = colony), data = cond_term.df, linewidth = 1) +
    geom_point(aes(x = !!sym(var_name.data), y = breeding_outcome.bin, col = colony), data = behaviour.bba_incub) +
    scale_fill_manual(values = c(bi_col, ker_col)) +
    scale_colour_manual(values = c(bi_col, ker_col)) +
    labs(x = var_label, y = expression("P|Breeding success"~(italic(P[i])))) +
    theme_bw() +
    theme(text = element_text(size = 16, family = "Calibri"),
          legend.position = "none")
  
  assign(paste0("cond_plot.", var_name.data, ".bba_incub"), cond_plot)
  
}

close(pb)


#### WAAL ------------------------------------------------------------
pb <- txtProgressBar(min = 0, max = length(effects), style = 3)

for (i in 1:length(effects)) {
  
  setTxtProgressBar(pb, i)
  
  cond_term <- conditional_effects(incub_brms.waal)[effects[i]]
  var_name <- names(cond_term)[1]
  var_name.data <- gsub(":colony|colony:", "", var_name) 
  var_label <- variable_labels[var_name]
  cond_term.df <- as.data.frame(cond_term[[1]])
  
  cond_plot <- ggplot() +
    geom_ribbon(aes(x = effect1__, ymin = lower__, ymax = upper__, fill = colony), 
                data = cond_term.df, alpha = 0.5) +
    geom_line(aes(x = effect1__, y = estimate__, col = colony), data = cond_term.df, linewidth = 1) +
    geom_point(aes(x = !!sym(var_name.data), y = breeding_outcome.bin, col = colony), data = behaviour.waal_incub) +
    scale_fill_manual(values = c(bi_col, ker_col)) +
    scale_colour_manual(values = c(bi_col, ker_col)) +
    labs(x = var_label, expression("P|Breeding success"~(italic(P[i])))) +
    theme_bw() +
    theme(text = element_text(size = 16, family = "Calibri"),
          legend.position = "none")
  
  assign(paste0("cond_plot.", var_name.data, ".waal_incub"), cond_plot)
  
}

close(pb)


# ............................................................ ####
# INTERPRETATION --------------------------------------------------------------

load("Data_inputs/all_pair_behaviour_incub.RData")

load("Data_outputs/bba_incub_brms_model.RData")
load("Data_outputs/waal_incub_brms_model.RData")

# Datasets
behaviour.bba_incub <- pair_behaviour.incub %>% filter(species == "BBA") 
behaviour.waal_incub <- pair_behaviour.incub %>% filter(species == "WAAL")


### BBAL ---------------------------------------------------------------------
summary(incub_brms.bba)
posterior_samples.incub_bba <- posterior_samples(incub_brms.bba)

##### Overall summary table ------------------------------

bbal_incub.fixef <- round(summary(incub_brms.bba)$fixed, digits = 2)
bbal_incub.fixef[,c(6,7)] <- round(bbal_incub.fixef[,c(6,7)], digits = 0)

# Get odds ratios
exp_coef <- round(exp(bbal_incub.fixef[,c(1, 3, 4)]), digits = 2)
colnames(exp_coef) <- c("Estimate_exp", "LCI_exp", "UCI_exp")
rownames(exp_coef) <- NULL

# Get proportion in ROPE
prop_rope <-  data.frame(apply(posterior_samples.incub_bba[,1:10], 2, prop_rope.func))
prop_rope <- prop_rope[,1]

## Bind everything together
bbal_incub.fixef <- cbind(bbal_incub.fixef, exp_coef)
bbal_incub.fixef$prop_rope <- round(prop_rope, digits = 2)

## Process each variable
ker_debt <- process_interaction_estimates("b_debt.days", bbal_incub.fixef, posterior_samples.incub_bba, "b_debt.days:colonykerguelen")
ker_median_trip <- process_interaction_estimates("b_median_trip.days", bbal_incub.fixef, posterior_samples.incub_bba, "b_colonykerguelen:median_trip.days")
ker_diff_trip <- process_interaction_estimates("b_diff_trip.days", bbal_incub.fixef, posterior_samples.incub_bba, "b_colonykerguelen:diff_trip.days")
ker_diff_var <- process_interaction_estimates("b_diff_var.days", bbal_incub.fixef, posterior_samples.incub_bba, "b_colonykerguelen:diff_var.days")

# Combine the results with the original fixef data
ker_fixef <- rbind(ker_debt, ker_median_trip, ker_diff_trip, ker_diff_var)
colnames(ker_fixef) <- colnames(bbal_incub.fixef)

bbal_incub.fixef <- rbind(bbal_incub.fixef, ker_fixef)



### WAAL ---------------------------------------------------------------------
summary(incub_brms.waal)
posterior_samples.incub_waal <- posterior_samples(incub_brms.waal)

### Overall summary table ------------------------------

waal_incub.fixef <- round(summary(incub_brms.waal)$fixed, digits = 2)
waal_incub.fixef[,c(6,7)] <- round(waal_incub.fixef[,c(6,7)], digits = 0)

# Get odds ratios
exp_coef <- round(exp(waal_incub.fixef[,c(1, 3, 4)]), digits = 2)
colnames(exp_coef) <- c("Estimate_exp", "LCI_exp", "UCI_exp")
rownames(exp_coef) <- NULL

# Get proportion in ROPE
prop_rope <-  data.frame(apply(posterior_samples.incub_waal[,1:10], 2, prop_rope.func))
prop_rope <- prop_rope[,1]

## Bind everything together
waal_incub.fixef <- cbind(waal_incub.fixef, exp_coef)
waal_incub.fixef$prop_rope <- round(prop_rope, digits = 2)

### Process each variable
cro_debt <- process_interaction_estimates("b_debt.days", waal_incub.fixef, posterior_samples.incub_waal, "b_debt.days:colonycrozet")
cro_median_trip <- process_interaction_estimates("b_median_trip.days", waal_incub.fixef, posterior_samples.incub_waal, "b_colonycrozet:median_trip.days")
cro_diff_trip <- process_interaction_estimates("b_diff_trip.days", waal_incub.fixef, posterior_samples.incub_waal, "b_colonycrozet:diff_trip.days")
cro_diff_var <- process_interaction_estimates("b_diff_var.days", waal_incub.fixef, posterior_samples.incub_waal, "b_colonycrozet:diff_var.days")

# Combine the results with the original fixef data
cro_fixef <- rbind(cro_debt, cro_median_trip, cro_diff_trip, cro_diff_var)
colnames(cro_fixef) <- colnames(waal_incub.fixef)

waal_incub.fixef <- rbind(waal_incub.fixef, cro_fixef)




# ______________________________ ####
# ~ * BROODING * ~ ###########################################################

# ANALYSIS =====================================================================

load("Data_inputs/all_pair_behaviour_brooding.RData")

# Datasets
behaviour.bba_brooding <- pair_behaviour.brooding %>% filter(species == "BBA") 
behaviour.waal_brooding <- pair_behaviour.brooding %>% filter(species == "WAAL")

#### Examine correlations ---------------------------------------------------------------

# BBA
GGally::ggpairs(
  data = behaviour.bba_brooding %>%
    select(colony, debt.days, median_trip.days, diff_trip.days, diff_var.days),
  columns = 2:5,
  mapping = ggplot2::aes(colour = colony),
  title = "Behavioural correlations by colony") +
  scale_colour_manual(values = c(bi_col, ker_col)) +
  scale_fill_manual(values = c(bi_col, ker_col))


# WAAL
GGally::ggpairs(
  data = behaviour.waal_brooding %>%
    select(colony, debt.days, median_trip.days, diff_trip.days, diff_var.days),
  columns = 2:5,
  mapping = ggplot2::aes(colour = colony),
  title = "Behavioural correlations by colony") +
  scale_colour_manual(values = c(bi_col, ker_col)) +
  scale_fill_manual(values = c(bi_col, ker_col))


#### Specify priors ---------------------------------------------------------------

# BBA 
priors.bba <- c( set_prior("normal(-0.5, 1)", class = "b", coef = "debt.days"),
                 set_prior("normal(-0.5, 0.5)", class = "b", coef = "median_trip.days"),       
                 set_prior("normal(-0.5, 0.5)", class = "b", coef = "diff_trip.days"),
                 set_prior("normal(-0.12, 0.5)", class = "b", coef = "diff_var.days"))

# WAAL
priors.waal <- c( set_prior("normal(-0.5, 1)", class = "b", coef = "debt.days"),
                  set_prior("normal(-0.5, 0.5)", class = "b", coef = "median_trip.days"),       
                  set_prior("normal(-0.5, 0.5)", class = "b", coef = "diff_trip.days"),
                  set_prior("normal(-0.12, 0.5)", class = "b", coef = "diff_var.days"))

#### Specify models ----------------------------------------------------------------

# For Kerguelen, I have very limited data, so I am using a partial pooling approach 

# Model structure
bf.brooding.partialpool <- brms::bf(
  breeding_outcome.bin ~ 1 +
    (1 + debt.days + median_trip.days + diff_trip.days + diff_var.days | colony) +
    (1 | season),
  family = bernoulli() )

brooding_brms.bba <- brm(bf.brooding.partialpool,
                         data = behaviour.bba_brooding, 
                         cores = 4, chains = 4, 
                         iter = 2000, warmup = 1000, thin = 10,
                         control = list(adapt_delta = 0.9999, max_treedepth = 14))

pp_check(brooding_brms.bba)
save(brooding_brms.bba, file = "Data_outputs/bba_brooding_brms_model.RData")



##### TESTING



####

pp_check(brooding_brms.bba)
save(brooding_brms.bba, file = "Data_outputs/bba_brooding_brms_model.RData")

# WAAL

bf.brooding <- brms::bf(breeding_outcome.bin ~
                          debt.days * colony +
                          median_trip.days * colony +
                          diff_trip.days * colony +
                          diff_var.days * colony +
                          (1|season), 
                        family = "bernoulli")

brooding_brms.waal <- brm(bf.brooding,
                       data = behaviour.waal_brooding, 
                       cores = 4, chains = 4, 
                       iter = 20000, warmup = 10000, thin = 10,
                       control = list(adapt_delta = 0.9999, max_treedepth = 14),
                       prior = priors.waal)

summary(brooding_brms.waal)
pp_check(brooding_brms.waal)
save(brooding_brms.waal, file = "Data_outputs/waal_brooding_brms_model.RData")

# ______________________________ ####
# VISUALISE --------------------------------------------------------------------

load("Data_inputs/all_pair_behaviour_brooding.RData")

load("Data_outputs/bba_brooding_brms_model.RData")
load("Data_outputs/waal_brooding_brms_model.RData")

# Datasets
behaviour.bba_brooding <- pair_behaviour.brooding %>% filter(species == "BBA") 
behaviour.waal_brooding <- pair_behaviour.brooding %>% filter(species == "WAAL")

### Posterior estimates --------------------------------------------------------

#### BBAL  ---------------------------------------------------------------------
##### Estimates ----------------------------------------------------------------

posterior_samples.brooding_bba <- as_draws_df(brooding_brms_partial.bba)

# Random effects for each colony and predictor
colony_effects <- posterior_samples.brooding_bba[, grep("^r_", names(posterior_samples.brooding_bba))]
posteriors.brooding_bba <- reshape2::melt(colony_effects) %>% filter(!grepl("Intercept", variable)) 

posteriors.brooding_bba %<>%
  mutate(colony = gsub("r_colony\\[([^,]+),.*\\]", "\\1", variable),  
    effect = gsub("r_colony\\[[^,]+,(.*)\\]", "\\1", variable)) %>%
  select(-variable) %>%
  rename(variable = effect) %>%
  relocate(variable, value, colony) %>%
  mutate(colony = ifelse(colony == "birdisland", "Bird Island", "Kerguelen"))

# Relevel factors
posteriors.brooding_bba %<>%
  mutate(variable = fct_rev(fct_relevel(variable, 
                                        c("debt.days", 
                                          "median_trip.days", 
                                          "diff_trip.days", 
                                          "diff_var.days"))))

# posterior_samples.brooding_bba <- posterior_samples(brooding_brms.bba)
# 
# # Extract posterior samples by colony
# predictors <- colnames(posterior_samples.brooding_bba)[2:6]
# predictors <- predictors[-2]
# interaction_terms <- colnames(posterior_samples.brooding_bba)[7:10]
# posteriors.brooding_bba <- data.frame()
# 
# for (i in 1:length(predictors)) {
#   
#   # Extract the main effect (Bird Island)
#   bird_est <- posterior_samples.brooding_bba[[predictors[i]]]
#   
#   # Extract the interaction effect and add it to the main effect for Kerguelen
#   kerg_est <- bird_est + posterior_samples.brooding_bba[[interaction_terms[i]]]
#   
#   temp_df <- data.frame(
#     variable = predictors[i],
#     value = c(bird_est, kerg_est),
#     colony = rep(c("Bird Island", "Kerguelen"), each = length(bird_est))
#   )
#   
#   posteriors.brooding_bba <- rbind(posteriors.brooding_bba, temp_df)
#   
# }
# 
# ## Relevel factors
# posteriors.brooding_bba %<>%
#   mutate(variable = fct_rev(fct_relevel(variable, 
#                                         c("b_debt.days", 
#                                           "b_median_trip.days", 
#                                           "b_diff_trip.days", 
#                                           "b_diff_var.days"))))


##### Plot ----------------------------------------------------------------


posteriors_plot.brooding_bba.horizontal <- 
  ggplot() +
  stat_halfeye(data = subset(posteriors.brooding_bba, colony == "Bird Island"), 
               side = "bottom", 
               aes(x = value, y = variable,
                   fill = ifelse(after_stat(x) < rope_lower, "bi_low", 
                                 ifelse(after_stat(x) > rope_upper, "bi_high", "grey80"))),
               height = 0.7) +
  stat_halfeye(data = subset(posteriors.brooding_bba, colony == "Kerguelen"), 
               side = "top", 
               aes(x = value, y = variable,
                   fill = ifelse(after_stat(x) < rope_lower, "ker_low", 
                                 ifelse(after_stat(x) > rope_upper, "ker_high", "grey80"))),
               height = 0.7) +
  scale_fill_manual(values = c(
    "bi_low" = bi_col.low, "bi_high" = bi_col.high,
    "ker_low" = ker_col.low, "ker_high" = ker_col.high ) ) +
  scale_y_discrete(labels = my_labels) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, col = "grey20") +
  geom_rect(aes(xmin = rope_lower, xmax = rope_upper, ymin = -Inf, ymax = Inf),
            fill = "lightcyan4", alpha = 0.25) +
  labs(y = "", 
       x = expression(logit(italic(p[i])))) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none")


# posteriors_plot.brooding_bba.horizontal <- 
#   ggplot() +
#   stat_halfeye(data = subset(posteriors.brooding_bba, colony == "Bird Island"), 
#                side = "bottom", 
#                aes(x = value, y = variable,
#                    fill = ifelse(after_stat(x) < rope_lower, "bi_low", 
#                                  ifelse(after_stat(x) > rope_upper, "bi_high", "grey80"))),
#                height = 0.7) +
#   stat_halfeye(data = subset(posteriors.brooding_bba, colony == "Kerguelen"), 
#                side = "top", 
#                aes(x = value, y = variable,
#                    fill = ifelse(after_stat(x) < rope_lower, "ker_low", 
#                                  ifelse(after_stat(x) > rope_upper, "ker_high", "grey80"))),
#                height = 0.7) +
#   scale_fill_manual(values = c(
#     "bi_low" = bi_col.low, "bi_high" = bi_col.high,
#     "ker_low" = ker_col.low, "ker_high" = ker_col.high ) ) +
#   scale_y_discrete(labels = my_labels) +
#   geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, col = "grey20") +
#   geom_rect(aes(xmin = rope_lower, xmax = rope_upper, ymin = -Inf, ymax = Inf),
#             fill = "lightcyan4", alpha = 0.25) +
#   labs(y = "", 
#        x = expression(logit(italic(p[i])))) +
#   theme_bw() +
#   theme(text = element_text(size = 16, family = "Calibri"),
#         legend.position = "none")




#### WAAL ----------------------------------------------------------------------
##### Estimates ----------------------------------------------------------------
posterior_samples.brooding_waal <- posterior_samples(brooding_brms.waal)

# Extract posterior samples by colony
predictors <- colnames(posterior_samples.brooding_waal)[2:6]
predictors <- predictors[-2]
interaction_terms <- colnames(posterior_samples.brooding_waal)[7:11]
posteriors.brooding_waal <- data.frame()

for (i in 1:length(predictors)) {
  
  # Extract the main effect (Bird Island)
  bird_est <- posterior_samples.brooding_waal[[predictors[i]]]
  
  # Extract the interaction effect and add it to the main effect for Kerguelen
  cro_est <- bird_est + posterior_samples.brooding_waal[[interaction_terms[i]]]
  
  temp_df <- data.frame(
    variable = predictors[i],
    value = c(bird_est, cro_est),
    colony = rep(c("Bird Island", "Crozet"), each = length(bird_est))
  )
  
  posteriors.brooding_waal <- rbind(posteriors.brooding_waal, temp_df)
  
}


## Relevel factors
posteriors.brooding_waal %<>%
  mutate(variable = fct_rev(fct_relevel(variable, 
                                        c("b_debt.days", 
                                          "b_median_trip.days", 
                                          "b_diff_trip.days", 
                                          "b_diff_var.days"))))

##### Plots ---------------------------------------------------------------------

posteriors_plot.brooding_waal.horizontal <- 
  ggplot() +
  stat_halfeye(data = subset(posteriors.brooding_waal, colony == "Bird Island"), 
               side = "bottom", 
               aes(x = value, y = variable,
                   fill = ifelse(after_stat(x) < rope_lower, "bi_low", 
                                 ifelse(after_stat(x) > rope_upper, "bi_high", "grey80"))),
               height = 0.7) +
  stat_halfeye(data = subset(posteriors.brooding_waal, colony == "Crozet"), 
               side = "top", 
               aes(x = value, y = variable,
                   fill = ifelse(after_stat(x) < rope_lower, "cro_low", 
                                 ifelse(after_stat(x) > rope_upper, "cro_high", "grey80"))),
               height = 0.7) +
  scale_fill_manual(values = c(
    "bi_low" = bi_col.low, "bi_high" = bi_col.high,
    "cro_low" = ker_col.low, "cro_high" = ker_col.high ) ) +
  scale_y_discrete(labels = my_labels) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, col = "grey20") +
  geom_rect(aes(xmin = rope_lower, xmax = rope_upper, ymin = -Inf, ymax = Inf),
            fill = "lightcyan4", alpha = 0.25) +
  labs(y = "", 
       x = expression(logit(italic(p[i])))) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none")


## Conditional effects  --------------------------------------------------------

#### BBAL  ---------------------------------------------------------------------

pb <- txtProgressBar(min = 0, max = length(effects), style = 3)

effects <- c("debt.days", "median_trip.days", "diff_trip.days", "diff_var.days")

for (i in 1:length(effects)) {
  
  setTxtProgressBar(pb, i)
  
  var_name <- effects[i]
  var_label <- variable_labels2[var_name]
 
  # Make the prediction DF
  fixed_values <- list(
    colony = c("kerguelen", "birdisland"),
    season = NA)
  
  # Dynamically set the values for the non-varying variables (set to 0, unless it's the one being predicted)
  for (fixed_var in c("debt.days", "diff_trip.days", "diff_var.days", "median_trip.days")) {
    if (fixed_var != var_name) {
      fixed_values[[fixed_var]] <- 0  # Set other variables to 0 if they are not the one being predicted
    } else {
      # If the current variable is the one being predicted, leave it out of the fixed values
      fixed_values[[fixed_var]] <- NULL
    }
  }
  
  ## Create a sequence for the variable currently being examined
  effect_seq <- seq(
    min(behaviour.bba_brooding[[var_name]], na.rm = TRUE),
    max(behaviour.bba_brooding[[var_name]], na.rm = TRUE),
    length.out = 100
  )
  
  ## Create a list of all variables to pass to expand.grid
  grid_list <- list(fixed_values)[[1]]
  
  # Overwrite the current variable with its actual sequence
  grid_list[[var_name]] <- effect_seq
  
  # Generate new data
  newdat <- expand.grid(grid_list) %>% arrange(colony)
  
  # Predict from model, including random effects
  preds <- fitted(
    brooding_brms_partial.bba,
    newdata = newdat,
    re_formula = NULL,    
    summary = TRUE,       
    scale = "response")
  
  # Combine predictions with the new data
  newdat$pred <- preds[, "Estimate"]
  newdat$lower <- preds[, "Q2.5"]
  newdat$upper <- preds[, "Q97.5"]
  
  cond_plot <- ggplot() +
    # geom_ribbon(aes(x = diff_trip.days, y = pred, ymin = lower, ymax = upper, fill = colony), 
    #             data = newdat, alpha = 0.5) +
    geom_line(aes(x = diff_trip.days, y = pred, colour = colony), data = newdat, linewidth = 1) +
    geom_point(aes(x = !!sym(var_name), y = breeding_outcome.bin, col = colony), data = behaviour.bba_brooding) +
    scale_fill_manual(values = c(ker_col, bi_col)) +
    scale_colour_manual(values = c(ker_col, bi_col)) +
    labs(x = var_label, y = expression("P|Breeding success"~(italic(P[i])))) +
    theme_bw() +
    theme(text = element_text(size = 16, family = "Calibri"),
          legend.position = "none")
  

}

effects <- 6:9


## BBA ##
pb <- txtProgressBar(min = 0, max = length(effects), style = 3)

for (i in 1:length(effects)) {
  
  setTxtProgressBar(pb, i)
  
  cond_term <- conditional_effects(brooding_brms.bba)[effects[i]]
  var_name <- names(cond_term)[1]
  var_name.data <- gsub(":colony|colony:", "", var_name) 
  var_label <- variable_labels[var_name]
  cond_term.df <- as.data.frame(cond_term[[1]])
  
  cond_plot <- ggplot() +
    geom_ribbon(aes(x = effect1__, ymin = lower__, ymax = upper__, fill = colony), 
                data = cond_term.df, alpha = 0.5) +
    geom_line(aes(x = effect1__, y = estimate__, col = colony), data = cond_term.df, linewidth = 1) +
    geom_point(aes(x = !!sym(var_name.data), y = breeding_outcome.bin, col = colony), data = behaviour.bba_brooding) +
    scale_fill_manual(values = c(bi_col, ker_col)) +
    scale_colour_manual(values = c(bi_col, ker_col)) +
    labs(x = var_label, y = expression("P|Breeding success"~(italic(P[i])))) +
    theme_bw() +
    theme(text = element_text(size = 16, family = "Calibri"),
          legend.position = "none")
  
  assign(paste0("cond_plot.", var_name.data, ".bba_brooding"), cond_plot)
  
}

close(pb)


## WAAL ##
pb <- txtProgressBar(min = 0, max = length(effects), style = 3)

for (i in 1:length(effects)) {
  
  setTxtProgressBar(pb, i)
  
  cond_term <- conditional_effects(brooding_brms.waal)[effects[i]]
  var_name <- names(cond_term)[1]
  var_name.data <- gsub(":colony|colony:", "", var_name) 
  var_label <- variable_labels[var_name]
  cond_term.df <- as.data.frame(cond_term[[1]])
  
  cond_plot <- ggplot() +
    geom_ribbon(aes(x = effect1__, ymin = lower__, ymax = upper__, fill = colony), 
                data = cond_term.df, alpha = 0.5) +
    geom_line(aes(x = effect1__, y = estimate__, col = colony), data = cond_term.df, linewidth = 1) +
    geom_point(aes(x = !!sym(var_name.data), y = breeding_outcome.bin, col = colony), data = behaviour.waal_brooding) +
    scale_fill_manual(values = c(bi_col, ker_col)) +
    scale_colour_manual(values = c(bi_col, ker_col)) +
    labs(x = var_label, y = expression("P|Breeding success"~(italic(P[i])))) +
    theme_bw() +
    theme(text = element_text(size = 16, family = "Calibri"),
          legend.position = "none")
  
  assign(paste0("cond_plot.", var_name.data, ".waal_brooding"), cond_plot)
 
}

close(pb)


# ______________________________ ####
# INTERPRETATION --------------------------------------------------------------

load("Data_inputs/all_pair_behaviour_brooding.RData")

load("Data_outputs/bba_brooding_brms_model.RData")
load("Data_outputs/waal_brooding_brms_model.RData")

### BBAL ---------------------------------------------------------------------

summary(brooding_brms.bba)
posterior_samples.brooding_bba <- posterior_samples(brooding_brms.bba)


### TEMP

## Get effects ##

mean(post[["r_colony[kerguelen,median_trip.days]"]])
exp(mean(post[["r_colony[kerguelen,median_trip.days]"]]))
exp(quantile(post[["r_colony[kerguelen,median_trip.days]"]]))

exp(mean(post[["r_colony[birdisland,median_trip.days]"]]))
exp(quantile(post[["r_colony[birdisland,median_trip.days]"]]))

###

##### Overall summary table ------------------------------

bbal_brooding.fixef <- round(summary(brooding_brms.bba)$fixed, digits = 2)
bbal_brooding.fixef[,c(6,7)] <- round(bbal_brooding.fixef[,c(6,7)], digits = 0)

# Get odds ratios
exp_coef <- round(exp(bbal_brooding.fixef[,c(1, 3, 4)]), digits = 2)
colnames(exp_coef) <- c("Estimate_exp", "LCI_exp", "UCI_exp")
rownames(exp_coef) <- NULL

# Get proportion in ROPE
prop_rope <-  data.frame(apply(posterior_samples.brooding_bba[,1:10], 2, prop_rope.func))
prop_rope <- prop_rope[,1]

## Bind everything together
bbal_brooding.fixef <- cbind(bbal_brooding.fixef, exp_coef)
bbal_brooding.fixef$prop_rope <- round(prop_rope, digits = 2)

## Process each variable
ker_debt <- process_interaction_estimates("b_debt.days", bbal_brooding.fixef, posterior_samples.brooding_bba, "b_debt.days:colonykerguelen")
ker_median_trip <- process_interaction_estimates("b_median_trip.days", bbal_brooding.fixef, posterior_samples.brooding_bba, "b_colonykerguelen:median_trip.days")
ker_diff_trip <- process_interaction_estimates("b_diff_trip.days", bbal_brooding.fixef, posterior_samples.brooding_bba, "b_colonykerguelen:diff_trip.days")
ker_diff_var <- process_interaction_estimates("b_diff_var.days", bbal_brooding.fixef, posterior_samples.brooding_bba, "b_colonykerguelen:diff_var.days")

# Combine the results with the original fixef data
ker_fixef <- rbind(ker_debt, ker_median_trip, ker_diff_trip, ker_diff_var)
colnames(ker_fixef) <- colnames(bbal_brooding.fixef)

bbal_brooding.fixef <- rbind(bbal_brooding.fixef, ker_fixef)



### WAAL ---------------------------------------------------------------------
summary(brooding_brms.waal)
posterior_samples.brooding_waal <- posterior_samples(brooding_brms.waal)

### Overall summary table ------------------------------

waal_brooding.fixef <- round(summary(brooding_brms.waal)$fixed, digits = 2)
waal_brooding.fixef[,c(6,7)] <- round(waal_brooding.fixef[,c(6,7)], digits = 0)

# Get odds ratios
exp_coef <- round(exp(waal_brooding.fixef[,c(1, 3, 4)]), digits = 2)
colnames(exp_coef) <- c("Estimate_exp", "LCI_exp", "UCI_exp")
rownames(exp_coef) <- NULL

# Get proportion in ROPE
prop_rope <-  data.frame(apply(posterior_samples.brooding_waal[,1:10], 2, prop_rope.func))
prop_rope <- prop_rope[,1]

## Bind everything together
waal_brooding.fixef <- cbind(waal_brooding.fixef, exp_coef)
waal_brooding.fixef$prop_rope <- round(prop_rope, digits = 2)

### Process each variable
cro_debt <- process_interaction_estimates("b_debt.days", waal_brooding.fixef, posterior_samples.brooding_waal, "b_debt.days:colonycrozet")
cro_median_trip <- process_interaction_estimates("b_median_trip.days", waal_brooding.fixef, posterior_samples.brooding_waal, "b_colonycrozet:median_trip.days")
cro_diff_trip <- process_interaction_estimates("b_diff_trip.days", waal_brooding.fixef, posterior_samples.brooding_waal, "b_colonycrozet:diff_trip.days")
cro_diff_var <- process_interaction_estimates("b_diff_var.days", waal_brooding.fixef, posterior_samples.brooding_waal, "b_colonycrozet:diff_var.days")

# Combine the results with the original fixef data
cro_fixef <- rbind(cro_debt, cro_median_trip, cro_diff_trip, cro_diff_var)
colnames(cro_fixef) <- colnames(waal_brooding.fixef)

waal_brooding.fixef <- rbind(waal_brooding.fixef, cro_fixef)




# ______________________________ ####
# * FIGURE 1 * COMBINED POSTERIOR PLOTS ========================================

# Add albatross silhouettes
posteriors_plot.incub_bba.horizontal2<- ggdraw() +
  draw_plot(posteriors_plot.incub_bba.horizontal + labs(y = "Incubation") +
              annotate("text", x = 0.9, y = 4.5, label = "Increases success") +
              annotate("text", x = -0.6, y = 4.5, label = "Decreases success") +
                theme(axis.title.y = element_text(margin = margin(r = 20)),
                    axis.title.x = element_blank(),
                    text = element_text(size = 16, family = "Calibri"))) +
  draw_image(file.path("Figures/bba_standing_silhouette.png"),
             scale = 0.25, x = 0.35, y = -0.28) 
  

posteriors_plot.incub_waal.horizontal2 <- ggdraw() +
  draw_plot(posteriors_plot.incub_waal.horizontal +
              annotate("text", x = 0.5, y = 4.5, label = "Increases success") +
              annotate("text", x = -0.4, y = 4.5, label = "Decreases success") +
              theme(axis.title.x = element_blank(),
                    axis.text.y = element_blank(),
                    text = element_text(size = 16, family = "Calibri"))) +
  draw_image(file.path("Figures/waal_standing_silhouette.png"),
             scale = 0.25, x = 0.38, y = -0.28)


png(file = "Figures/FIGURE1.png", width = 10, height = 10, units = "in", res = 100)
ggarrange(posteriors_plot.incub_bba.horizontal2,
          posteriors_plot.incub_waal.horizontal2,
          posteriors_plot.brooding_bba.horizontal + 
            labs(y = "Brooding") + 
            theme(axis.title.y = element_text(margin = margin(r = 20)),
                  plot.title = element_blank()),
          posteriors_plot.brooding_waal.horizontal + 
            theme(plot.title = element_blank(),
                  axis.text.y = element_blank()),
          ncol = 2,
          nrow = 2,
          widths = c(1, 0.95))
dev.off()


# * FIGURE 2 * BBAL CONDITIONAL PLOTS ==========================================

cond_plot.median_trip.days.bba_incub2 <- ggdraw() +
  draw_plot(cond_plot.median_trip.days.bba_incub + labs(tag = "A",
                                                        y = "Incubation \nP|Breeding success") +
              theme(axis.title.y = element_text(margin = margin(r = 35)),
                    plot.tag = element_text(size = 16, face = "bold"),
                    plot.tag.position = c(0.05, 1.01),
                    plot.margin = unit(c(0.6, 0.2, 0.1, 0.5), "cm"),
                    text = element_text(size = 16, family = "Calibri"))) +
  draw_image(file.path("Figures/bba_standing_silhouette.png"),
             scale = 0.25, x = 0.32, y = -0.23) 

png(file = "Figures/FIGURE2.png", width = 10, height = 10, units = "in", res = 300)
ggarrange(
  cond_plot.median_trip.days.bba_incub2,
  NULL,
  cond_plot.median_trip.days.bba_brooding + 
    labs(tag = "B", y = "Brooding \nP|Breeding success") + 
    xlim(0.5, 5.8) +
    theme(plot.tag = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(margin = margin(r = 35)),
          plot.tag.position = c(0.05, 1.01),
          plot.margin = unit(c(0.6, 0.2, 0.1, 0.5), "cm")),
  cond_plot.diff_trip.days.bba_brooding + 
    labs(tag = "C") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag = element_text(size = 16, face = "bold"),
          plot.margin = unit(c(0.3, 0.2, 0.1, 0.5), "cm")),
  ncol = 2,
  nrow = 2,
  widths = c(1, 0.85)
)
dev.off()


# * FIGURE 3 * WAAL CONDITIONAL PLOTS ==========================================

# Add albatross silhouette
cond_plot.diff_var.days.waal_incub2 <- ggdraw() +
  draw_plot(cond_plot.diff_var.days.waal_incub +
              labs(tag = "A", y = "Incubation \nP|Breeding success") +
              theme(axis.title.y = element_text(margin = margin(r = 35)),
                    plot.tag = element_text(size = 16, face = "bold"),
                    plot.tag.position = c(0.05, 1),
                    plot.margin = unit(c(0.6, 0.2, 0.1, 0.5), "cm"),
                    text = element_text(size = 16, family = "Calibri"))) +
  draw_image(file.path("Figures/waal_standing_silhouette.png"),
             scale = 0.3, x = -0.12, y = -0.18)


png(file = "Figures/FIGURE3.png", width = 10, height = 10, units = "in", res = 300)
ggarrange(
  cond_plot.diff_var.days.waal_incub2,
  NULL,
  cond_plot.median_trip.days.waal_brooding + 
    labs(tag = "B", y = "Brooding \nP|Breeding success") + 
    theme(plot.tag = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(margin = margin(r = 35)),
           plot.tag.position = c(0.05, 1),
          plot.margin = unit(c(0.6, 0.2, 0.1, 0.5), "cm")),
  cond_plot.diff_trip.days.waal_brooding + 
    labs(tag = "C") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag = element_text(size = 16, face = "bold"),
          plot.margin = unit(c(0.3, 0.2, 0.1, 0.5), "cm")),
  ncol = 2,
  nrow = 2,
  widths = c(1, 0.85)
)
dev.off()

# * TABLE 1 * BBAL COEFFICIENTS ================================================

summary(incub_brms.bba); summary(brooding_brms.bba)

bbal_incub.fixef_MS <- bbal_incub.fixef %>%
  mutate(est_error = paste0(Estimate, " \u00B1 ", Est.Error),
         cis = paste0('[', `l-95% CI`, ", ", `u-95% CI`, ']')) %>%
  select(est_error, cis, below_zero, above_zero)

bbal_brooding.fixef_MS <- bbal_brooding.fixef %>%
  mutate(est_error = paste0(Estimate, " \u00B1 ", Est.Error),
         cis = paste0('[', `l-95% CI`, ", ", `u-95% CI`, ']')) %>%
  select(est_error, cis, below_zero, above_zero)

table1.csv <- cbind(bbal_incub.fixef_MS, bbal_brooding.fixef_MS)
readr::write_excel_csv(table1.csv, "Data_outputs/bbal_coefficients.csv")


# * TABLE 2 * BBAL COEFFICIENTS ================================================

summary(incub_brms.waal); summary(brooding_brms.waal)

waal_incub.fixef_MS <- waal_incub.fixef %>%
  mutate(est_error = paste0(Estimate, " ± ", Est.Error),
         cis = paste0('[', `l-95% CI`, ", ", `u-95% CI`, ']')) %>%
  select(est_error, cis, below_zero, above_zero)

waal_brooding.fixef_MS <- waal_brooding.fixef %>%
  mutate(est_error = paste0(Estimate, " ± ", Est.Error),
         cis = paste0('[', `l-95% CI`, ", ", `u-95% CI`, ']')) %>%
  select(est_error, cis, below_zero, above_zero)

table2.csv <- cbind(waal_incub.fixef_MS, waal_brooding.fixef_MS)
readr::write_excel_csv(table2.csv, "Data_outputs/waal_coefficients.csv")


