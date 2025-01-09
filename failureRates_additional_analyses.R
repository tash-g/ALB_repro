## ---------------------------
##
## Script name: failureRates_additional_analyses.R
##
## Purpose of script: Supplementary analyses for the failure rates project
##
## Author: Dr. Natasha Gillies
##
## Created: 2024-08-10
##
## Email: gilliesne@gmail.com
##
## ---------------------------


# Load functions, packages, & data ---------------------------------------------

# Functions for GPS processing and plotting
source("ALB_FOR_functions.R")

# Define the packages
packages <- c("dplyr", "magrittr", "ggplot2", "brms", "tidyverse", "tidybayes",
              "ggpubr", "extrafont", "lme4", "lmerTest", "cowplot")

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())
# 
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Suppress dplyr summarise warning
options(dplyr.summarise.inform = FALSE)

# Make sure using dplyr select
select <- dplyr::select


### Set colour parameters ------------------------------------------------------

incub_col <- "#4682B4"
brooding_col <- "#FFB347"
brooding_col.low <- "#ffd69c"
brooding_col.high <- "#fc980a"

ker_col <-  "#156082" 
cro_col <- "#0F9ED5"
bi_col <- "#E97132" 

fail_chick <- "#ED7A6E"
fail_egg <- "#FFA07A"
fail_unknown <- "#D3D3D3"

# LOAD DATA ==========================================================

## BBA - BIRD ISLAND (baBI) ##
load("Data_inputs/bba_birdis_breedingTrips.RData")
baBI_trips <- all_trips %>% distinct(); rm(all_trips)

# Remove implausible trips
baBI_trips %<>% filter(duration.mins > 360 & duration.days < 31)

# Remove birds where sex or breeding outcome unknown
baBI_trips <- subset(baBI_trips, !is.na(rs))
baBI_trips.sex <- subset(baBI_trips, !is.na(sex))


## BBA - KERGUELEN (baKer) ##
load("Data_inputs/bba_kerguelen_breedingTrips.RData")
baKer_trips <- all_trips %>% distinct(); rm(all_trips)

# Remove implausible trips
baKer_trips %<>% filter(duration.mins > 360 & duration.days < 31)

# Remove birds where sex or breeding outcome unknown
baKer_trips <- subset(baKer_trips, !is.na(rs) & rs != "UNKNOWN")
baKer_trips.sex <- subset(baKer_trips, !is.na(sex))


## WAAL - BIRD ISLAND (waBI) ##
load("Data_inputs/waal_birdis_breedingTrips.RData")
waBI_trips <- all_trips %>% distinct(); rm(all_trips)

# Remove implausible trips
waBI_trips %<>% filter(duration.mins > 360 & duration.days < 31)

# Remove birds where sex or breeding outcome unknown
waBI_trips <- subset(waBI_trips, !is.na(rs))
waBI_trips.sex <- subset(waBI_trips, !is.na(sex))



## WAAL - CROZET (waCro) ##
load("Data_inputs/waal_crozet_breedingTrips.RData")
waCro_trips <- all_trips %>% distinct(); rm(all_trips)

# Remove implausible trips
waCro_trips %<>% filter(duration.mins > 360 & duration.days < 31)

# Remove birds where sex or breeding outcome unknown
waCro_trips <- subset(waCro_trips, !is.na(rs) & rs != "UNKNOWN" & rs != "NON-BREEDER")
waCro_trips.sex <- subset(waCro_trips, !is.na(sex))


# Just Bird Island
bas_trips <- rbind(baBI_trips, waBI_trips)

# ______________________________ ####
# ~ * SEX ANALYSIS - INCUBATION * ~ ###########################################################

# ............................................................ ####
# ANALYSIS =====================================================================

load("Data_inputs/all_pair_behaviour_incub_withSex.RData")

#### Specify priors ---------------------------------------------------------------

# BBA 
priors.bba <- c(set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.F"),
                set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.M"),
                set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.M"),
                set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.F"),
                set_prior("normal(-0.12, 0.5)", class = "b", coef = "debt.days")#,
                #set_prior("normal(-0.5, 1)", class = "b", coef = "male_female_diff.var")
                )

# WAAL
priors.waal <- c(set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.F"),
                 set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.M"),
                 set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.M"),
                 set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.F"),
                 set_prior("normal(-0.12, 0.5)", class = "b", coef = "debt.days")#,
                 #set_prior("normal(-0.5, 1)", class = "b", coef = "male_female_diff.var")
                )

#### Specify models ----------------------------------------------------------------

# Datasets
behaviour.bba_incub <- pair_behaviour.incub_withSex %>% filter(species == "BBA") 
behaviour.waal_incub <- pair_behaviour.incub_withSex %>% filter(species == "WAAL")

# Model structure
bf.incub <- brms::bf(breeding_outcome.bin ~
                       trip_days.F * colony +
                       trip_days.M * colony +
                       sd_trip.F * colony +
                       sd_trip.M * colony + 
                       debt.days * colony +
                       (1|season), 
                     family = "bernoulli")


## Fit models
incub_brms.bba.sex <- brm(bf.incub,
                      data = behaviour.bba_incub, 
                      cores = 4, chains = 4, 
                      iter = 20000, warmup = 10000, thin = 10,
                      control = list(adapt_delta = 0.9999, max_treedepth = 14),
                      prior = priors.bba)

summary(incub_brms.bba.sex)
pp_check(incub_brms.bba.sex, ndraw = 50)
save(incub_brms.bba.sex, file = "Data_outputs/bba_incub_brms_model_withSex.RData")


incub_brms.waal.sex <- brm(bf.incub,
                       data = behaviour.waal_incub, 
                       cores = 4, chains = 4, 
                       iter = 20000, warmup = 10000, thin = 10,
                       control = list(adapt_delta = 0.9999, max_treedepth = 14),
                       prior = priors.waal)

summary(incub_brms.waal.sex)
pp_check(incub_brms.waal.sex, ndraw = 50)
save(incub_brms.waal.sex, file = "Data_outputs/waal_incub_brms_model_withSex.RData")

# ............................................................ ####
# VISUALISE --------------------------------------------------------------------

load("Data_inputs/all_pair_behaviour_incub_withSex.RData")

load("Data_outputs/bba_incub_brms_model_withSex.RData")
load("Data_outputs/waal_incub_brms_model_withSex.RData")

# Datasets
behaviour.bba_incub <- pair_behaviour.incub_withSex %>% filter(species == "BBA") 
behaviour.waal_incub <- pair_behaviour.incub_withSex %>% filter(species == "WAAL")

### Posterior estimates --------------------------------------------------------

#### BBAL  ---------------------------------------------------------------------
##### Estimates ----------------------------------------------------------------
posterior_samples.incub_bba <- posterior_samples(incub_brms.bba.sex)

# Extract posterior samples by colony
predictors <- colnames(posterior_samples.incub_bba)[2:7]
predictors <- predictors[-2]
interaction_terms <- colnames(posterior_samples.incub_bba)[8:12]
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


## Rename and relevel factors
posteriors.incub_bba %<>%
  mutate(variable = case_when(
    variable == "b_sd_trip.M" ~ "Male trip variability",
    variable == "b_sd_trip.F" ~ "Female trip variability",
    variable == "b_trip_days.M" ~ "Male trip duration (median)",
    variable == "b_trip_days.F" ~ "Female trip duration (median)",
    variable == "b_debt.days" ~ "Pair debt (male - female)"))


##### Plot ----------------------------------------------------------------

posteriors.incub_bba %<>%
  mutate(variable = fct_relevel(variable, c(
    "Male trip variability",
    "Female trip variability",
    "Male trip duration (median)",
    "Female trip duration (median)",
    "Pair debt (male - female)")))

my_labels <- c("Trip variability \n(male)", "Trip variability \n(female)",
               "Trip duration \n(male; median)", "Trip duration \n(female, median)", 
               "Pair debt \n(male-female)")

posteriors_plot.incub_bba.horizontal <- 
  ggplot(posteriors.incub_bba, aes(x = value, y = variable)) +
  stat_halfeye(data = subset(posteriors.incub_bba, colony == "Bird Island"), 
               side = "bottom", 
               aes(fill = ifelse(after_stat(x) < 0, "bi_low", "bi_high")),
               height = 0.8) +
  stat_halfeye(data = subset(posteriors.incub_bba, colony == "Kerguelen"), 
               side = "top", 
               aes(fill = ifelse(after_stat(x) < 0, "ker_low", "ker_high")),
               height = 0.8) +
  scale_fill_manual(values = c(
    "bi_low" = bi_col.low, "bi_high" = bi_col.high,
    "ker_low" = ker_col.low, "ker_high" = ker_col.high ) ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, col = "grey20") +
  xlim(-1.5, 1.5) +
  scale_y_discrete(labels = my_labels) +
  labs(y = "", x = "Posterior estimate",
       title = "Black-browed albatrosses") +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none")

#### WAAL  ---------------------------------------------------------------------
##### Estimates ---------------------------------------------------------------------
posterior_samples.incub_waal <- posterior_samples(incub_brms.waal.sex)

# Extract posterior samples by colony
predictors <- colnames(posterior_samples.incub_waal)[2:7]
predictors <- predictors[-2]
interaction_terms <- colnames(posterior_samples.incub_waal)[8:12]
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


## Rename and relevel factors
posteriors.incub_waal %<>%
  mutate(variable = case_when(
    variable == "b_sd_trip.M" ~ "Male trip variability",
    variable == "b_sd_trip.F" ~ "Female trip variability",
    variable == "b_trip_days.M" ~ "Male trip duration (median)",
    variable == "b_trip_days.F" ~ "Female trip duration (median)",
    variable == "b_debt.days" ~ "Pair debt (male - female)"))


##### Plot ---------------------------------------------------------------

posteriors.incub_waal %<>%
  mutate(variable = fct_relevel(variable, c(
    "Male trip variability",
    "Female trip variability",
    "Male trip duration (median)",
    "Female trip duration (median)",
    "Pair debt (male - female)")))

posteriors_plot.incub_waal.horizontal <- 
  ggplot(posteriors.incub_waal, aes(x = value, y = variable)) +
  stat_halfeye(data = subset(posteriors.incub_waal, colony == "Bird Island"), 
               side = "bottom", 
               aes(fill = ifelse(after_stat(x) < 0, "bi_low", "bi_high")),
               height = 0.7) +
  stat_halfeye(data = subset(posteriors.incub_waal, colony == "Crozet"), 
               side = "top", 
               aes(fill = ifelse(after_stat(x) < 0, "cro_low", "cro_high")),
               height = 0.7) +
  scale_fill_manual(values = c(
    "bi_low" = bi_col.low, "bi_high" = bi_col.high,
    "cro_low" = cro_col.low, "cro_high" = cro_col.high ) ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, col = "grey20") +
  labs(y = "", x = "Posterior estimate",
       title = "Wandering albatrosses") +
 # xlim(-1.3, 1) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none",
        axis.text.y = element_blank())



#### COMBINED posterior plots -----------------------------------------------------------
png(file = "Figures/rs_failure/posterior_estimates_incub_horizontal.bySex.png", 
    width = 15, height = 9, units = "in", res = 600)
ggarrange(posteriors_plot.incub_bba.horizontal + theme(plot.margin = unit(c(1,0.05,1,4), "cm")), 
          posteriors_plot.incub_waal.horizontal + theme(plot.margin = unit(c(1,1,1,0.75), "cm")), 
          ncol = 2,
          widths = c(1, 0.9))
dev.off()

# Laptop
# png(file = "Figures/rs_failure/posterior_estimates_incub_horizontal2.png", 
#     width = 15, height = 9, units = "in", res = 100)
# ggarrange(posteriors_plot.incub_bba.horizontal, posteriors_plot.incub_waal.horizontal,
#           ncol = 2, widths = c(1, 0.85))
# dev.off()

# ............................................................ ####
# INTERPRETATION --------------------------------------------------------------

load("Data_inputs/all_pair_behaviour_incub_withSex.RData")

load("Data_outputs/bba_incub_brms_model_withSex.RData")
load("Data_outputs/waal_incub_brms_model_withSex.RData")

# Datasets
behaviour.bba_incub <- pair_behaviour.incub_withSex %>% filter(species == "BBA") 
behaviour.waal_incub <- pair_behaviour.incub_withSex %>% filter(species == "WAAL")

above_zero.func <- function (variable) { mean(variable > 0) }
below_zero.func <- function (variable) { mean(variable < 0) }

### BBAL ---------------------------------------------------------------------
summary(incub_brms.bba.sex)
posterior_samples.incub_bba <- posterior_samples(incub_brms.bba.sex)

##### Overall summary table ------------------------------

bbal_incub.fixef <- round(summary(incub_brms.bba.sex)$fixed, digits = 2)
bbal_incub.fixef[,c(6,7)] <- round(bbal_incub.fixef[,c(6,7)], digits = 0)

# Get odds ratios
exp_coef <- round(exp(bbal_incub.fixef[,c(1, 3, 4)]), digits = 2)
colnames(exp_coef) <- c("Estimate_exp", "LCI_exp", "UCI_exp")
rownames(exp_coef) <- NULL

# Get prop above/below zero
below_zero <- data.frame(apply(posterior_samples.incub_bba[,1:12], 2, below_zero.func))
below_zero <- below_zero[,1]

above_zero <- data.frame(apply(posterior_samples.incub_bba[,1:12], 2, above_zero.func))
above_zero <- above_zero[,1]

## Bind everything together
bbal_incub.fixef <- cbind(bbal_incub.fixef, exp_coef)
bbal_incub.fixef$below_zero <- round(below_zero, digits = 2)
bbal_incub.fixef$above_zero <- round(above_zero, digits = 2)

## Get Kerguelen estimates separately
bbal_incub.fixef %<>% 
  mutate(Est.extra = NA, LCI.extra = NA, UCI.extra = NA, below.extra = NA, above.extra = NA)

### Process each variable for Kerguelen estimates
ker_trip_F <- process_interaction_estimates("b_trip_days.F", bbal_incub.fixef, posterior_samples.incub_bba, "b_trip_days.F:colonykerguelen")
ker_trip_M <- process_interaction_estimates("b_trip_days.M", bbal_incub.fixef, posterior_samples.incub_bba, "b_colonykerguelen:trip_days.M")
ker_sd_trip_F <- process_interaction_estimates("b_sd_trip.F", bbal_incub.fixef, posterior_samples.incub_bba, "b_colonykerguelen:sd_trip.F")
ker_sd_trip_M <- process_interaction_estimates("b_sd_trip.M", bbal_incub.fixef, posterior_samples.incub_bba, "b_colonykerguelen:sd_trip.M")
ker_debt <- process_interaction_estimates("b_debt.days", bbal_incub.fixef, posterior_samples.incub_bba, "b_colonykerguelen:debt.days")

# Combine the results with the original fixef data
ker_fixef <- rbind(ker_trip_F, ker_trip_M, ker_sd_trip_F, ker_sd_trip_M, ker_debt)
colnames(ker_fixef) <- colnames(bbal_incub.fixef)

bbal_incub.fixef <- rbind(bbal_incub.fixef, ker_fixef)
bbal_incub.fixef


### WAAL ---------------------------------------------------------------------
summary(incub_brms.waal.sex)
posterior_samples.incub_waal <- posterior_samples(incub_brms.waal.sex)

### Overall summary table ------------------------------

waal_incub.fixef <- round(summary(incub_brms.waal.sex)$fixed, digits = 2)
waal_incub.fixef[,c(6,7)] <- round(waal_incub.fixef[,c(6,7)], digits = 0)

# Get odds ratios
exp_coef <- round(exp(waal_incub.fixef[,c(1, 3, 4)]), digits = 2)
colnames(exp_coef) <- c("Estimate_exp", "LCI_exp", "UCI_exp")
rownames(exp_coef) <- NULL

# Get prop above/below zero
below_zero <- data.frame(apply(posterior_samples.incub_waal[,1:12], 2, below_zero.func))
below_zero <- below_zero[,1]

above_zero <- data.frame(apply(posterior_samples.incub_waal[,1:12], 2, above_zero.func))
above_zero <- above_zero[,1]

## Bind everything together
waal_incub.fixef <- cbind(waal_incub.fixef, exp_coef)
waal_incub.fixef$below_zero <- round(below_zero, digits = 2)
waal_incub.fixef$above_zero <- round(above_zero, digits = 2)

## Get Crozet estimates separately
waal_incub.fixef %<>% 
  mutate(Est.extra = NA, LCI.extra = NA, UCI.extra = NA, below.extra = NA, above.extra = NA)

### Process each variable for crozet estimates
cro_trip_F <- process_interaction_estimates("b_trip_days.F", waal_incub.fixef, posterior_samples.incub_waal, "b_trip_days.F:colonycrozet")
cro_trip_M <- process_interaction_estimates("b_trip_days.M", waal_incub.fixef, posterior_samples.incub_waal, "b_colonycrozet:trip_days.M")
cro_sd_trip_F <- process_interaction_estimates("b_sd_trip.F", waal_incub.fixef, posterior_samples.incub_waal, "b_colonycrozet:sd_trip.F")
cro_sd_trip_M <- process_interaction_estimates("b_sd_trip.M", waal_incub.fixef, posterior_samples.incub_waal, "b_colonycrozet:sd_trip.M")
cro_debt <- process_interaction_estimates("b_debt.days", waal_incub.fixef, posterior_samples.incub_waal, "b_colonycrozet:debt.days")

# Combine the results with the original fixef data
cro_fixef <- rbind(cro_trip_F, cro_trip_M, cro_sd_trip_F, cro_sd_trip_M, cro_debt)
colnames(cro_fixef) <- colnames(waal_incub.fixef)

waal_incub.fixef <- rbind(waal_incub.fixef, cro_fixef)
waal_incub.fixef

# ______________________________ ####
# ~ * SEX ANALYSIS - BROODING * ~ ###########################################################

# ............................................................ ####
# ANALYSIS =====================================================================

load("Data_inputs/all_pair_behaviour_brooding_withSex.RData")

#### Specify priors ---------------------------------------------------------------

# BBA 
priors.bba <- c(set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.F"),
                set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.M"),
                set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.M"),
                set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.F"),
                set_prior("normal(-0.12, 0.5)", class = "b", coef = "debt.days"))

# WAAL
priors.waal <- c(set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.F"),
                 set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.M"),
                 set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.M"),
                 set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.F"),
                 set_prior("normal(-0.12, 0.5)", class = "b", coef = "debt.days"))

#### Specify models ----------------------------------------------------------------

# Datasets
behaviour.bba_brooding <- pair_behaviour.brooding_withSex %>% filter(species == "BBA") 
behaviour.waal_brooding <- pair_behaviour.brooding_withSex %>% filter(species == "WAAL")

# Model structure
bf.brooding <- brms::bf(breeding_outcome.bin ~
                       trip_days.F * colony +
                       trip_days.M * colony +
                       sd_trip.F * colony +
                       sd_trip.M * colony + 
                       debt.days * colony +
                       (1|season), 
                     family = "bernoulli")

bf.brooding_bba <- brms::bf(breeding_outcome.bin ~
                          trip_days.F +
                          trip_days.M +
                          sd_trip.F +
                          sd_trip.M + 
                          debt.days +
                            colony +
                          (1|season), 
                        family = "bernoulli")

## Fit models
brooding_brms.bba.sex <- brm(bf.brooding_bba,
                          data = behaviour.bba_brooding, 
                          cores = 4, chains = 4, 
                          iter = 20000, warmup = 10000, thin = 10,
                          control = list(adapt_delta = 0.9999, max_treedepth = 14),
                          prior = priors.bba)

summary(brooding_brms.bba.sex)
pp_check(brooding_brms.bba.sex, ndraw = 50)
save(brooding_brms.bba.sex, file = "Data_outputs/bba_brooding_brms_model_withSex.RData")


brooding_brms.waal.sex <- brm(bf.brooding,
                           data = behaviour.waal_brooding, 
                           cores = 4, chains = 4, 
                           iter = 20000, warmup = 10000, thin = 10,
                           control = list(adapt_delta = 0.9999, max_treedepth = 14),
                           prior = priors.waal)

summary(brooding_brms.waal.sex)
pp_check(brooding_brms.waal.sex, ndraw = 50)
save(brooding_brms.waal.sex, file = "Data_outputs/waal_brooding_brms_model_withSex.RData")

# ............................................................ ####
# VISUALISE --------------------------------------------------------------------

load("Data_inputs/all_pair_behaviour_brooding_withSex.RData")

load("Data_outputs/bba_brooding_brms_model_withSex.RData")
load("Data_outputs/waal_brooding_brms_model_withSex.RData")

# Datasets
behaviour.bba_brooding <- pair_behaviour.brooding_withSex %>% filter(species == "BBA") 
behaviour.waal_brooding <- pair_behaviour.brooding_withSex %>% filter(species == "WAAL")

### Posterior estimates --------------------------------------------------------

#### BBAL  ---------------------------------------------------------------------
##### Estimates ----------------------------------------------------------------
posterior_samples.brooding_bba <- posterior_samples(brooding_brms.bba.sex)

# Extract posterior samples by colony
predictors <- colnames(posterior_samples.brooding_bba)[2:6]
posteriors.brooding_bba <- data.frame()

for (i in 1:length(predictors)) {
  
  # Extract the main effect (Bird Island)
  bird_est <- posterior_samples.brooding_bba[[predictors[i]]]
  
  # Extract the interaction effect and add it to the main effect for Kerguelen
  temp_df <- data.frame(
    variable = predictors[i],
    value = bird_est,
    colony = rep("Bird Island", length(bird_est))
  )
  
  posteriors.brooding_bba <- rbind(posteriors.brooding_bba, temp_df)
  
}


## Rename and relevel factors
posteriors.brooding_bba %<>%
  mutate(variable = case_when(
    variable == "b_sd_trip.M" ~ "Male trip variability",
    variable == "b_sd_trip.F" ~ "Female trip variability",
    variable == "b_trip_days.M" ~ "Male trip duration (median)",
    variable == "b_trip_days.F" ~ "Female trip duration (median)",
    variable == "b_debt.days" ~ "Pair debt (male - female)"))


##### Plot ----------------------------------------------------------------

posteriors.brooding_bba %<>%
  mutate(variable = fct_relevel(variable, c(
    "Male trip variability",
    "Female trip variability",
    "Male trip duration (median)",
    "Female trip duration (median)",
    "Pair debt (male - female)")))

my_labels <- c("Trip variability \n(male)", "Trip variability \n(female)",
               "Trip duration \n(male; median)", "Trip duration \n(female, median)", 
               "Pair debt \n(male-female)")

posteriors_plot.brooding_bba.horizontal <- 
  ggplot(posteriors.brooding_bba, aes(x = value, y = variable)) +
  stat_halfeye(data = posteriors.brooding_bba, 
               side = "top", 
               aes(fill = ifelse(after_stat(x) < 0, "low", "high")),
               height = 0.8) +
  scale_fill_manual(values = c(
    "low" = brooding_col.low, "high" = brooding_col.high ) ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, col = "grey20") +
  #xlim(-1.2, 1.2) +
  scale_y_discrete(labels = my_labels) +
  labs(y = "", x = "Posterior estimate",
       title = "Black-browed albatrosses") +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none")

#### WAAL  ---------------------------------------------------------------------
##### Estimates ---------------------------------------------------------------------
posterior_samples.brooding_waal <- posterior_samples(brooding_brms.waal.sex)

# Extract posterior samples by colony
predictors <- colnames(posterior_samples.brooding_waal)[2:7]
predictors <- predictors[-2]
interaction_terms <- colnames(posterior_samples.brooding_waal)[8:12]
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


## Rename and relevel factors
posteriors.brooding_waal %<>%
  mutate(variable = case_when(
    variable == "b_sd_trip.M" ~ "Male trip variability",
    variable == "b_sd_trip.F" ~ "Female trip variability",
    variable == "b_trip_days.M" ~ "Male trip duration (median)",
    variable == "b_trip_days.F" ~ "Female trip duration (median)",
    variable == "b_debt.days" ~ "Pair debt (male - female)"))


##### Plot ---------------------------------------------------------------

posteriors.brooding_waal %<>%
  mutate(variable = fct_relevel(variable, c(
    "Male trip variability",
    "Female trip variability",
    "Male trip duration (median)",
    "Female trip duration (median)",
    "Pair debt (male - female)")))

posteriors_plot.brooding_waal.horizontal <- 
  ggplot(posteriors.brooding_waal, aes(x = value, y = variable)) +
  stat_halfeye(data = subset(posteriors.brooding_waal, colony == "Bird Island"), 
               side = "bottom", 
               aes(fill = ifelse(after_stat(x) < 0, "bi_low", "bi_high")),
               height = 0.7) +
  stat_halfeye(data = subset(posteriors.brooding_waal, colony == "Crozet"), 
               side = "top", 
               aes(fill = ifelse(after_stat(x) < 0, "cro_low", "cro_high")),
               height = 0.7) +
  scale_fill_manual(values = c(
    "bi_low" = bi_col.low, "bi_high" = bi_col.high,
    "cro_low" = cro_col.low, "cro_high" = cro_col.high ) ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, col = "grey20") +
  labs(y = "", x = "Posterior estimate",
       title = "Wandering albatrosses") +
  # xlim(-1.3, 1) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none",
        axis.text.y = element_blank())



#### COMBINED posterior plots -----------------------------------------------------------
png(file = "Figures/rs_failure/posterior_estimates_brooding_horizontal.bySex.png", 
    width = 15, height = 9, units = "in", res = 600)
ggarrange(posteriors_plot.brooding_bba.horizontal + theme(plot.margin = unit(c(1,0.05,1,4), "cm")), 
          posteriors_plot.brooding_waal.horizontal + theme(plot.margin = unit(c(1,1,1,0.75), "cm")), 
          ncol = 2,
          widths = c(1, 0.9))
dev.off()

# Laptop
# png(file = "Figures/rs_failure/posterior_estimates_brooding_horizontal2.png", 
#     width = 15, height = 9, units = "in", res = 100)
# ggarrange(posteriors_plot.brooding_bba.horizontal, posteriors_plot.brooding_waal.horizontal,
#           ncol = 2, widths = c(1, 0.85))
# dev.off()

# ............................................................ ####
# INTERPRETATION --------------------------------------------------------------

load("Data_inputs/all_pair_behaviour_brooding_withSex.RData")

load("Data_outputs/bba_brooding_brms_model_withSex.RData")
load("Data_outputs/waal_brooding_brms_model_withSex.RData")

# Datasets
behaviour.bba_brooding <- pair_behaviour.brooding %>% filter(species == "BBA") 
behaviour.waal_brooding <- pair_behaviour.brooding %>% filter(species == "WAAL")

above_zero.func <- function (variable) { mean(variable > 0) }
below_zero.func <- function (variable) { mean(variable < 0) }

### BBAL ---------------------------------------------------------------------
summary(brooding_brms.bba.sex)
posterior_samples.brooding_bba <- posterior_samples(brooding_brms.bba.sex)

##### Overall summary table ------------------------------

bbal_brooding.fixef <- round(summary(brooding_brms.bba.sex)$fixed, digits = 2)
bbal_brooding.fixef[,c(6,7)] <- round(bbal_brooding.fixef[,c(6,7)], digits = 0)

# Get odds ratios
exp_coef <- round(exp(bbal_brooding.fixef[,c(1, 3, 4)]), digits = 2)
colnames(exp_coef) <- c("Estimate_exp", "LCI_exp", "UCI_exp")
rownames(exp_coef) <- NULL

# Get prop above/below zero
below_zero <- data.frame(apply(posterior_samples.brooding_bba[,1:7], 2, below_zero.func))
below_zero <- below_zero[,1]

above_zero <- data.frame(apply(posterior_samples.brooding_bba[,1:7], 2, above_zero.func))
above_zero <- above_zero[,1]

## Bind everything together
bbal_brooding.fixef <- cbind(bbal_brooding.fixef, exp_coef)
bbal_brooding.fixef$below_zero <- round(below_zero, digits = 2)
bbal_brooding.fixef$above_zero <- round(above_zero, digits = 2)

bbal_brooding.fixef


### WAAL ---------------------------------------------------------------------
summary(brooding_brms.waal.sex)
posterior_samples.brooding_waal <- posterior_samples(brooding_brms.waal.sex)

### Overall summary table ------------------------------

waal_brooding.fixef <- round(summary(brooding_brms.waal.sex)$fixed, digits = 2)
waal_brooding.fixef[,c(6,7)] <- round(waal_brooding.fixef[,c(6,7)], digits = 0)

# Get odds ratios
exp_coef <- round(exp(waal_brooding.fixef[,c(1, 3, 4)]), digits = 2)
colnames(exp_coef) <- c("Estimate_exp", "LCI_exp", "UCI_exp")
rownames(exp_coef) <- NULL

# Get prop above/below zero
below_zero <- data.frame(apply(posterior_samples.brooding_waal[,1:12], 2, below_zero.func))
below_zero <- below_zero[,1]

above_zero <- data.frame(apply(posterior_samples.brooding_waal[,1:12], 2, above_zero.func))
above_zero <- above_zero[,1]

## Bind everything together
waal_brooding.fixef <- cbind(waal_brooding.fixef, exp_coef)
waal_brooding.fixef$below_zero <- round(below_zero, digits = 2)
waal_brooding.fixef$above_zero <- round(above_zero, digits = 2)

## Get Crozet estimates separately
waal_brooding.fixef %<>% 
  mutate(Est.extra = NA, LCI.extra = NA, UCI.extra = NA, below.extra = NA, above.extra = NA)

### Process each variable for crozet estimates
cro_trip_F <- process_interaction_estimates("b_trip_days.F", waal_brooding.fixef, posterior_samples.brooding_waal, "b_trip_days.F:colonycrozet")
cro_trip_M <- process_interaction_estimates("b_trip_days.M", waal_brooding.fixef, posterior_samples.brooding_waal, "b_colonycrozet:trip_days.M")
cro_sd_trip_F <- process_interaction_estimates("b_sd_trip.F", waal_brooding.fixef, posterior_samples.brooding_waal, "b_colonycrozet:sd_trip.F")
cro_sd_trip_M <- process_interaction_estimates("b_sd_trip.M", waal_brooding.fixef, posterior_samples.brooding_waal, "b_colonycrozet:sd_trip.M")
cro_debt <- process_interaction_estimates("b_debt.days", waal_brooding.fixef, posterior_samples.brooding_waal, "b_colonycrozet:debt.days")

# Combine the results with the original fixef data
cro_fixef <- rbind(cro_trip_F, cro_trip_M, cro_sd_trip_F, cro_sd_trip_M, cro_debt)
colnames(cro_fixef) <- colnames(waal_brooding.fixef)

waal_brooding.fixef <- rbind(waal_brooding.fixef, cro_fixef)
waal_brooding.fixef

# +++++++++++++++++++++++++++++ ####

# * FIGURE SX * COMBINED POSTERIOR PLOTS ========================================

# Add albatross silhouettes
posteriors_plot.incub_bba.horizontal2 <- ggdraw() +
  draw_plot(posteriors_plot.incub_bba.horizontal + labs(y = "Incubation") +
              theme(axis.title.y = element_text(margin = margin(r = 90)),
                    axis.title.x = element_blank(),
                    text = element_text(size = 16, family = "Calibri"))) +
  draw_image(file.path("Figures/rs_failure/bba_standing_silhouette.png"),
             scale = 0.15, x = 0.39, y = 0.28) 


posteriors_plot.incub_waal.horizontal2 <- ggdraw() +
  draw_plot(posteriors_plot.incub_waal.horizontal +
              theme(axis.title.x = element_blank(),
                    text = element_text(size = 16, family = "Calibri"))) +
  draw_image(file.path("Figures/rs_failure/waal_standing_silhouette.png"),
             scale = 0.2, x = 0.38, y = 0.32)


png(file = "Figures/rs_failure/FIGURESX.png", width = 10, height = 10, units = "in", res = 100)
ggarrange(posteriors_plot.incub_bba.horizontal2,
          posteriors_plot.incub_waal.horizontal2,
          posteriors_plot.brooding_bba.horizontal + 
            labs(y = "Brooding") + 
            theme(axis.title.y = element_text(margin = margin(r = 90)),
                  plot.title = element_blank()),
          posteriors_plot.brooding_waal.horizontal + 
            theme(plot.title = element_blank()),
          ncol = 2,
          nrow = 2,
          widths = c(1, 0.8))
dev.off()

# +++++++++++++++++++++++++++++ ####

# ______________________________ ####
# * VISUALISATION * ----------------------------------------------------------------


## Age of failure for each stage -----------------------------------------------
fail_dates <- bas_trips %>% 
  distinct() %>%
  mutate(fail_age = as.numeric(fail_date - lay_date),
         chick_age = as.numeric(fail_date - hatch_date))

## Set species names and colours
species_names <- c( bba = "Black-browed albatross", waal = "Wandering albatross")

## Set line for typical hatch date
vline_data <- data.frame(
  species = c("bba", "waal"),  
  vline_position = c(68, 78) )

## Plot it out
failure_ages.plot <- ggplot(data = fail_dates, aes(x = fail_age)) + 
  geom_histogram(aes(fill = rs), alpha = 0.8) +
  facet_grid(species ~ ., labeller = as_labeller(species_names)) +
  labs(x = "Days since laying", y = "Frequency") + 
  scale_fill_manual(values = c(fail_unknown, fail_chick, fail_egg),
                    name = "",
                    labels = c("Unknown failure", "Chick failure", "Egg failure")) +
  geom_vline(data = vline_data, aes(xintercept = vline_position), 
             linetype = "dashed", color = "black", size = 1) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = c(0.89, 0.95),
        legend.background = element_rect(fill = alpha("white", 0)))

## * FIGURE SX : FAILURE AGES * ===============================================================
png(file = "Figures/rs_failure/FIGURESX_failure_ages.png", 
    width = 9, height = 7, units = "in", res = 600)
failure_ages.plot
dev.off()

# ............................................................ ####
## Trip duration relative to age -----------------------------------------------

species_names <- c( bba = "Black-browed albatross", waal = "Wandering albatross")

trip_duration_age.plot <-ggplot(data = bas_trips, aes(x = chick_age, y = duration.days, col = phase)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(y = "Trip duration (days)", x = "Chick age (days)") +
  facet_grid(species ~ ., labeller = as_labeller(species_names)) +
  scale_colour_manual(values = c(brooding_col, incub_col),
                    name = "",
                    labels = c("Brooding", "Incubation")) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = c(0.89, 0.95),
        legend.background = element_rect(fill = alpha("white", 0))) +
  guides(color=guide_legend(override.aes=list(fill=NA)))


## * FIGURE SX : TRIP DURATION BY AGE* ===============================================================
png(file = "Figures/rs_failure/FIGURESX_trip_duration_age.png", 
    width = 9, height = 7, units = "in", res = 600)
trip_duration_age.plot
dev.off()



# ______________________________ ####
# ADDITIONAL ANALYSES ----------------------------------------------------------

## Trip duration comparison  ---------------------------------------------------

### Process datasets -----------------------------------------------------------

spec_col.codes <- c("bba_birdis", "waal_birdis", "waal_crozet", "bba_kerguelen")
all_data.list <- list()

for (i in 1:length(spec_col.codes)) {

  print(paste0("Processing ", spec_col.codes[i]))
  
  load(paste0("Data_inputs/", spec_col.codes[i], "_gls_labelled.RData"))
  load(paste0("Data_inputs/", spec_col.codes[i], "_gps_labelled.RData"))
  
  devices <- c("gls", "gps")
  
  for (d in 1:length(devices)) { 
    
    load(paste0("Data_inputs/", spec_col.codes[i], "_", devices[d], "_labelled.RData"))
    labelled.data <- get(paste0(devices[d], "_labelled.df")) 
    labelled.data %<>% mutate(ringYr = paste0(ring, "_", season))
    
    my_birds <- unique(labelled.data$ringYr)
    
    device_level.list <- list()
    
    print(paste0("Processing ", devices[d]))
    
    pb <- txtProgressBar(min = 0, max = length(my_birds), style = 3)
    
    for (b in 1:length(my_birds)) {
    
      setTxtProgressBar(pb, b)
      
      subset.data <- subset(labelled.data, ringYr == my_birds[b])
      subset.data %<>% arrange(datetime)
      subset.data$tripID = rep(1:length(rle(subset.data$loc)$values), rle(subset.data$loc)$lengths)
  
      my_trips <- subset.data %>% 
        group_by(ring, season, tripID) %>%
        mutate(start = datetime[1],
             end = datetime[n()],
             duration.mins = as.numeric(difftime(end, start, unit = "mins"))) %>%
        arrange(start) %>%
        filter(!is.na(datetime) & !is.na(end))
  
       ## Get rid of very short trips/col visits
      if (all(my_trips$duration.mins < 480)) next
       
      max_tries <- 5
      tries <- 0
       
      while(any(my_trips$duration.mins <= 480) & tries <= max_tries) {
         
         tries = tries + 1
         
         my_trips$loc <- ifelse(my_trips$duration.mins <= 480, 
                                ifelse(my_trips$loc == "trip", "col", "trip"),
                                my_trips$loc)
         my_trips$tripID <- rep(1:length(rle(my_trips$loc)$values), rle(my_trips$loc)$length)
         
         my_trips %<>%
           group_by(ring, season, tripID) %>%
           mutate(start = datetime[1],
                  end = datetime[n()],
                  duration.mins = as.numeric(difftime(end, start, unit = "mins"))) %>%
           arrange(start)
         
      }
       
      # Summarise the trips and get species/colony
      which_species <- ifelse(grepl("bba", spec_col.codes[i]), "bba", "waal")
      which_col <- tolower(gsub(paste0(which_species, "_"), "", spec_col.codes[i]))
      
      my_trips %<>%
         group_by(ring, season, tripID) %>%
         dplyr::summarise(start = datetime[1],
                          end = datetime[n()],
                          loc = loc[1]) %>%
         mutate(duration.mins = as.numeric(difftime(end, start, unit = "mins")),
                duration.days = duration.mins/60/24,
                device = devices[d],
                species = which_species,
                col = which_col) %>%
         arrange(start) %>%
         filter(duration.mins >= 480)
      
      all_data.list[[length(all_data.list) + 1]] <- my_trips
  
    }
    
  }

}

close(pb)

all_data.list <- lapply(all_data.list, function(df) {
  df$start <- as.POSIXct(df$start, tz = "UTC")  
  df$end <- as.POSIXct(df$end, tz = "UTC")      
  return(df)
})

all_data <- do.call("rbind", all_data.list)
all_data %<>% filter(duration.mins > 360 & duration.days < 31)

save(all_data, file = "Data_inputs/gls_gps_comparison_data.RData")

#### Cut data to breeding ------------------------------------------------------

load("Data_inputs/gls_gps_comparison_data.RData")

breeding_table <- data.frame(species = c("bba", "bba", "waal", "waal"),
                             col = c("birdis", "kerguelen", "birdis", "crozet"),
                             stt_breed = c("10-27", "10-09", "12-21", "12-25"),
                             end_breed = c("01-25", "01-04", "04-19", "04-23"))

# Merge into all_data
all_data <- merge(all_data, breeding_table, by = c("species", "col"), all.x = T)
all_data %<>% 
  mutate(stt_breed = as.Date(paste0(season, "-", stt_breed), format = "%Y-%m-%d"),
         end_breed = as.Date(paste0(season, "-", end_breed), format = "%Y-%m-%d")) %>%
  mutate(stt_month = as.numeric(format(as.Date(stt_breed), "%m")),
         stt_breed = ifelse(stt_month > 9,
                    as.character(as.Date(stt_breed) - 365),
                    as.character(stt_breed)))

all_data %<>% filter(as.Date(start) >= stt_breed & as.Date(start) <= end_breed) %>%
  select(-c(stt_breed, end_breed,stt_month)) %>%
  mutate(device = toupper(device))

### Look at comparisons -----------------------------------------------------

# Split by species
tripDur.bba <- subset(all_data, species == "bba")
tripDur.waal <- subset(all_data, species == "waal")

# Build models
## BBAL
dur_lmm.bba <- lmer(duration.days ~ device * col + (1|ring) + (1|season), data = tripDur.bba)
summary(dur_lmm.bba)

## WAAL
dur_lmm.waal <- lmer(duration.days ~ device * col + (1|ring) + (1|season), data = tripDur.waal)
summary(dur_lmm.waal)


# Make plots
## BBAL
summary.bba <- summarySE(tripDur.bba, measurevar = "duration.days", groupvars = c("device", "col"))

tripDur_bba.plot <- ggplot() + 
  geom_point(data = bind_rows(
    sample_n(filter(tripDur.bba, device == "GLS"), 1000),  
    filter(tripDur.bba, device == "GPS")),
    aes(x = device, y = duration.days, col = col),
    position = position_jitter(width = 0.2),
    alpha = 0.2) +
  geom_errorbar(data = summary.bba, aes(x=device, y=duration.days, group = col,
                                        ymin=duration.days-se, ymax=duration.days+se,
                                        col = col), width=0.5) +
  geom_point(data = summary.bba, 
             aes(x=device, y=duration.days, group = col, col = col),
             size = 1.5) +
  scale_colour_manual(values = c(bi_col, ker_col), 
                      labels = c("Bird Island", "Kerguelen"),
                      name = "") +
  ylim(0, 20) +
  labs(y = "Trip duration (days)", x = "Device", 
       title = "Black-browed albatrosses",
       caption = "") +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = c(0.75, 0.97),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))


## WAAL
summary.waal <- summarySE(tripDur.waal, measurevar = "duration.days", groupvars = c("device", "col"))

tripDur_waal.plot <- ggplot() + 
  geom_point(data = bind_rows(
    sample_n(filter(tripDur.waal, device == "GLS"), 1000),  
    sample_n(filter(tripDur.waal, device == "GPS"), 1000)),
    aes(x = device, y = duration.days, col = col),
    position = position_jitter(width = 0.2),
    alpha = 0.2) +
  geom_errorbar(data = summary.waal, aes(x=device, y=duration.days, group = col,
                                        ymin=duration.days-se, ymax=duration.days+se,
                                        col = col), width=0.5) +
  geom_point(data = summary.waal, 
             aes(x=device, y=duration.days, group = col, col = col),
             size = 1.5) +
  scale_colour_manual(values = c(bi_col, cro_col), 
                      labels = c("Bird Island", "Crozet"),
                      name = "") +
  labs(x = "Device", 
      title = "Wandering albatrosses") +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        plot.caption = element_text(size = 12, face = "italic"),
        axis.title.y = element_blank(),
        legend.position = c(0.75, 0.97),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.key = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5)) 



## * FIGURE SX : TRIP DURATION BY DEVICE * ===============================================================

png(file = "Figures/rs_failure/FIGURES1_tripDur_deviceEffects.png", 
    width = 10, height = 6, units = "in", res = 300)
ggarrange(tripDur_bba.plot +
            theme(axis.title.y = element_text(margin = margin(r = 15)),
              plot.margin = unit(c(0, 1, 0, 0), "cm")), 
          tripDur_waal.plot +
            theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm")),
          ncol = 2,
          nrow = 1,
          widths = c(1, 0.875))
dev.off()
