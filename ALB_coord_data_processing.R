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

species <- c("baBI", "baKer", "waBI", "waCro")

# ______________________________ ####
# LOAD DATA ==========================================================

# Remove implausible trips and birds where breeding outcome is unknown

## BBAL - BIRD ISLAND (baBI) ##
load("Data_inputs/bbal_birdis_breeding_trips.RData")
baBI_trips %<>% distinct() %>% 
  filter(duration.mins > 360 & duration.days < 31 & !is.na(rs))

## BBAL - KERGUELEN (baKer) ##
load("Data_inputs/bbal_ker_breeding_trips.RData")
baKer_trips %<>% distinct() %>% 
  filter(duration.mins > 360 & duration.days < 31 & !is.na(rs) & rs != "UNKNOWN")

## WAAL - BIRD ISLAND (waBI) ##
load("Data_inputs/waal_birdis_breeding_trips.RData")
waBI_trips %<>% distinct() %>% 
  filter(duration.mins > 360 & duration.days < 31 & !is.na(rs))

## WAAL - CROZET (waCro) ##
load("Data_inputs/waal_crozet_breeding_trips.RData")
waCro_trips %<>% distinct() %>%
  filter(duration.mins > 360 & duration.days < 31 & !is.na(rs) & rs != "UNKNOWN" & rs != "NON-BREEDER")


# How many breeding dates were estimated? ---------------------------------

breeding_dates.df <- rbind(baBI_trips %>% select(c(pairID, lay_code, hatch_code)), 
                   baKer_trips %>% select(c(pairID, lay_code, hatch_code)),
                   waBI_trips %>% select(c(pairID, lay_code, hatch_code)), 
                   waCro_trips %>% select(c(pairID, lay_code, hatch_code)))

(nrow(subset(breeding_dates.df, lay_code == "EST" & hatch_code == "EST"))/nrow(breeding_dates.df))*100


# Data for sex differences ---------------------------------

all_trips <- rbind(baBI_trips %>% mutate(colony = "birdisland", species = "BBAL", speCol = "baBI"), 
                   baKer_trips %>% mutate(colony = "kerguelen", species = "BBAL", speCol = "baKer"), 
                   waBI_trips %>% mutate(colony = "birdisland", species = "WAAL", speCol = "waBI"),
                   waCro_trips %>% mutate(colony = "crozet", species = "WAAL", speCol = "waCro"))

save(all_trips, file = "Data_outputs/all_trips.RData")



# Overall breeding success ------------------------------------------------

all_trips %>% 
  select(pairID, season, rs, species) %>%
  filter(rs != "NON-BREEDER") %>% 
  mutate(outcome = ifelse(rs == "SUCCESSFUL", 1, 0)) %>%
  group_by(species) %>%
  summarise(mean(outcome))


# ______________________________ ####
# DATA PROCESSING --------------------------------------------------------------

### Incubation -----------------------------------------------------------------

for (i in 1:length(species)) {
  
  my_species = species[i]
  species_trips <- get(paste0(my_species, "_trips"))
  species_trips %<>% distinct() %>% arrange(pairID, start)
  
  species_trips %<>%
    filter(phase == "incubation") %>%
    group_by(pairID, season) %>%
    filter(n() >= 4) %>%
    arrange(start) %>%
    mutate(prev_trip.mins = lag(duration.mins)) %>%
    group_by(ring) %>%
    mutate(fasting_debt.hours = (duration.mins - prev_trip.mins)/60,
           duration.hours = duration.mins/60) %>%
    select(-prev_trip.mins) %>%
    distinct()

  ## Trip duration & variability
  
  pair_behaviour <- species_trips %>%  
    group_by(pairID, season) %>%
    mutate(bird1 = ring[1],
           bird2 = partner[1]) %>%
    summarise(breeding_outcome = rs[1],
              colony = colony[1],
              species = species[1],
              # Trip variability
              sd_trip.1 = sd(duration.hours[ring == bird1], na.rm = T),
              sd_trip.2 = sd(duration.hours[ring == bird2], na.rm = T),
              diff_var_abs = abs(sd_trip.1 - sd_trip.2),
              # Median trip duration
              median_trip = median(na.omit(duration.hours)),
              median_trip.1 = median(na.omit(duration.hours[ring == bird1])),
              median_trip.2 = median(na.omit(duration.hours[ring == bird2])),
              median_trip.days = median(na.omit(duration.hours)),
              diff_trip_abs = abs(median_trip.1 - median_trip.2)) %>%
    # Binary breeding outcome
    mutate(breeding_outcome.bin = ifelse(breeding_outcome == "SUCCESSFUL", 1, 0)) %>%
    select(-c(median_trip.1, median_trip.2, sd_trip.1, sd_trip.2))

  pair_behaviour %<>% 
    mutate(median_trip.days = median_trip/24,
           diff_trip_abs.days = diff_trip_abs/24,
           diff_var_abs.days = diff_var_abs/24) %>%
    select(-c(median_trip, diff_trip_abs, diff_var_abs))
  
  ## Output the data
  assign(paste0(my_species, "_pair_behaviour"), pair_behaviour)
  
}

## Make a combined dataframe with dummy species/col variable
pair_behaviour.incub <- rbind(baBI_pair_behaviour %>% mutate(colony = "birdisland", species = "BBAL", speCol = "baBI"), 
                              baKer_pair_behaviour %>% mutate(colony = "kerguelen", species = "BBAL", speCol = "baKer"), 
                              waBI_pair_behaviour %>% mutate(colony = "birdisland", species = "WAAL", speCol = "waBI"),
                              waCro_pair_behaviour %>% mutate(colony = "crozet", species = "WAAL", speCol = "waCro"))

save(pair_behaviour.incub, file = "Data_inputs/all_pair_behaviour_incub.RData")



### Brooding ------------------------------------------------------------------

for (i in 1:length(species)) {
  
  my_species = species[i]
  species_trips <- get(paste0(my_species, "_trips"))
  species_trips %<>% distinct() %>% arrange(pairID, start)
  
  species_trips %<>%
    filter(phase == "brooding" & !is.na(sex)) %>%
    group_by(pairID, season) %>%
    filter(n() >= 4) %>%
    arrange(start) %>%
    mutate(prev_trip.mins = lag(duration.mins)) %>%
    group_by(ring) %>%
    mutate(fasting_debt.hours = (duration.mins - prev_trip.mins)/60,
           duration.hours = duration.mins/60) %>%
    select(-prev_trip.mins) %>%
    distinct()
  
  ## Trip duration & variability
  pair_behaviour <- species_trips %>%  
    group_by(pairID, season) %>%
    mutate(bird1 = ring[1],
           bird2 = partner[1]) %>%
    summarise(breeding_outcome = rs[1],
              colony = colony[1],
              species = species[1],
              # Trip variability
              sd_trip.1 = sd(duration.hours[ring == bird1], na.rm = T),
              sd_trip.2 = sd(duration.hours[ring == bird2], na.rm = T),
              diff_var_abs = abs(sd_trip.1 - sd_trip.2),
              # Median trip duration
              median_trip = median(na.omit(duration.hours)),
              median_trip.1 = median(na.omit(duration.hours[ring == bird1])),
              median_trip.2 = median(na.omit(duration.hours[ring == bird2])),
              median_trip.days = median(na.omit(duration.hours)),
              diff_trip_abs = abs(median_trip.1 - median_trip.2)) %>%
    # Binary breeding outcome
    mutate(breeding_outcome.bin = ifelse(breeding_outcome == "SUCCESSFUL", 1, 0)) %>%
    select(-c(median_trip.1, median_trip.2, sd_trip.1, sd_trip.2))
  
  pair_behaviour %<>% 
    mutate(
      median_trip.days = median_trip/24,
      diff_trip_abs.days = diff_trip_abs/24,
      diff_var_abs.days = diff_var_abs/24) %>%
    select(-c(median_trip, diff_trip_abs, diff_var_abs))
  
  ## Output the data
  assign(paste0(my_species, "_pair_behaviour"), pair_behaviour)
  
}

## Make a combined dataframe with dummy species/col variable
pair_behaviour.brooding <- rbind(baBI_pair_behaviour %>% mutate(colony = "birdisland", species = "BBAL", speCol = "baBI"), 
                                 baKer_pair_behaviour %>% mutate(colony = "kerguelen", species = "BBAL", speCol = "baKer"), 
                                 waBI_pair_behaviour %>% mutate(colony = "birdisland", species = "WAAL", speCol = "waBI"),
                                 waCro_pair_behaviour %>% mutate(colony = "crozet", species = "WAAL", speCol = "waCro"))

save(pair_behaviour.brooding, file = "Data_inputs/all_pair_behaviour_brooding.RData")

