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

## BBA - BIRD ISLAND (baBI) ##
load("Data_inputs/ANON_bba_birdis_breedingTrips.RData")
baBI_trips %<>% distinct() %>% 
  filter(duration.mins > 360 & duration.days < 31 & !is.na(rs))

## BBA - KERGUELEN (baKer) ##
load("Data_inputs/ANON_bba_kerguelen_breedingTrips.RData")
baKer_trips %<>% distinct() %>% 
  filter(duration.mins > 360 & duration.days < 31 & !is.na(rs) & rs != "UNKNOWN")

## WAAL - BIRD ISLAND (waBI) ##
load("Data_inputs/ANON_waal_birdis_breedingTrips.RData")
waBI_trips %<>% distinct() %>% 
  filter(duration.mins > 360 & duration.days < 31 & !is.na(rs))

## WAAL - CROZET (waCro) ##
load("Data_inputs/ANON_waal_crozet_breedingTrips.RData")
waCro_trips %<>% distinct() %>%
  filter(duration.mins > 360 & duration.days < 31 & !is.na(rs) & rs != "UNKNOWN" & rs != "NON-BREEDER")


# How many breeding dates were estimated? ---------------------------------

all_trips <- rbind(baBI_trips %>% select(c(pairID, lay_code, hatch_code)), 
                   baKer_trips %>% select(c(pairID, lay_code, hatch_code)),
                   waBI_trips %>% select(c(pairID, lay_code, hatch_code)), 
                   waCro_trips %>% select(c(pairID, lay_code, hatch_code)))

(nrow(subset(all_trips, lay_code == "EST" & hatch_code == "EST"))/nrow(all_trips))*100

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
  
  ## Fasting debt
  fasting_metrics <- species_trips %>%
    group_by(pairID, season) %>%
    mutate(bird1 = ring[1],
           bird2 = partner[1],
           n_trips = n(),
           trip_rank = row_number()) %>%
    filter(!(n_trips %% 2 != 0 & trip_rank == max(trip_rank))) %>%
    summarise(inc_time.1 = sum(duration.hours[ring == bird2], na.rm = TRUE),
              inc_time.2 = sum(duration.hours[ring == bird1], na.rm = TRUE),
              debt.days = abs(inc_time.1 - inc_time.2)) %>%
    select(-c(inc_time.1, inc_time.2))
  
  ## Trip duration & variability
  trip_metrics <- species_trips %>%
    group_by(pairID, season) %>%
    mutate(bird1 = ring[1],
           bird2 = partner[1]) %>%
    summarise(breeding_outcome = rs[1],
              colony = colony[1],
              species = species[1],
              # Trip variability
              sd_trip.1 = sd(duration.hours[ring == bird1], na.rm = T),
              sd_trip.2 = sd(duration.hours[ring == bird2], na.rm = T),
              diff_var.days = abs(sd_trip.1 - sd_trip.2),
              # Median trip duration
              median_trip.1 = median(na.omit(duration.hours[ring == bird1])),
              median_trip.2 = median(na.omit(duration.hours[ring == bird2])),
              median_trip.days = median(na.omit(duration.hours)),
              diff_trip.days = abs(median_trip.1 - median_trip.2)) %>%
    # Binary breeding outcome
    mutate(breeding_outcome.bin = ifelse(breeding_outcome == "SUCCESSFUL", 1, 0)) %>%
    select(-c(median_trip.1, median_trip.2, sd_trip.1, sd_trip.2))
  
  # Merge behaviour
  pair_behaviour <- merge(trip_metrics, fasting_metrics, by = c("pairID", "season"))
  pair_behaviour %<>% 
    mutate(debt.days = debt.days/24,
           median_trip.days = median_trip.days/24,
           diff_trip.days = diff_trip.days/24,
           diff_var.days = diff_var.days/24)
  
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
    filter(phase == "brooding") %>%
    group_by(pairID, season) %>%
    filter(n() >= 4) %>%
    arrange(start) %>%
    mutate(prev_trip.mins = lag(duration.mins)) %>%
    group_by(ring) %>%
    mutate(fasting_debt.hours = (duration.mins - prev_trip.mins)/60,
           duration.hours = duration.mins/60) %>%
    select(-prev_trip.mins) %>%
    distinct()
  
  ## Fasting debt
  fasting_metrics <- species_trips %>%
    group_by(pairID, season) %>%
    mutate(bird1 = ring[1],
           bird2 = partner[1],
           n_trips = n(),
           trip_rank = row_number()) %>%
    filter(!(n_trips %% 2 != 0 & trip_rank == max(trip_rank))) %>%
    summarise(inc_time.1 = sum(duration.hours[ring == bird2], na.rm = TRUE),
              inc_time.2 = sum(duration.hours[ring == bird1], na.rm = TRUE),
              debt.days = abs(inc_time.1 - inc_time.2)) %>%
    select(-c(inc_time.1, inc_time.2))
  
  ## Trip duration & variability
  trip_metrics <- species_trips %>%
    group_by(pairID, season) %>%
    mutate(bird1 = ring[1],
           bird2 = partner[1]) %>%
    summarise(breeding_outcome = rs[1],
              colony = colony[1],
              species = species[1],
              # Trip variability
              sd_trip.1 = sd(duration.hours[ring == bird1], na.rm = T),
              sd_trip.2 = sd(duration.hours[ring == bird2], na.rm = T),
              diff_var.days = abs(sd_trip.1 - sd_trip.2),
              # Median trip duration
              median_trip.1 = median(na.omit(duration.hours[ring == bird1])),
              median_trip.2 = median(na.omit(duration.hours[ring == bird2])),
              median_trip.days = median(na.omit(duration.hours)),
              diff_trip.days = abs(median_trip.1 - median_trip.2)) %>%
    # Binary breeding outcome
    mutate(breeding_outcome.bin = ifelse(breeding_outcome == "SUCCESSFUL", 1, 0)) %>%
    select(-c(median_trip.1, median_trip.2, sd_trip.1, sd_trip.2))
  
  # Merge behaviour
  pair_behaviour <- merge(trip_metrics, fasting_metrics, by = c("pairID", "season"))
  pair_behaviour %<>% 
    mutate(debt.days = debt.days/24,
           median_trip.days = median_trip.days/24,
           diff_trip.days = diff_trip.days/24,
           diff_var.days = diff_var.days/24)
  
  ## Output the data
  assign(paste0(my_species, "_pair_behaviour"), pair_behaviour)
  
}

## Make a combined dataframe with dummy species/col variable
pair_behaviour.brooding <- rbind(baBI_pair_behaviour %>% mutate(colony = "birdisland", species = "BBAL", speCol = "baBI"), 
                                 baKer_pair_behaviour %>% mutate(colony = "kerguelen", species = "BBAL", speCol = "baKer"), 
                                 waBI_pair_behaviour %>% mutate(colony = "birdisland", species = "WAAL", speCol = "waBI"),
                                 waCro_pair_behaviour %>% mutate(colony = "crozet", species = "WAAL", speCol = "waCro"))

save(pair_behaviour.brooding, file = "Data_inputs/all_pair_behaviour_brooding.RData")

