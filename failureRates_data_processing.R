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

## BBA - BIRD ISLAND (baBI) ##
load("Data_inputs/bba_birdis_breedingTrips.RData")
baBI_trips <- all_trips %>% distinct(); rm(all_trips)
baBI_trips %<>% mutate(tripID = NA)

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


# ______________________________ ####
# DATA PROCESSING --------------------------------------------------------------

### Incubation -----------------------------------------------------------------

for (i in 1:length(species)) {
  
  my_species = species[i]
  species_trips <- get(paste0(my_species, "_trips"))
  species_trips %<>% select(-tripID) %>% distinct()
  species_trips %<>% arrange(pairID, start)
  
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
  
  ## For birds where sex is known ##
  fasting_metrics.knownSex <- species_trips %>%
    filter(!is.na(sex)) %>%
    group_by(pairID, season) %>%
    mutate(bird1 = ring[1],
           bird2 = partner[1],
           n_trips = n(),
           trip_rank = row_number()) %>%
    filter(!(n_trips %% 2 != 0 & trip_rank == max(trip_rank))) %>%
    summarise(debt.F = sum(fasting_debt.hours[sex == "F"], na.rm = T),
              debt.M = sum(fasting_debt.hours[sex == "M"], na.rm = T),
              debt.days = (debt.M - debt.F)/24) %>%
    select(-c(debt.F, debt.M))
  
  pair_behaviour.known_sex <- species_trips %>%
    filter(!is.na(sex)) %>%
    group_by(pairID, season) %>%
    summarise(trip_days.F = median(duration.days[sex == "F"], na.rm = T),
              trip_days.M = median(duration.days[sex == "M"], na.rm = T),
              sd_trip.F = sd(duration.days[sex == "F"], na.rm = T),
              sd_trip.M = sd(duration.days[sex == "M"], na.rm = T),
              male_female_diff.var = sd_trip.M - sd_trip.F, 
              breeding_outcome = rs[1],
              colony = colony[1],
              species = species[1]) %>%
    # Binary breeding outcome
    mutate(breeding_outcome.bin = ifelse(breeding_outcome == "SUCCESSFUL", 1, 0))
  
  pair_behaviour.known_sex <- merge(pair_behaviour.known_sex, fasting_metrics.knownSex,
                                    by = c("pairID", "season"), all.x = TRUE)
  
  ## Output the data
  assign(paste0(my_species, "_pair_behaviour"), pair_behaviour)
  assign(paste0(my_species, "_pair_behaviour.known_sex"), pair_behaviour.known_sex)
  
}

## Make a combined dataframe with dummy species/col variable
pair_behaviour.incub <- rbind(baBI_pair_behaviour %>% mutate(colony = "birdisland", species = "BBA", speCol = "baBI"), 
                              baKer_pair_behaviour %>% mutate(colony = "kerguelen", species = "BBA", speCol = "baKer"), 
                              waBI_pair_behaviour %>% mutate(colony = "birdisland", species = "WAAL", speCol = "waBI"),
                              waCro_pair_behaviour %>% mutate(colony = "crozet", species = "WAAL", speCol = "waCro"))

pair_behaviour.incub_withSex <- rbind(baBI_pair_behaviour.known_sex %>% mutate(colony = "birdisland", species = "BBA", speCol = "baBI"), 
                                      baKer_pair_behaviour.known_sex %>% mutate(colony = "kerguelen", species = "BBA", speCol = "baKer"), 
                                      waBI_pair_behaviour.known_sex %>% mutate(colony = "birdisland", species = "WAAL", speCol = "waBI"),
                                      waCro_pair_behaviour.known_sex %>% mutate(colony = "crozet", species = "WAAL", speCol = "waCro"))


save(pair_behaviour.incub, file = "Data_inputs/all_pair_behaviour_incub.RData")
save(pair_behaviour.incub_withSex, file = "Data_inputs/all_pair_behaviour_incub_withSex.RData")




### Brooding ------------------------------------------------------------------

## lots of NA values for differences at Kerguelen?


for (i in 1:length(species)) {
  
  my_species = species[i]
  species_trips <- get(paste0(my_species, "_trips"))
  species_trips %<>% select(-tripID) %>% distinct()
  species_trips %<>% arrange(pairID, start)
  
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
    summarise(bird1 = bird1[1], bird2 = bird2[1],
              breeding_outcome = rs[1],
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
  
  ## For birds where sex is known ##
  fasting_metrics.knownSex <- species_trips %>%
    filter(!is.na(sex)) %>%
    group_by(pairID, season) %>%
    mutate(bird1 = ring[1],
           bird2 = partner[1],
           n_trips = n(),
           trip_rank = row_number()) %>%
    filter(!(n_trips %% 2 != 0 & trip_rank == max(trip_rank))) %>%
    summarise(debt.F = sum(fasting_debt.hours[sex == "F"], na.rm = T),
              debt.M = sum(fasting_debt.hours[sex == "M"], na.rm = T),
              debt.days = (debt.M - debt.F)/24) %>%
    select(-c(debt.F, debt.M))
  
  pair_behaviour.known_sex <- species_trips %>%
    filter(!is.na(sex)) %>%
    group_by(pairID, season) %>%
    summarise(trip_days.F = median(duration.days[sex == "F"], na.rm = T),
              trip_days.M = median(duration.days[sex == "M"], na.rm = T),
              sd_trip.F = sd(duration.days[sex == "F"], na.rm = T),
              sd_trip.M = sd(duration.days[sex == "M"], na.rm = T),
              male_female_diff.var = sd_trip.M - sd_trip.F, 
              breeding_outcome = rs[1],
              colony = colony[1],
              species = species[1]) %>%
    # Binary breeding outcome
    mutate(breeding_outcome.bin = ifelse(breeding_outcome == "SUCCESSFUL", 1, 0))
  
  pair_behaviour.known_sex <- merge(pair_behaviour.known_sex, fasting_metrics.knownSex,
                                    by = c("pairID", "season"), all.x = TRUE)
  
  ## Output the data
  assign(paste0(my_species, "_pair_behaviour"), pair_behaviour)
  assign(paste0(my_species, "_pair_behaviour.known_sex"), pair_behaviour.known_sex)
  
}

## Make a combined dataframe with dummy species/col variable
pair_behaviour.brooding <- rbind(baBI_pair_behaviour %>% mutate(colony = "birdisland", species = "BBA", speCol = "baBI"), 
                                 baKer_pair_behaviour %>% mutate(colony = "kerguelen", species = "BBA", speCol = "baKer"), 
                                 waBI_pair_behaviour %>% mutate(colony = "birdisland", species = "WAAL", speCol = "waBI"),
                                 waCro_pair_behaviour %>% mutate(colony = "crozet", species = "WAAL", speCol = "waCro"))

pair_behaviour.brooding_withSex <- rbind(baBI_pair_behaviour.known_sex %>% mutate(colony = "birdisland", species = "BBA", speCol = "baBI"), 
                                         baKer_pair_behaviour.known_sex %>% mutate(colony = "kerguelen", species = "BBA", speCol = "baKer"), 
                                         waBI_pair_behaviour.known_sex %>% mutate(colony = "birdisland", species = "WAAL", speCol = "waBI"),
                                         waCro_pair_behaviour.known_sex %>% mutate(colony = "crozet", species = "WAAL", speCol = "waCro"))


save(pair_behaviour.brooding, file = "Data_inputs/all_pair_behaviour_brooding.RData")
save(pair_behaviour.brooding_withSex, file = "Data_inputs/all_pair_behaviour_brooding_withSex.RData")


