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


### Set colour parameters ------------------------------------------------------

incub_col <- "#4682B4"
brooding_col <- "#FFB347"

ker_col <-  "#156082" 
ker_col.high <- "#155470"
ker_col.low <- "#96d0eb"

cro_col <- "#0F9ED5"
cro_col.low <- "#7bc6e3"
cro_col.high <- "#0080b3"

bi_col <- "#E97132" 
bi_col.low <- "#f59b6c"
bi_col.high <- "#b54b14"

### Set species ----------------------------------------------------------------

species <- c("baBI", "baKer", "waBI", "waCro")


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


### Make some space ---------------------------------------------------------

rm(baBI_trips); rm(baBI_trips.sex); rm(baKer_trips); rm(baKer_trips.sex);
rm(waBI_trips); rm(waBI_trips.sex); rm(waCro_trips); rm(waCro_trips.sex);
rm(baBI_pair_behaviour); rm(baKer_pair_behaviour); rm(waBI_pair_behaviour);
rm(waCro_pair_behaviour); rm(baBI_pair_behaviour.known_sex); 
rm(baKer_pair_behaviour.known_sex); rm(waBI_pair_behaviour.known_sex); 
rm(waCro_pair_behaviour.known_sex); rm(species_trips)

gc()

# SAMPLE SIZES =================================================================

load("Data_inputs/all_pair_behaviour_incub.RData")
load("Data_inputs/all_pair_behaviour_brooding.RData")
load("Data_inputs/all_pair_behaviour_incub_withSex.RData")
load("Data_inputs/all_pair_behaviour_brooding_withSex.RData")

## INCUBATION ##

# Pairs
tapply(pair_behaviour.incub$pairID, pair_behaviour.incub$speCol, n_distinct)
tapply(pair_behaviour.incub_withSex$pairID, pair_behaviour.incub_withSex$speCol, n_distinct)

# Observations
table(pair_behaviour.incub$speCol)
table(pair_behaviour.incub_withSex$speCol)

# Years
tapply(pair_behaviour.incub$season, pair_behaviour.incub$speCol, unique)
tapply(pair_behaviour.incub_withSex$season, pair_behaviour.incub_withSex$speCol, n_distinct)


## BROODING ##

# Pairs
tapply(pair_behaviour.brooding$pairID, pair_behaviour.brooding$speCol, n_distinct)
tapply(pair_behaviour.brooding_withSex$pairID, pair_behaviour.brooding_withSex$speCol, n_distinct)

# Observations
table(pair_behaviour.brooding$speCol)
table(pair_behaviour.brooding_withSex$speCol)

# Years
tapply(pair_behaviour.brooding$season, pair_behaviour.brooding$speCol, unique)
tapply(pair_behaviour.brooding_withSex$season, pair_behaviour.brooding_withSex$speCol, n_distinct)


# CORRELATIONS -----------------------------------------------------------------

# Colony look-up
col_lookup <- data.frame(colID = c("bi", "ker", "cro"),
                         colExp = c("Bird Island", "Kerguelen", "Crozet"))

pair_behaviour.all <- rbind(pair_behaviour.incub %>% mutate(phase = "incubation"), 
                            pair_behaviour.brooding %>% mutate(phase = "brooding"))

dummy_cols <- unique(pair_behaviour.all$speCol)

for (i in 1:length(dummy_cols)) {
  
  which_species <- ifelse(grepl("ba", dummy_cols[i]), "ba", "wa")
  which_col <- tolower(gsub(which_species, "", dummy_cols[i]))
  
  col_name <- col_lookup %>% filter(colID == which_col) %>% pull(colExp)
  spec_name <- ifelse(which_species == "ba", "BBAL", "WAAL")
  
  my_data <- subset(pair_behaviour.all, speCol == dummy_cols[i])
  correlations_plot <- GGally::ggpairs(my_data[,c(6, 7, 8, 10)], title = paste(spec_name, col_name))
  
  
  png(filename = paste0("Figures/rs_failure/correlations_", dummy_cols[i], ".png"), width = 7, height = 6, units = "in", res = 300)
  print(correlations_plot)
  dev.off()
 
}

# BEHAVIOURAL VARIABILTIY ------------------------------------------------------

pair_behaviour.all <- rbind(pair_behaviour.incub %>% mutate(phase = "incubation"), 
                            pair_behaviour.brooding %>% mutate(phase = "brooding"))

dummy_cols <- unique(pair_behaviour.all$speCol)

for (i in 1:length(dummy_cols)) {

  which_species <- ifelse(grepl("ba", dummy_cols[i]), "ba", "wa")
  which_col <- tolower(gsub(which_species, "", dummy_cols[i]))
  
  col_name <- col_lookup %>% filter(colID == which_col) %>% pull(colExp)
  spec_name <- ifelse(which_species == "ba", "BBAL", "WAAL")
  
  my_pairs <- subset(pair_behaviour.all, speCol == dummy_cols[i])
  all_pairs <- subset(pair_behaviour.all, grepl(which_species, speCol))

  if(which_species == "ba") {
  
    maxDebt <- ggplot(aes(x = max_debt.hrs, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$max_debt.hrs), max(all_pairs$max_debt.hrs))
    sumDebt <- ggplot(aes(x = sum_debt.hrs, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$sum_debt.hrs), max(all_pairs$sum_debt.hrs))
    medTrip <- ggplot(aes(x = median_trip.days, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$median_trip.days), max(all_pairs$median_trip.days))
    tripDiff <- ggplot(aes(x = diff_trip.days, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$diff_trip.days), max(all_pairs$diff_trip.days))
    tripVar <- ggplot(aes(x = pair_var.days, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$pair_var.days), max(all_pairs$pair_var.days))
    varDiff <- ggplot(aes(x = diff_var.days, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$diff_var.days), max(all_pairs$diff_var.days))
    
    myplot <- ggarrange(maxDebt, sumDebt, medTrip, tripDiff, tripVar, varDiff,
              ncol = 3, nrow = 2)
    myplot <- annotate_figure(myplot, top = text_grob(paste(spec_name, col_name), 
                                          face = "bold", size = 14))
    png(filename = paste0("Figures/rs_failure/", dummy_cols[i], "_behaviours.png"), width = 9, height = 6, units = "in", res = 300)
    print(myplot)
    dev.off()

  } else {
    
    maxDebt <- ggplot(aes(x = max_debt.hrs, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$max_debt.hrs), max(all_pairs$max_debt.hrs))
    sumDebt <- ggplot(aes(x = sum_debt.hrs, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$sum_debt.hrs), max(all_pairs$sum_debt.hrs))
    medTrip <- ggplot(aes(x = median_trip.days, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$median_trip.days), max(all_pairs$median_trip.days))
    tripDiff <- ggplot(aes(x = diff_trip.days, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$diff_trip.days), max(all_pairs$diff_trip.days))
    tripVar <- ggplot(aes(x = pair_var.days, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$pair_var.days), max(all_pairs$pair_var.days))
    varDiff <- ggplot(aes(x = diff_var.days, fill = phase), data = my_pairs) + geom_histogram(alpha = 0.5, col = "black") + 
      scale_fill_manual(values = c(brooding_col, incub_col)) + theme_classic() + theme(legend.position = "none") +
      xlim(min(all_pairs$diff_var.days), max(all_pairs$diff_var.days))
    
    myplot <- ggarrange(maxDebt, sumDebt, medTrip, tripDiff, tripVar, varDiff,
                        ncol = 3, nrow = 2)
    myplot <- annotate_figure(myplot, top = text_grob(paste(spec_name, col_name), 
                                                      face = "bold", size = 14))
    png(filename = paste0("Figures/rs_failure/", dummy_cols[i], "_behaviours.png"), width = 9, height = 6, units = "in", res = 300)
    print(myplot)
    dev.off()
    
  }
  
  # Summarise variation
  col_maxDebt <- rbind("max debt", 
                        mean(na.omit(my_pairs$max_debt.hrs)),
                        sd(na.omit(my_pairs$max_debt.hrs)),
                        range(na.omit(my_pairs$max_debt.hrs))[1],
                        range(na.omit(my_pairs$max_debt.hrs))[2],
                        IQR(na.omit(my_pairs$max_debt.hrs)) ) 
  
  col_sumDebt <- rbind("summed debt", 
                       mean(na.omit(my_pairs$sum_debt.hrs)),
                       sd(na.omit(my_pairs$sum_debt.hrs)),
                       range(na.omit(my_pairs$sum_debt.hrs))[1],
                       range(na.omit(my_pairs$sum_debt.hrs))[2],
                       IQR(na.omit(my_pairs$sum_debt.hrs)) ) 
  
  col_medTrip <- rbind("median trip", 
                       mean(na.omit(my_pairs$median_trip.days)),
                       sd(na.omit(my_pairs$median_trip.days)),
                       range(na.omit(my_pairs$median_trip.days))[1],
                       range(na.omit(my_pairs$median_trip.days))[2],
                       IQR(na.omit(my_pairs$median_trip.days)) ) 
  
  col_medTripDiff <- rbind("median trip diff", 
                        mean(na.omit(my_pairs$diff_trip.days)),
                          sd(na.omit(my_pairs$diff_trip.days)),
                       range(na.omit(my_pairs$diff_trip.days))[1],
                       range(na.omit(my_pairs$diff_trip.days))[2],
                         IQR(na.omit(my_pairs$diff_trip.days)) ) 
  
  col_tripVar <- rbind("trip variability", 
                            mean(na.omit(my_pairs$pair_var.days)),
                              sd(na.omit(my_pairs$pair_var.days)),
                           range(na.omit(my_pairs$pair_var.days))[1],
                           range(na.omit(my_pairs$pair_var.days))[2],
                             IQR(na.omit(my_pairs$pair_var.days)) ) 
  
  col_tripVarDiff <- rbind("trip variability diff", 
                          mean(na.omit(my_pairs$diff_var.days)),
                            sd(na.omit(my_pairs$diff_var.days)),
                         range(na.omit(my_pairs$diff_var.days))[1],
                         range(na.omit(my_pairs$diff_var.days))[2],
                           IQR(na.omit(my_pairs$diff_var.days)) ) 
  
  
  variables <- cbind(col_maxDebt, col_sumDebt, col_medTrip, 
                     col_medTripDiff, col_tripVar, col_tripVarDiff)
  
  colnames(variables) <- variables[1,]
  variables <- variables[-1,]
  variables <- data.frame(variables)
  variables <- sapply(variables, as.numeric)
  rownames(variables) <- c("mean", "sd", "min", "max", "IQR")
  
  assign(paste0(dummy_cols[i], "_variables"), variables)
 
}

# ______________________________ ####
# ~ * INCUBATION * ~ ###########################################################


# ANALYSIS =====================================================================
    
load("Data_inputs/all_pair_behaviour_incub.RData")

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

# Datasets
behaviour.bba_incub <- pair_behaviour.incub %>% filter(species == "BBA") 
behaviour.waal_incub <- pair_behaviour.incub %>% filter(species == "WAAL")

# Model structure
bf.incub <- brms::bf(breeding_outcome.bin ~
                       debt.days * colony +
                       median_trip.days * colony +
                       diff_trip.days * colony +
                       diff_var.days * colony +
                       (1|season), 
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


# ............................................................ ####
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

my_labels <- c("Trip variability \n(difference)",
               "Trip duration \n(difference)", "Trip duration \n(median)", 
               "Pair debt \n(summed)")

posteriors_plot.incub_bba.horizontal <- 
  ggplot(posteriors.incub_bba, aes(x = value, y = variable)) +
  stat_halfeye(data = subset(posteriors.incub_bba, colony == "Bird Island"), 
               side = "bottom", 
               aes(fill = ifelse(after_stat(x) < 0, "bi_low", "bi_high")),
               height = 0.7) +
  stat_halfeye(data = subset(posteriors.incub_bba, colony == "Kerguelen"), 
               side = "top", 
               aes(fill = ifelse(after_stat(x) < 0, "ker_low", "ker_high")),
               height = 0.7) +
  scale_fill_manual(values = c(
    "bi_low" = bi_col.low, "bi_high" = bi_col.high,
    "ker_low" = ker_col.low, "ker_high" = ker_col.high ) ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, col = "grey20") +
  xlim(-1.5, 2.5) +
  scale_y_discrete(labels = my_labels) +
  labs(y = "", x = "Posterior estimate",
       title = "Black-browed albatrosses") +
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
  xlim(-1.5, 1.2) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none",
        axis.text.y = element_blank())



#### COMBINED posterior plots -----------------------------------------------------------
png(file = "Figures/rs_failure/posterior_estimates_incub_horizontal.png", 
    width = 15, height = 9, units = "in", res = 600)
ggarrange(posteriors_plot.incub_bba.horizontal + theme(plot.margin = unit(c(1,0.05,1,3), "cm")), 
          posteriors_plot.incub_waal.horizontal + theme(plot.margin = unit(c(1,1,1,0.75), "cm")), 
          ncol = 2,
          widths = c(1, 0.9))
dev.off()

# Laptop
png(file = "Figures/rs_failure/posterior_estimates_incub_horizontal2.png", 
    width = 15, height = 9, units = "in", res = 100)
ggarrange(posteriors_plot.incub_bba.horizontal, posteriors_plot.incub_waal.horizontal,
          ncol = 2, widths = c(1, 0.85))
dev.off()



### Conditional effects  -------------------------------------------------------

effects <- 6:9

# Set up variable labeller
variable_labels <- c(
  "debt.days:colony" = "Summed pair debt (days)",
  "median_trip.days:colony" = "Median trip duration (pair; days)",
  "diff_trip.days:colony" = "Median trip duration (pair difference; days)",
  "diff_var.days:colony" = "Trip variability (pair difference; days)"
)

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
    labs(x = var_label, y = "P|Breeding success") +
    theme_bw() +
    theme(text = element_text(size = 16, family = "Calibri"),
          legend.position = "none")
  
  assign(paste0("cond_plot.", var_name.data, ".bba_incub"), cond_plot)
  
  # png(file = paste0("Figures/rs_failure/conditional_estimate_bba_incub_", var_name.data, ".png"), 
  #     width = 7, height = 6, units = "in", res = 300)
  # print(cond_plot)
  # dev.off()
  
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
    scale_fill_manual(values = c(bi_col, cro_col)) +
    scale_colour_manual(values = c(bi_col, cro_col)) +
    labs(x = var_label, y = "P|Breeding success") +
    theme_bw() +
    theme(text = element_text(size = 16, family = "Calibri"),
          legend.position = "none")
  
  assign(paste0("cond_plot.", var_name.data, ".waal_incub"), cond_plot)
  
  # png(file = paste0("Figures/rs_failure/conditional_estimate_waal_incub_", var_name.data, ".png"), 
  #     width = 7, height = 6, units = "in", res = 300)
  # print(cond_plot)
  # dev.off()
  
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

# Get prop above/below zero
below_zero <- data.frame(apply(posterior_samples.incub_bba[,1:10], 2, below_zero.func))
below_zero <- below_zero[,1]

above_zero <- data.frame(apply(posterior_samples.incub_bba[,1:10], 2, above_zero.func))
above_zero <- above_zero[,1]

## Bind everything together
bbal_incub.fixef <- cbind(bbal_incub.fixef, exp_coef)
bbal_incub.fixef$below_zero <- round(below_zero, digits = 2)
bbal_incub.fixef$above_zero <- round(above_zero, digits = 2)

## Get Kerguelen estimates separately
bbal_incub.fixef_extra <- bbal_incub.fixef %>% 
  mutate(Est.extra = NA, LCI.extra = NA, UCI.extra = NA, below.extra = NA, above.extra = NA)

### Process each variable
ker_debt <- process_interaction_estimates("b_debt.days", bbal_incub.fixef, posterior_samples.incub_bba, "b_debt.days:colonykerguelen")
ker_median_trip <- process_interaction_estimates("b_median_trip.days", bbal_incub.fixef, posterior_samples.incub_bba, "b_colonykerguelen:median_trip.days")
ker_diff_trip <- process_interaction_estimates("b_diff_trip.days", bbal_incub.fixef, posterior_samples.incub_bba, "b_colonykerguelen:diff_trip.days")
ker_diff_var <- process_interaction_estimates("b_diff_var.days", bbal_incub.fixef, posterior_samples.incub_bba, "b_colonykerguelen:diff_var.days")

# Combine the results with the original fixef data
ker_fixef <- rbind(ker_debt, ker_median_trip, ker_diff_trip, ker_diff_var)
colnames(ker_fixef) <- colnames(bbal_incub.fixef_extra)

bbal_incub.fixef_extra <- rbind(bbal_incub.fixef_extra, ker_fixef)
bbal_incub.fixef_extra


# Look at mean effects
incub_brms.bba %>%
  spread_draws(b_Intercept, b_colonykerguelen, b_sum_debt.hrs,
               b_median_trip.days, b_diff_trip.days, b_pair_var.days,
               b_diff_var.days, `b_colonykerguelen:sum_debt.hrs`,
               `b_colonykerguelen:median_trip.days`, `b_colonykerguelen:diff_trip.days`,
               `b_colonykerguelen:pair_var.days`, `b_colonykerguelen:diff_var.days`) %>%
  gather(variable, value, b_Intercept:`b_colonykerguelen:diff_var.days`) %>%
  group_by(variable) %>%
  summarize(
    mean = mean(value),
    sd = sd(value),
    lwr = quantile(value, .05),
    upr = quantile(value, .95) )

# Get posterior estimates
brms.incub_posterior <- posterior_samples(brms.incub, "^b")

brms.incub_posterior.df <- stack(brms.incub_posterior) %>%
  rename(variable = ind, value = values) %>%
  relocate(variable, .before = value) %>%
  filter(variable != "b_Intercept")


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

# Get prop above/below zero
below_zero <- data.frame(apply(posterior_samples.incub_waal[,1:10], 2, below_zero.func))
below_zero <- below_zero[,1]

above_zero <- data.frame(apply(posterior_samples.incub_waal[,1:10], 2, above_zero.func))
above_zero <- above_zero[,1]

## Bind everything together
waal_incub.fixef <- cbind(waal_incub.fixef, exp_coef)
waal_incub.fixef$below_zero <- round(below_zero, digits = 2)
waal_incub.fixef$above_zero <- round(above_zero, digits = 2)

## Get Crozet estimates separately
waal_incub.fixef_extra <- waal_incub.fixef %>% 
  mutate(Est.extra = NA, LCI.extra = NA, UCI.extra = NA, below.extra = NA, above.extra = NA)

### Process each variable
cro_debt <- process_interaction_estimates("b_debt.days", waal_incub.fixef, posterior_samples.incub_waal, "b_debt.days:colonycrozet")
cro_median_trip <- process_interaction_estimates("b_median_trip.days", waal_incub.fixef, posterior_samples.incub_waal, "b_colonycrozet:median_trip.days")
cro_diff_trip <- process_interaction_estimates("b_diff_trip.days", waal_incub.fixef, posterior_samples.incub_waal, "b_colonycrozet:diff_trip.days")
cro_diff_var <- process_interaction_estimates("b_diff_var.days", waal_incub.fixef, posterior_samples.incub_waal, "b_colonycrozet:diff_var.days")

# Combine the results with the original fixef data
cro_fixef <- rbind(cro_debt, cro_median_trip, cro_diff_trip, cro_diff_var)
colnames(cro_fixef) <- colnames(waal_incub.fixef_extra)

waal_incub.fixef_extra <- rbind(waal_incub.fixef, cro_fixef)
waal_incub.fixef_extra



# Look at mean effects
incub_brms.waal %>%
  spread_draws(b_Intercept, b_max_debt.hrs, b_colonycrozet, b_sum_debt.hrs,
               b_median_trip.days, b_diff_trip.days, b_pair_var.days,
               b_diff_var.days, `b_max_debt.hrs:colonycrozet`, `b_colonycrozet:sum_debt.hrs`,
               `b_colonycrozet:median_trip.days`, `b_colonycrozet:diff_trip.days`,
               `b_colonycrozet:pair_var.days`, `b_colonycrozet:diff_var.days`) %>%
  gather(variable, value, b_Intercept:`b_colonycrozet:diff_var.days`) %>%
  group_by(variable) %>%
  summarize(
    mean = mean(value),
    sd = sd(value),
    lwr = quantile(value, .05),
    upr = quantile(value, .95) )


# ______________________________ ####
# ~ * BROODING * ~ ###########################################################

# ANALYSIS =====================================================================

load("Data_inputs/all_pair_behaviour_brooding.RData")

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

# Datasets
behaviour.bba_brooding <- pair_behaviour.brooding %>% filter(species == "BBA") #%>% filter(pair_var.days < 4 & max_debt.hrs > 0)
behaviour.waal_brooding <- pair_behaviour.brooding %>% filter(species == "WAAL")

# Model structure
bf.brooding <- brms::bf(breeding_outcome.bin ~
                       debt.days * colony +
                       median_trip.days * colony +
                       diff_trip.days * colony +
                       diff_var.days * colony +
                       (1|season), 
                     family = "bernoulli")

## Fit models
brooding_brms.bba <- brm(bf.brooding,
                            data = behaviour.bba_brooding, 
                            cores = 4, chains = 4, 
                            iter = 20000, warmup = 10000, thin = 10,
                            control = list(adapt_delta = 0.9999, max_treedepth = 14),
                            prior = priors.bba)

pp_check(brooding_brms.bba)
save(brooding_brms.bba, file = "Data_outputs/bba_brooding_brms_model.RData")


brooding_brms.waal <- brm(bf.brooding,
                       data = behaviour.waal_brooding, 
                       cores = 4, chains = 4, 
                       iter = 20000, warmup = 10000, thin = 10,
                       control = list(adapt_delta = 0.9999, max_treedepth = 14),
                       prior = priors.waal)

summary(brooding_brms.waal)
pp_check(brooding_brms.waal)
save(brooding_brms.waal, file = "Data_outputs/waal_brooding_brms_model.RData")

# ............................................................ ####
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
posterior_samples.brooding_bba <- posterior_samples(brooding_brms.bba)

# Extract posterior samples by colony
predictors <- colnames(posterior_samples.brooding_bba)[2:6]
predictors <- predictors[-2]
interaction_terms <- colnames(posterior_samples.brooding_bba)[7:10]
posteriors.brooding_bba <- data.frame()

for (i in 1:length(predictors)) {
  
  # Extract the main effect (Bird Island)
  bird_est <- posterior_samples.brooding_bba[[predictors[i]]]
  
  # Extract the interaction effect and add it to the main effect for Kerguelen
  kerg_est <- bird_est + posterior_samples.brooding_bba[[interaction_terms[i]]]
  
  temp_df <- data.frame(
    variable = predictors[i],
    value = c(bird_est, kerg_est),
    colony = rep(c("Bird Island", "Kerguelen"), each = length(bird_est))
  )
  
  posteriors.brooding_bba <- rbind(posteriors.brooding_bba, temp_df)
  
}

## Relevel factors
posteriors.brooding_bba %<>%
  mutate(variable = fct_rev(fct_relevel(variable, 
                                        c("b_debt.days", 
                                          "b_median_trip.days", 
                                          "b_diff_trip.days", 
                                          "b_diff_var.days"))))


##### Plots ----------------------------------------------------------------

my_labels <- c("Trip variability \n(difference)",
               "Trip duration \n(difference)", "Trip duration \n(median)", 
               "Pair debt \n(summed)")

posteriors_plot.brooding_bba.horizontal <- 
  ggplot(posteriors.brooding_bba, aes(x = value, y = variable)) +
  stat_halfeye(data = subset(posteriors.brooding_bba, colony == "Bird Island"), 
               side = "bottom", 
               aes(fill = ifelse(after_stat(x) < 0, "bi_low", "bi_high")),
               height = 0.7) +
  stat_halfeye(data = subset(posteriors.brooding_bba, colony == "Kerguelen"), 
               side = "top", 
               aes(fill = ifelse(after_stat(x) < 0, "ker_low", "ker_high")),
               height = 0.7) +
  scale_fill_manual(values = c(
    "bi_low" = bi_col.low, "bi_high" = bi_col.high,
    "ker_low" = ker_col.low, "ker_high" = ker_col.high ) ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, col = "grey20") +
  scale_y_discrete(labels = my_labels) +
  labs(y = "", x = "Posterior estimate",
       title = "Black-browed albatrosses") +
  xlim(-16, 15) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none")




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

#### Plots ---------------------------------------------------------------------

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
  #xlim(-3.2, 1.5) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none",
        axis.text.y = element_blank())


### Combined posterior plots -----------------------------------------------------------

png(file = "Figures/rs_failure/posterior_estimates_brooding_horizontal.png", 
    width = 15, height = 9, units = "in", res = 600)
ggarrange(posteriors_plot.brooding_bba.horizontal + theme(plot.margin = unit(c(1,0.05,1,3), "cm")), 
          posteriors_plot.brooding_waal.horizontal + theme(plot.margin = unit(c(1,1,1,0.75), "cm")), 
          ncol = 2,
          widths = c(1, 0.9))
dev.off()

# Laptop #
# png(file = "Figures/rs_failure/posterior_estimates_brooding_horizontal.png", width = 15, height = 9, units = "in", res = 600)
# ggarrange(posteriors_plot.brooding_bba.horizontal, posteriors_plot.brooding_waal.horizontal, ncol = 2,
#           widths = c(1, 0.85))
# dev.off()



## Conditional effects  -------------------------------------------------------

effects <- 6:9

# Set up variable labeller
variable_labels <- c(
  "debt.days:colony" = "Summed pair debt (days)",
  "median_trip.days:colony" = "Median trip duration (pair; days)",
  "diff_trip.days:colony" = "Median trip duration (pair difference; days)",
  "diff_var.days:colony" = "Trip variability (pair difference; days)"
)

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
    labs(x = var_label, y = "P|Breeding success") +
    theme_bw() +
    theme(text = element_text(size = 16, family = "Calibri"),
          legend.position = "none")
  
  assign(paste0("cond_plot.", var_name.data, ".bba_brooding"), cond_plot)
  
  # png(file = paste0("Figures/rs_failure/conditional_estimate_bba_brooding_", var_name.data, ".png"), 
  #     width = 7, height = 6, units = "in", res = 300)
  # print(cond_plot)
  # dev.off()
  
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
    scale_fill_manual(values = c(bi_col, cro_col)) +
    scale_colour_manual(values = c(bi_col, cro_col)) +
    labs(x = var_label, y = "P|Breeding success") +
    theme_bw() +
    theme(text = element_text(size = 16, family = "Calibri"),
          legend.position = "none")
  
  assign(paste0("cond_plot.", var_name.data, ".waal_brooding"), cond_plot)
  
  # png(file = paste0("Figures/rs_failure/conditional_estimate_waal_brooding_", var_name.data, ".png"), 
  #     width = 7, height = 6, units = "in", res = 300)
  # print(cond_plot)
  # dev.off()
  
}

close(pb)


# ............................................................ ####
# INTERPRETATION --------------------------------------------------------------

load("Data_inputs/all_pair_behaviour_brooding.RData")

load("Data_outputs/bba_brooding_brms_model.RData")
load("Data_outputs/waal_brooding_brms_model.RData")

### BBAL ---------------------------------------------------------------------
summary(brooding_brms.bba)
posterior_samples.brooding_bba <- posterior_samples(brooding_brms.bba)

#### Overall summary table ------------------------------
bbal_brooding.fixef <- round(summary(brooding_brms.bba)$fixed, digits = 2)
bbal_brooding.fixef[,c(6,7)] <- round(bbal_brooding.fixef[,c(6,7)], digits = 0)

# Get odds ratios
exp_coef <- round(exp(bbal_brooding.fixef[,c(1, 3, 4)]), digits = 2)
colnames(exp_coef) <- c("Estimate_exp", "LCI_exp", "UCI_exp")
rownames(exp_coef) <- NULL

# Get prop above/below zero
below_zero <- data.frame(apply(posterior_samples.brooding_bba[,1:10], 2, below_zero.func))
below_zero <- below_zero[,1]

above_zero <- data.frame(apply(posterior_samples.brooding_bba[,1:10], 2, above_zero.func))
above_zero <- above_zero[,1]

## Bind everything together
bbal_brooding.fixef <- cbind(bbal_brooding.fixef, exp_coef)
bbal_brooding.fixef$below_zero <- round(below_zero, digits = 2)
bbal_brooding.fixef$above_zero <- round(above_zero, digits = 2)

## Get Kerguelen estimates separately
bbal_brooding.fixef_extra <- bbal_brooding.fixef %>% 
  mutate(Est.extra = NA, LCI.extra = NA, UCI.extra = NA, below.extra = NA, above.extra = NA)

### Process each variable
ker_debt <- process_interaction_estimates("b_debt.days", bbal_brooding.fixef, posterior_samples.brooding_bba, "b_debt.days:colonykerguelen")
ker_median_trip <- process_interaction_estimates("b_median_trip.days", bbal_brooding.fixef, posterior_samples.brooding_bba, "b_colonykerguelen:median_trip.days")
ker_diff_trip <- process_interaction_estimates("b_diff_trip.days", bbal_brooding.fixef, posterior_samples.brooding_bba, "b_colonykerguelen:diff_trip.days")
ker_diff_var <- process_interaction_estimates("b_diff_var.days", bbal_brooding.fixef, posterior_samples.brooding_bba, "b_colonykerguelen:diff_var.days")

# Combine the results with the original fixef data
ker_fixef <- rbind(ker_debt, ker_median_trip, ker_diff_trip, ker_diff_var)
colnames(ker_fixef) <- colnames(bbal_brooding.fixef_extra)

bbal_brooding.fixef_extra <- rbind(bbal_brooding.fixef_extra, ker_fixef)
bbal_brooding.fixef_extra


# Look at mean effects
knitr::kable(brms.brooding %>%
               spread_draws(b_Intercept, b_sd_trip.F, b_sd_trip.M, b_median_trip.M, b_median_trip.F,
                            b_mean_debt.F, b_mean_debt.M) %>%
               gather(variable, value, b_Intercept:b_mean_debt.M) %>%
               group_by(variable) %>%
               summarize(mean = mean(value),
                         sd = sd(value),
                         lwr = quantile(value, .05),
                         upr = quantile(value, .95)), "simple")

# Get posterior estimates
brms.brooding_posterior <- posterior_samples(brms.brooding, "^b")

brms.brooding_posterior.df <- stack(brms.brooding_posterior) %>%
  rename(variable = ind, value = values) %>%
  relocate(variable, .before = value) %>%
  filter(variable != "b_Intercept")


### WAAL ---------------------------------------------------------------------
summary(brooding_brms.waal)
posterior_samples.brooding_waal <- posterior_samples(brooding_brms.waal)

#### Overall summary table ------------------------------

waal_brooding.fixef <- round(summary(brooding_brms.waal)$fixed, digits = 2)
waal_brooding.fixef[,c(6,7)] <- round(waal_brooding.fixef[,c(6,7)], digits = 0)

# Get odds ratios
exp_coef <- round(exp(waal_brooding.fixef[,c(1, 3, 4)]), digits = 2)
colnames(exp_coef) <- c("Estimate_exp", "LCI_exp", "UCI_exp")
rownames(exp_coef) <- NULL

# Get prop above/below zero
below_zero <- data.frame(apply(posterior_samples.brooding_waal[,1:10], 2, below_zero.func))
below_zero <- below_zero[,1]

above_zero <- data.frame(apply(posterior_samples.brooding_waal[,1:10], 2, above_zero.func))
above_zero <- above_zero[,1]

## Bind everything together
waal_brooding.fixef <- cbind(waal_brooding.fixef, exp_coef)
waal_brooding.fixef$below_zero <- round(below_zero, digits = 2)
waal_brooding.fixef$above_zero <- round(above_zero, digits = 2)

## Get Crozet estimates separately
waal_brooding.fixef_extra <- waal_brooding.fixef %>% 
  mutate(Est.extra = NA, LCI.extra = NA, UCI.extra = NA, below.extra = NA, above.extra = NA)

### Process each variable
cro_debt <- process_interaction_estimates("b_debt.days", waal_brooding.fixef, posterior_samples.brooding_waal, "b_debt.days:colonycrozet")
cro_median_trip <- process_interaction_estimates("b_median_trip.days", waal_brooding.fixef, posterior_samples.brooding_waal, "b_colonycrozet:median_trip.days")
cro_diff_trip <- process_interaction_estimates("b_diff_trip.days", waal_brooding.fixef, posterior_samples.brooding_waal, "b_colonycrozet:diff_trip.days")
cro_diff_var <- process_interaction_estimates("b_diff_var.days", waal_brooding.fixef, posterior_samples.brooding_waal, "b_colonycrozet:diff_var.days")

# Combine the results with the original fixef data
cro_fixef <- rbind(cro_debt, cro_median_trip, cro_diff_trip, cro_diff_var)
colnames(cro_fixef) <- colnames(waal_brooding.fixef_extra)

waal_brooding.fixef_extra <- rbind(waal_brooding.fixef_extra, cro_fixef)
waal_brooding.fixef_extra

# ++++++++++++++++++++++++ ####
# * FIGURE 1 * COMBINED POSTERIOR PLOTS ========================================

# Add albatross silhouettes
posteriors_plot.incub_bba.horizontal2<- ggdraw() +
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


png(file = "Figures/rs_failure/FIGURE1.png", width = 10, height = 10, units = "in", res = 100)
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


# Laptop #
# posteriors_plot.incub_bba.horizontal2 <- ggdraw() +
#   draw_plot(posteriors_plot.incub_bba.horizontal + labs(y = "Incubation") +
#               theme(axis.title.y = element_text(margin = margin(r = 10)),
#                     axis.title.x = element_blank(),
#                     text = element_text(size = 16, family = "Calibri"))) +
#   draw_image(file.path("Figures/rs_failure/bba_standing_silhouette.png"),
#              scale = 0.15, x = 0.39, y = 0.28) 
# 
# 
# posteriors_plot.incub_waal.horizontal2 <- ggdraw() +
#   draw_plot(posteriors_plot.incub_waal.horizontal +
#               theme(axis.title.x = element_blank(),
#                     text = element_text(size = 16, family = "Calibri"))) +
#   draw_image(file.path("Figures/rs_failure/waal_standing_silhouette.png"),
#              scale = 0.2, x = 0.38, y = 0.32)
# 
# 
# png(file = "Figures/rs_failure/FIGURE1.png", width = 10, height = 10, units = "in", res = 100)
# ggarrange(posteriors_plot.incub_bba.horizontal2,
#           posteriors_plot.incub_waal.horizontal2,
#           posteriors_plot.brooding_bba.horizontal + 
#             labs(y = "Brooding") + 
#             theme(axis.title.y = element_text(margin = margin(r = 10)),
#                   plot.title = element_blank()),
#           posteriors_plot.brooding_waal.horizontal + 
#             theme(plot.title = element_blank()),
#           ncol = 2,
#           nrow = 2,
#           widths = c(1, 0.8))
# dev.off()


# * FIGURE 2 * BBAL INCUB CONDITIONAL PLOTS ==========================================

png(file = "Figures/rs_failure/FIGURE2.png", width = 15, height = 10, units = "in", res = 300)
ggarrange(
  cond_plot.debt.days.bba_incub +
    labs(tag = "A", y = "Incubation \nP|Breeding success") +
    theme(axis.title.y = element_text(margin = margin(r = 35)), # pc: r = 35 ; laptop 10
          plot.tag = element_text(size = 16, face = "bold"),
          plot.tag.position = c(0.1, 1.01),
          plot.margin = unit(c(0.6, 0.2, 0.1, 0.5), "cm")),
  cond_plot.median_trip.days.bba_incub + 
    labs(tag = "B") + 
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag = element_text(size = 16, face = "bold")),
  cond_plot.diff_trip.days.bba_incub + 
    labs(tag = "C") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.tag = element_text(size = 16, face = "bold")),
  cond_plot.diff_var.days.bba_incub + 
    labs(tag = "D", y = "Brooding \nP|Breeding success") +
    theme(axis.title.y = element_text(margin = margin(r = 35)),
          plot.tag = element_text(size = 16, face = "bold"),
          plot.tag.position = c(0.1, 1.01),
          plot.margin = unit(c(0.6, 0.2, 0.1, 0.5), "cm")),
  ncol = 2,
  nrow = 2,
  widths = c(1, 0.85, 0.85)# pc: c(1, 0.825, 0.825) laptop =0.85
)
dev.off()


# LAPTOP #
cond_plot.debt.days.bba_incub2<- ggdraw() +
  draw_plot(cond_plot.debt.days.bba_incub + labs(tag = "A",y = "P|Breeding success") +
              theme(axis.title.y = element_text(margin = margin(r = 10)),
                    plot.tag = element_text(size = 16, face = "bold"),
                    plot.tag.position = c(0.1, 1.01),
                    plot.margin = unit(c(0.6, 0.2, 0.1, 0.5), "cm"),
                    text = element_text(size = 16, family = "Calibri"))) +
  draw_image(file.path("Figures/rs_failure/bba_standing_silhouette.png"),
             scale = 0.15, x = 0.39, y = 0.28) 

png(file = "Figures/rs_failure/FIGURE2.png", width = 10, height = 10, units = "in", res = 300)
ggarrange(
  cond_plot.debt.days.bba_incub2,# +
    # labs(tag = "A", y = "P|Breeding success") +
    # theme(axis.title.y = element_text(margin = margin(r = 10)),
    #       plot.tag = element_text(size = 16, face = "bold"),
    #       plot.tag.position = c(0.1, 1.01),
    #       plot.margin = unit(c(0.6, 0.2, 0.1, 0.5), "cm")),
  cond_plot.median_trip.days.bba_incub + 
    labs(tag = "B") + 
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag = element_text(size = 16, face = "bold")),
  cond_plot.diff_trip.days.bba_incub + 
    labs(tag = "C", y = "P|Breeding success") +
    theme(axis.title.y = element_text(margin = margin(r = 10)),
          plot.tag = element_text(size = 16, face = "bold"),
          plot.tag.position = c(0.1, 1.01),
          plot.margin = unit(c(0.6, 0.2, 0.1, 0.5), "cm")),
  cond_plot.diff_var.days.bba_incub + 
    labs(tag = "D") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.tag = element_text(size = 16, face = "bold")),
  ncol = 2,
  nrow = 2,
  widths = c(1, 0.85, 0.85) # pc: c(1, 0.825, 0.825) laptop =0.85
)
dev.off()

# * FIGURE 3 * MEDIAN TRIP DURATIONG BROODING ==========================================

# Add albatross silhouettes
median_trip.bba_brooding.plot <- ggdraw() +
  draw_plot(cond_plot.median_trip.days.bba_brooding +
              labs(tag = "A") +
              theme(axis.title.y = element_text(margin = margin(r = 10)),
                    plot.tag = element_text(size = 16, face = "bold"),
                    plot.tag.position = c(0, 1.01),
                    plot.margin = unit(c(0.6, 0.2, 0.1, 0.5), "cm"))) +
  draw_image(file.path("Figures/rs_failure/bba_standing_silhouette.png"),
             scale = 0.15, x = 0.39, y = 0.28) 


median_trip.waal_brooding.plot <- ggdraw() +
  draw_plot(cond_plot.median_trip.days.waal_brooding +
              labs(tag = "B") +
              theme(axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    text = element_text(size = 16, family = "Calibri"),
                    plot.tag = element_text(size = 16, face = "bold"),
                    plot.tag.position = c(-0.025, 1.01),
                    plot.margin = unit(c(0.6, 0.2, 0.1, 0.5), "cm"))) +
  draw_image(file.path("Figures/rs_failure/waal_standing_silhouette.png"),
             scale = 0.2, x = 0.38, y = 0.32)



png(file = "Figures/rs_failure/FIGURE3.png", width = 12, height = 6, units = "in", res = 300)
ggarrange(
  median_trip.bba_brooding.plot,
  median_trip.waal_brooding.plot,
  nrow = 1, ncol = 2,
  widths = c(1, 0.85)
)

dev.off()



# * TABLE 1 * BBAL COEFFICIENTS ================================================

summary(incub_brms.bba); summary(brooding_brms.bba)

bbal_incub.fixef_MS <- bbal_incub.fixef %>%
  mutate(est_error = paste0(Estimate, " ± ", Est.Error),
         cis = paste0('[', `l-95% CI`, ", ", `u-95% CI`, ']')) %>%
  select(est_error, cis, below_zero, above_zero)

bbal_brooding.fixef_MS <- bbal_brooding.fixef %>%
  mutate(est_error = paste0(Estimate, " ± ", Est.Error),
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

# ++++++++++++++++++++++++ ####
# ______________________________ ####
# ~ * APPENDIX * ~ ###########################################################

# INCUBATION =====================================================================

load("Data_inputs/all_pair_behaviour_incub_withSex.RData")

#### Specify priors ---------------------------------------------------------------

# BBA 
priors.bba <- c(set_prior("normal(-0.5, 0.5)", class = "b", coef = "male_female_diff.days"),
                set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.F"),
                set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.M"),       
                set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.F"),
                set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.M"),
                set_prior("normal(-0.5, 0.5)", class = "b", coef = "male_female_diff.var") )

# WAAL
priors.waal <-  c(set_prior("normal(-0.5, 0.5)", class = "b", coef = "male_female_diff.days"),
                  set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.F"),
                  set_prior("normal(-0.25, 0.5)", class = "b", coef = "trip_days.M"),       
                  set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.F"),
                  set_prior("normal(-0.12, 0.5)", class = "b", coef = "sd_trip.M"),
                  set_prior("normal(-0.5, 0.5)", class = "b", coef = "male_female_diff.var") )

#### Specify models ----------------------------------------------------------------

# Datasets
behaviour.bba.sex_incub <- pair_behaviour.incub_withSex %>% filter(species == "BBA")
behaviour.waal.sex_incub <- pair_behaviour.incub_withSex %>% filter(species == "WAAL")

# Model structure
bf.incub.sex <- brms::bf(breeding_outcome.bin ~
               male_female_diff.days * colony +
                         trip_days.F * colony +
                         trip_days.M * colony +
                           sd_trip.F * colony + 
                           sd_trip.M * colony +
                male_female_diff.var * colony +
                       (1|season), 
                     family = "bernoulli")


## Fit models
incub_brms.bba.sex <- brm(bf.incub.sex,
                      data = behaviour.bba.sex_incub, 
                      cores = 3, chains = 3, 
                      iter = 6000, warmup = 2000, thin = 10,
                      control = list(adapt_delta = 0.9999, max_treedepth = 14),
                      prior = priors.bba)

summary(incub_brms.bba.sex)
pp_check(incub_brms.bba.sex, ndraw = 50)
save(incub_brms.bba.sex, file = "Data_outputs/bba_incub_sex_brms_model.RData")

round(fixef(incub_brms.bba.sex), digits = 2)
round(exp(fixef(incub_brms.bba.sex)), digits = 2)


incub_brms.waal.sex <- brm(bf.incub,
                       data = behaviour.waal.sex_incub, 
                       cores = 3, chains = 3, 
                       iter = 6000, warmup = 2000, thin = 10,
                       control = list(adapt_delta = 0.9999, max_treedepth = 14),
                       prior = priors.waal)

summary(incub_brms.waal.sex)
pp_check(incub_brms.waal.sex, ndraw = 50)
save(incub_brms.waal.sex, file = "Data_outputs/waal_incub_sex_brms_model.RData")

round(fixef(incub_brms.waal.sex), digits = 2)
round(exp(fixef(incub_brms.waal.sex)), digits = 2)





# OLD CODE ----------------------------------------------------------------



## Split by sex
posteriors_plot <- brms.incub_posterior.df %>%
  mutate(sex = ifelse(grepl("F", variable), "F", "M")) %>%
  mutate(variable = case_when(variable == "b_sd_trip.M" ~ "Trip variability",
                              variable == "b_sd_trip.F" ~ "Trip variability",
                              variable == "b_median_trip.M" ~ "Median trip duration",
                              variable == "b_median_trip.F" ~ "Median trip duration",
                              variable == "b_mean_debt.F" ~ "Fasting debt",
                              variable == "b_mean_debt.M" ~ "Fasting debt"),
         variable = fct_relevel(variable, c("Trip variability", 
                                            "Median trip duration", 
                                            "Fasting debt"))) %>%
  ggplot(aes(x = value, y = variable, fill = after_stat(x < 0))) + 
  stat_halfeye() +
  facet_wrap(.~sex,
             labeller = as_labeller(c("F" = "Female",
                                      "M" = "Male"))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#FFD2C8", "#E5F7FF")) +
  labs(y = "", x = "Posterior estimate",
       caption = "Data: pair.df in {analyse_behavioural_impacts.R}. Intercept omitted for readability.") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.x.top = element_text(face = "bold"),
        plot.caption = element_text(face = "italic",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 17))

png(file = paste0("Figures/", species, "_posteriors.png"), width = 9, height = 6, units = "in", res = 300)
posteriors_plot
dev.off()


# Coefficient range
brms.incub %>%
  tidy_draws() %>%
  select(.chain:b_sd_trip.M) %>%
  group_by(.chain, .iteration, .draw) %>%
  gather(var, val, b_Intercept:ncol(.)) %>%
  ungroup() %>%
  group_by(var) %>%
  summarize(mean = mean(val),
            sd = sd(val),
            lwr = quantile(val, .05),
            upr = quantile(val, .95)) %>%
  filter(var != "b_Intercept") %>%
  mutate(var = case_when(var == "b_sd_trip.M" ~ "Male trip variability",
                         var == "b_sd_trip.F" ~ "Female trip variability",
                         var == "b_median_trip.M" ~ "Male median trip duration",
                         var == "b_median_trip.F" ~ "Female median trip duration",
                         var == "b_mean_debt.F" ~ "Female debt",
                         var == "b_mean_debt.M" ~ "Male debt")) %>%
  mutate(var = fct_relevel(var, c("Male trip variability", "Female trip variability", 
                            "Male median trip duration", "Female median trip duration",
                            "Female debt", "Male debt"))) %>%
  ggplot(.,aes(var, mean, ymin=lwr, ymax=upr)) +
  coord_flip() +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(x = "", y = "Mean coefficient (with 95% Intervals)",
       caption = "Data: pair.df in {analyse_behavioural_impacts}. Intercept omitted for readability.") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "italic", hjust = 1),
        plot.caption = element_text(face = "italic"))


## Split by sex
coefficient_plot <- brms.incub %>%
  tidy_draws() %>%
  select(.chain:b_sd_trip.M) %>%
  group_by(.chain, .iteration, .draw) %>%
  gather(var, val, b_Intercept:ncol(.)) %>%
  ungroup() %>%
  group_by(var) %>%
  summarize(mean = mean(val),
            sd = sd(val),
            lwr = quantile(val, .05),
            upr = quantile(val, .95)) %>%
  filter(var != "b_Intercept") %>%
  mutate(sex = ifelse(grepl("F", var), "Female", "Male")) %>%
  mutate(var = case_when(var == "b_sd_trip.M" ~ "Trip variability",
                              var == "b_sd_trip.F" ~ "Trip variability",
                              var == "b_median_trip.M" ~ "Median trip duration",
                              var == "b_median_trip.F" ~ "Median trip duration",
                              var == "b_mean_debt.F" ~ "Fasting debt",
                              var == "b_mean_debt.M" ~ "Fasting debt"),
         var = fct_relevel(var, c("Trip variability", 
                                            "Median trip duration", 
                                            "Fasting debt"))) %>%
  ggplot(.,aes(var, mean, ymin=lwr, ymax=upr)) +
  facet_wrap(.~sex) +
  coord_flip() +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(x = "", y = "Mean coefficient (with 95% credible intervals)",
       caption = "Data: pair.df in {analyse_behavioural_impacts}. Intercept omitted for readability.") +
  theme_bw() +
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        plot.caption = element_text(face = "italic",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        strip.text.x.top = element_text(face = "bold"),
        axis.text.y = element_text(size = 17))

png(file = paste0("Figures/", species, "_coefficients.png"), width = 9, height = 6, units = "in", res = 300)
coefficient_plot
dev.off()


#### Non-interactive posteriors --------------

brms.incub_posterior.waal <- posterior_samples(incub_brms.waal, "^b")

brms.incub_posterior.waal.df <- stack(brms.incub_posterior.waal) %>%
  rename(variable = ind, value = values) %>%
  relocate(variable, .before = value) %>%
  filter(variable != "b_Intercept")

# Pair-level effects
posteriors_incub.waal <- brms.incub_posterior.waal.df %>%
  mutate(variable = case_when(
    variable == "b_max_debt.hrs" ~ "Max pair debt",
    variable == "b_sum_debt.hrs" ~ "Summed pair debt",
    variable == "b_max_debt.hrs" ~ "Max trip duration",
    variable == "b_median_trip.days" ~ "Median trip duration (pair)",
    variable == "b_Imedian_trip.daysE2" ~ "Median trip duration (pair, quadratic)",
    variable == "b_diff_trip.days" ~ "Median trip duration (difference)",
    variable == "b_pair_var.days" ~ "Trip variability (pair)",
    variable == "b_diff_var.days" ~ "Trip variability (difference)",
    variable == "b_colonykerguelen" ~ "Kerguelen"),
    variable = fct_relevel(variable, c("Trip variability (difference)",
                                       "Trip variability (pair)",
                                       "Median trip duration (difference)",
                                       "Median trip duration (pair, quadratic)",
                                       "Median trip duration (pair)",
                                       "Max trip duration",
                                       "Summed pair debt",
                                       "Max pair debt"))) %>%
  filter(variable != "Kerguelen" & variable != "Median trip duration (pair, quadratic)") %>%
  ggplot(aes(x = value, y = variable, fill = after_stat(x < 0))) + 
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#FFD2C8", "#E5F7FF")) +
  xlim(-0.025, 0.04) +
  labs(y = "", x = "Posterior estimate",
       #caption = "Data: behaviour.waal in {analyse_behavioural_impacts.R}. Intercept omitted for readability.",
       title = "Wandering albatrosses") +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Calibri"),
        legend.position = "none",
        strip.text.x.top = element_text(face = "bold"),
        plot.caption = element_text(face = "italic",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_blank(),
        title = element_blank())




### Plot conditional effects  --------------------------------------------------

#### Fasting debt --------------------------------------------------------------

data_term.female_debt <- conditional_effects(brms.incub)$mean_debt.F

female_debt.plot <- ggplot() +
  geom_ribbon(aes(x = mean_debt.F, ymin = lower__, ymax = upper__), data = data_term.female_debt, fill = "grey80", alpha = 0.5) +
  geom_line(aes(x = mean_debt.F, y = estimate__), data = data_term.female_debt, col = "orange", linewidth = 1) +
  geom_point(aes(x = mean_debt.F, y = breeding_outcome.bin), data = pair_behaviour) +
  labs(x = "Female fasting debt (hours)", y = "P|Breeding success") +
  theme_bw()

data_term.male_debt <- conditional_effects(brms.incub)$mean_debt.M

male_debt.plot <- ggplot() +
  geom_ribbon(aes(x = mean_debt.M, ymin = lower__, ymax = upper__), data = data_term.male_debt, fill = "grey80", alpha = 0.5) +
  geom_line(aes(x = mean_debt.M, y = estimate__), data = data_term.male_debt, col = "orange", linewidth = 1) +
  geom_point(aes(x = mean_debt.M, y = breeding_outcome.bin), data = pair_behaviour) +
  labs(x = "Male fasting debt (hours)", y = "P|Breeding success") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())


png(file = paste0("Figures/RSFailure_ConditionalEffects/", species, "_fastingDebt.png"), width = 9, height = 6, units = "in", res = 300)
gridExtra::grid.arrange(female_debt.plot, male_debt.plot,
                        ncol = 2)
dev.off()

#### Median trip duration ----------------------------------------------------------

data_term.female_medianTrip <- conditional_effects(brms.incub)$median_trip.F

female_medianTrip.plot <- ggplot() +
  geom_ribbon(aes(x = median_trip.F, ymin = lower__, ymax = upper__), data = data_term.female_medianTrip, fill = "grey80", alpha = 0.5) +
  geom_line(aes(x = median_trip.F, y = estimate__), data = data_term.female_medianTrip, col = "orange", linewidth = 1) +
  geom_point(aes(x = median_trip.F, y = breeding_outcome.bin), data = pair_behaviour) +
  labs(x = "Female median trip duration (hours)", y = "P|Breeding success") +
  theme_bw()

data_term.male_medianTrip <- conditional_effects(brms.incub)$median_trip.M

male_medianTrip.plot <- ggplot() +
  geom_ribbon(aes(x = median_trip.M, ymin = lower__, ymax = upper__), data = data_term.male_medianTrip, fill = "grey80", alpha = 0.5) +
  geom_line(aes(x = median_trip.M, y = estimate__), data = data_term.male_medianTrip, col = "orange", linewidth = 1) +
  geom_point(aes(x = median_trip.M, y = breeding_outcome.bin), data = pair_behaviour) +
  labs(x = "Male median trip duration (hours)", y = "P|Breeding success") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())


png(file = paste0("Figures/RSFailure_ConditionalEffects/", species, "_medianTrip.png"), width = 9, height = 6, units = "in", res = 300)
gridExtra::grid.arrange(female_variability.plot, male_variability.plot,
                        ncol = 2)
dev.off()


#### Trip variability ----------------------------------------------------------

data_term.female_variability <- conditional_effects(brms.incub)$sd_trip.F

female_variability.plot <- ggplot() +
  geom_ribbon(aes(x = sd_trip.F, ymin = lower__, ymax = upper__), data = data_term.female_variability, fill = "grey80", alpha = 0.5) +
  geom_line(aes(x = sd_trip.F, y = estimate__), data = data_term.female_variability, col = "orange", linewidth = 1) +
  geom_point(aes(x = sd_trip.F, y = breeding_outcome.bin), data = pair_behaviour) +
  labs(x = "Female trip variability (hours)", y = "P|Breeding success") +
  theme_bw()

data_term.male_variability <- conditional_effects(brms.incub)$sd_trip.M

male_variability.plot <- ggplot() +
  geom_ribbon(aes(x = sd_trip.M, ymin = lower__, ymax = upper__), data = data_term.male_variability, fill = "grey80", alpha = 0.5) +
  geom_line(aes(x = sd_trip.M, y = estimate__), data = data_term.male_variability, col = "orange", linewidth = 1) +
  geom_point(aes(x = sd_trip.M, y = breeding_outcome.bin), data = pair_behaviour) +
  labs(x = "Male trip variability (hours)", y = "P|Breeding success") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())


png(file = paste0("Figures/RSFailure_ConditionalEffects/", species, "_tripVariability.png"), width = 9, height = 6, units = "in", res = 300)
gridExtra::grid.arrange(female_variability.plot, male_variability.plot,
                        ncol = 2)
dev.off()





