# Parental coordination and environment interact to drive breeding success in albatrosses.
Natasha Gillies, Richard A. Phillips, Henri Weimerskirch, Jonathan Potts, Denis Réale, Alastair J. Wilson, Frédéric Angelier, Christophe Barbraud, Ashley Bennison, Karine Delord, Prescillia Lemesle, Samuel Peroteau, Andrew G. Wood, José C. Xavier, Samantha C. Patrick.

## Overview
This repository contains scripts and data to recreate the main results and figures of the following manuscript:

Parental coordination and environment interact to drive breeding success in albatrosses. Gillies N, Phillips RA, Weimerskirch H, Potts J, Réale D, Wilson AJ, Angelier F, Barbraud C, Bennison A, Delord K, Lemesle P, Peroteau S, Wood AG, Xavier JC, Patrick SC

## Scripts
A short description of each script is given below.

- **ALB_FOR_functions**: Contains custom functions required for data processing and analysis.
- **failureRates_data_processing.Rmd**: *Run this first*. This script processes the trip data into the format required for analysis. Foraging trips tracked during the breeding period for each population are grouped by pair, and pair-level behavioural metrics are derived.
- **failureRates_analyses.Rmd** Runs all analyses presented in the paper and generates the figures and tables for both the main text and supplementary materials. This markdown file also produces the text for the Results section based on the analysis output.

## Data inputs

These data are used in the above scripts. Note that all rings/bird identities have been anonymised and cannot be linked to existing datasets. Please contact the authors if you wish to use these datasets, as we may be able to offer additional information, data, or advice.

- **ANON_[species]_[colony]_breedingTrips.RData** These are the main datasets for analysis; there are four files, one for each population. Each row corresponds to a single foraging trip by the focal bird. The columns are as follows:
  -  _ring_: Factor encoding the unique ID of the focal bird;
  -  _season_: Breeding season during which the foraging trip was recorded. Because albatrosses breed over the austral summer—spanning two calendar years—'season' refers to the year in which chicks are expected to hatch;
  -  _partner_: Factor encoding the unique ID of the partner bird;
  -  _pairID_:  Unique identifier encoding the pair (i.e., combination of ring and partner, in descending alphanumeric order);
  -  -start_: Start time of the foraging trip;
  -  _end_: End time of the foraging trip;
  -  _duration.mins_: Duration of foraging trip in minutes;
  -  _duration.hours_: Duration in hours;
  -  _duration.days_: Duration in days;
  -  _rs_: Factor encoding reproductive success for the pair in that year. One of: FAILED_EGG, FAILED_CHICK, FAILED (*unknown failure*), SUCCESSFUL;
  -  _colony_: Factor encoding the colony of the bird. One of: birdisland, crozet, kerguelen;
  -  _species_: Factor encoding the species of bird. One of: BBAL (black-browed albatross), WAAL (wandering albatross);
  -  _lay_date_: Date on which egg was laid;
  -  _hatch_date_: Date on which chick hatched, if applicable;
  -  _fail_date_: Date on which breeding attempt failed (egg or chick died), if applicable;
  -  _sex_: Sex of bird, if known. One of F (female) or M (male);
  -  _lay_code_: Factor encoding whether lay date was measured (MEAS) or estimated (EST);
  -  _hatch_code_: Factor encoding whether hatch date was measured (MEAS) or estimated (EST);
  -  _phase_: Phase of breeding attempt during which trip fell. One of brooding or incubation;
  -  _days_since_lay_:  Numeric variable representing the number of days since the egg was laid. During incubation, negative values may occur if the male began his foraging trip just before the egg was laid.
  
- **all_pair_behaviour_[incub/brooding].RData** This is the output dataset from the data processing script. It forms the primary dataset used for analysis. Each row represents a pair’s breeding attempt in a given season and summarises their incubation or brooding behaviour during that period. The columns are as follows:
  -  _pairID_:  Unique identifier encoding the pair (i.e., combination of ring and partner, in descending alphanumeric order);
  -  _season_: Breeding season during which the foraging trip was recorded. Because albatrosses breed over the austral summer—spanning two calendar years—'season' refers to the year in which chicks are expected to hatch;
  -  _breeding_outcome_: Factor encoding reproductive success for the pair in that year. One of: FAILED_EGG, FAILED_CHICK, FAILED (*unknown failure*), SUCCESSFUL;
  -  _colony_: Factor encoding the colony of the bird. One of: birdisland, crozet, kerguelen;
  -  _species_: Factor encoding the species of bird. One of: BBAL (black-browed albatross), WAAL (wandering albatross);
  -  _diff_var.days_: Difference in variability of pair members' foraging trip durations;
  -  _median_trip.days_: Median foraging trip duration of both pair members;
  -  _diff_trip.days_: Difference in pair members' foraging trip durations;
  -  _breeding_outcome.bin_: Binary variable encoding breeding success: 1 for successful, 0 for failed;
  -  _debt.days_: Measure of the overall disparity in investment between the pair members;
  -  _speCol_: Unique identity encoding each population. One of baBI (BBAL at Bird Island), baKer (BBAL at Kerguelen), waBI (WAAL at Bird Island), waCro (WAAL at Crozet).