# Contrasting associations between parental coordination and breeding success in albatrosses across differing environmental contexts.

Natasha Gillies, Richard A. Phillips, Henri Weimerskirch, Jonathan Potts, Denis Réale, Alastair J. Wilson, Frédéric Angelier, Christophe Barbraud, Ashley Bennison, Karine Delord, Prescillia Lemesle, Samuel Peroteau, Andrew G. Wood, José C. Xavier, Samantha C. Patrick.

## Overview

This repository contains scripts and data to recreate the main results and figures of the following manuscript:

Contrasting associations between parental coordination and breeding success in albatrosses across differing environmental contexts. Gillies N, Phillips RA, Weimerskirch H, Potts J, Réale D, Wilson AJ, Angelier F, Barbraud C, Bennison A, Delord K, Lemesle P, Peroteau S, Wood AG, Xavier JC, Patrick SC

## Scripts

A short description of each script is given below.

-   `ALB_coord_functions`: Contains custom functions required for data processing and analysis.

-   `ALB_coord_data_processing.Rmd*`: *Run this first*. This script processes the trip data into the format required for analysis. Foraging trips tracked during the breeding period for each population are grouped by pair, and pair-level behavioural metrics are derived.

-   `ALB_coord_main_analyses.Rmd`: Runs all analyses presented in the paper and generates the figures and tables for both the main text and supplementary materials. This markdown file also produces the text for the Results section based on the analysis output.

-   `ALB_coord_SOM_compare_variability.Rmd`: Calculates temporal SST variability within the foraging ranges of the four populations, and produces datasets used to create the panels **(B)** and **(C)** in **Figure 1** of the main text.

-   `ALB_coord_SOM_prior_sensitivity.Rmd`: Runs the sensitivity analysis in the SOM, designed to test sensitivity of model outputs to prior specification.

## Data inputs

These data are used in the above scripts. Note that all rings/bird identities have been anonymised and cannot be linked to existing datasets. Please contact the authors if you wish to use these datasets, as we may be able to offer additional information, data, or advice.

### Behaviour/

-   `{species}_{colony}_breeding_trips.RData:` These are the main datasets for analysis; there are four files, one for each population. Each row corresponds to a single foraging trip by the focal bird. The columns are as follows:

    -   *ring*: Factor encoding the unique ID of the focal bird;
    -   *season*: Breeding season during which the foraging trip was recorded. Because albatrosses breed over the austral summer (spanning two calendar years), 'season' refers to the year in which chicks are expected to hatch;
    -   *partner*: Factor encoding the unique ID of the partner bird;
    -   *pairID*: Unique identifier encoding the pair (i.e., combination of ring and partner, in descending alphanumeric order);
    -   *start*: Start time of the foraging trip;
    -   *end*: End time of the foraging trip;
    -   *duration.mins*: Duration of foraging trip in minutes;
    -   *duration.hours*: Duration in hours;
    -   *duration.days*: Duration in days;
    -   *rs*: Factor encoding reproductive success for the pair in that year. One of: FAILED_EGG, FAILED_CHICK, FAILED (*unknown failure*), SUCCESSFUL;
    -   *colony*: Factor encoding the colony of the bird. One of: birdisland, crozet, kerguelen;
    -   *species*: Factor encoding the species of bird. One of: BBAL (black-browed albatross), WAAL (wandering albatross);
    -   *lay_date*: Date on which egg was laid;
    -   *hatch_date*: Date on which chick hatched, if applicable;
    -   *fail_date*: Date on which breeding attempt failed (egg or chick died), if applicable;
    -   *sex*: Sex of bird, if known. One of F (female) or M (male);
    -   *lay_code*: Factor encoding whether lay date was measured (MEAS) or estimated (EST);
    -   *hatch_code*: Factor encoding whether hatch date was measured (MEAS) or estimated (EST);
    -   *phase*: Phase of breeding attempt during which trip fell. One of brooding or incubation;
    -   *days_since_lay*: Numeric variable representing the number of days since the egg was laid. During incubation, negative values may occur if the male began his foraging trip just before the egg was laid.

-   `all_pair_behaviour_{incub/brooding}.RData` This is the output dataset from the data processing script. It forms the primary dataset used for analysis. Each row represents a pair's breeding attempt in a given season and summarises their incubation or brooding behaviour during that period. The columns are as follows:

    -   *pairID*: Unique identifier encoding the pair (i.e., combination of ring and partner, in descending alphanumeric order);
    -   *season*: Breeding season during which the foraging trip was recorded. Because albatrosses breed over the austral summer (spanning two calendar years) 'season' refers to the year in which chicks are expected to hatch;
    -   *breeding_outcome*: Factor encoding reproductive success for the pair in that year. One of: FAILED_EGG, FAILED_CHICK, FAILED (*unknown failure*), SUCCESSFUL;
    -   *colony*: Factor encoding the colony of the bird. One of: birdisland, crozet, kerguelen;
    -   *species*: Factor encoding the species of bird. One of: BBAL (black-browed albatross), WAAL (wandering albatross);
    -   *diff_var.days*: Difference in variability of pair members' foraging trip durations;
    -   *median_trip.days*: Median foraging trip duration of both pair members;
    -   *diff_trip.days*: Difference in pair members' foraging trip durations;
    -   *breeding_outcome.bin*: Binary variable encoding breeding success: 1 for successful, 0 for failed;
    -   *debt.days*: Measure of the overall disparity in investment between the pair members;
    -   *speCol*: Unique identity encoding each population. One of baBI (BBAL at Bird Island), baKer (BBAL at Kerguelen), waBI (WAAL at Bird Island), waCro (WAAL at Crozet).

## General statement on data use

The files included in this repository contain data derived from the long term field projects monitoring individual black-browed albatrosses at Bird Island (South Georgia), Canyon des Sourcils Noirs (Kerguelen archipelago), and Possession Island (Crozet Islands archipelago).

**This is a request to please let us know if you use them.** Several people have spent the best part of their careers collecting the data.

If you plan to analyse the data, there are a number of reasons why it would be very helpful if you could contact Richard Phillips ([raphil\@bas.ac.uk](raphil@bas.ac.uk)) for Bird Island data and Christophe Barbraud ([christophe.barbraud\@cebc.cnrs.fr](christophe.barbraud@cebc.cnrs.fr)) and Karine Delord ([karine.delord\@cebc.cnrs.fr](karine.delord@cebc.cnrs.fr)) for Kerguelen/Crozet data before doing so:

1.  Occasionally we discover and correct errors in the data.

2.  The data are complex and workers who do not know the study system are likely to benefit from advice on interpretation.

3.  At any one time, quite a few people within the existing project collaboration are analysing data from this project. Someone else may already be conducting the analysis you have in mind and it is desirable to prevent duplication of effort.

4.  In order to maintain funding for the project and for further analyses, every few years we submit proposals to funding agencies. It is therefore very helpful for those running the project to know which data analyses are in progress elsewhere.

If you are interested in analysing the detailed project data in any depth, you may find it helpful to have access to the full relational database rather than the files available here. If so, we are always open to further collaboration.

In addition, the long term field project at Kerguelen can be found here: <https://cat.indores.fr/geonetwork/srv/eng/catalog.search#/metadata/44e56d59-a904-4a7d-91c7-fee29aeb0114>

The individual bird identities have been recoded and should therefore not be linked with data archived from other papers which include albatross data from this study site.
