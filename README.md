# Parental coordination and environment interact to drive breeding success in albatrosses.
Natasha Gillies, Richard A. Phillips, Henri Weimerskirch, Jonathan Potts, Denis Réale, Alastair J. Wilson, Frédéric Angelier, Christophe Barbraud, Ashley Bennison, Karine Delord, Prescillia Lemesle, Samuel Peroteau, Andrew G. Wood, José C. Xavier, Samantha C. Patrick.

## Overview
This repository contains scripts and data to recreate the main results and figures of the following manuscript:

Parental coordination and environment interact to drive breeding success in albatrosses. Gillies N, Phillips RA, Weimerskirch H, Potts J, Réale D, Wilson AJ, Angelier F, Barbraud C, Bennison A, Delord K, Lemesle P, Peroteau S, Wood AG, Xavier JC, Patrick SC

## Scripts
A short description of each script is given below.

- **ALB_FOR_functions** This script contains custom functions required to run the analyses and data porcessing.
- **failureRates_analyses.Rmd** This script runs all the analyses in the paper, and produces the figures and tables for the main text and supplementary materials. The markdown file will generate the text from the results section.
- **failureRates_data_processing.Rmd** This script processes the trip data into the format required for the analyses. Foraging trips tracked during the breeding period for each population are grouped by each pair, and the derived pair-level behavioural metrics are calculated.

## Data inputs - WIP

These data are used in the above scripts. Note that all Rings/BirdIDs have been recoded and so cannot be linked to existing datasets. Please contact the authors if you would like to make use of these datasets, as we may be able to offer additional information, data, or advice. 

- **BAS_demo_complete.RData** This is the main dataset for analysis. Each row corresponds to an individual GPS fix. The columns are as follows:
  -  _Ring_: Factor encoding unique ID of bird
