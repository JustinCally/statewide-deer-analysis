
# Abundance of Deer in Victoria: Statewide Analysis of Camera Trap Data  

This repository contains analysis of modelling deer abundance and distribution in Victoria using camera trap distance sampling (CTDS). The core of the analysis is available in the `statewide_deer_analysis.Rmd` and as an html at https://justincally.github.io/statewide-deer-analysis/. The data used in this analysis is primarily sourced from ouur database (private access). However, data used in the model has been pre-formatted and is available in `data/multispecies_data.rds`. Additional files and folders in this repository are:  

+ `data/`: Contains various pre-prepared data files including transect counts of signs (created by the `data/species_assign_wrangle.R` script). It also contains data and scripts for detwermining ecological vegetation communities (EVCs) for the model data. The `data/prediction_raster` has the raster we use to predict to across Victoria (public land). And `data/historic_ala_records.rds` is data used in `data/species_assign_wrangle.R`.      

+  `functions/` has a suite of functions used in data preparation and then model interegation.  

+  `outputs/` is where models and model outputs such as rasters and plots are located. Note that the stan model draws are not git tracked as they are very large files.  

+  `stan/` is where the stan models used in the analysis are stored.  

+  `docs/` is where the rendered Rmarkdown html and interactive maps are stored. 
