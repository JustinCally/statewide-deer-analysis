This repository contains analysis of modelling deer abundance and distribution in Victoria. The core of the analysis is available in the `statewide_deer_analysis.Rmd`. The data used in this analysis is primarily sourced from ouur database (private access). However, data used in the model has been pre-formatted and is available in `data/multispecies_data.rds`. Additional files and folders in this repository are:  

+ `data/`: Contains various pre-prepared data files including transect counts of signs (created by the `data/species_assign_wrangle.R` script). It also contains data and scripts for detwermining ecological vegetation communities (EVCs) for the model data. The `data/prediction_raster` has the raster we use to predict to across Victoria (public land). And `data/historic_ala_records.rds` is data used in `data/species_assign_wrangle.R`.      

+  
