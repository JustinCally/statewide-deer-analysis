# Model across species.
# Source internal functions
sapply(list.files("functions", full.names = T), source, verbose = F)
# Species to run model for.
deer_species_all <- c("Cervus unicolor", "Dama dama", "Cervus elaphus")
# Projects to select.
project_short_name <- c("hog_deer_2023", "StatewideDeer")
# Buffer for data extraction.
spatial_buffer <- 1000
# Covariate Rasters
raster_files <- "/Volumes/DeerVic\ Photos/Processed_Rasters"
prediction_raster <- "data/vic_model_data_resampled.tif"
n_max_no_det <- 10
n_max_det <- 30

#### Model formulas ####

transect_formula <- ~Survey
ab_formula <- ~ scale(BIO12)  + scale(BIO04)  + scale(BIO01) + scale(TreeDensity) + scale(sqrt(PastureDistance)) + scale(TWIND) + scale(BIO15) + scale(SLOPE) + scale(MRVBF) + scale(sqrt(ForestEdge))
det_formula <- ~ scale(HerbaceousUnderstoryCover)

#### Database connection ####
con <- weda::weda_connect(password = keyring::key_get(service = "ari-dev-weda-psql-01",
                                                username = "psql_user"), username = "psql_user")

model_data <- lapply(deer_species_all[1],
                     FUN = prepare_model_data,
                     projects = project_short_name,
                     buffer = spatial_buffer,
                     detection_formula = det_formula,
                     abundance_formula = ab_formula,
                     transect_formula = transect_formula,
                     con = con,
                     raster_dir = raster_files,
                     prediction_raster = prediction_raster,
                     n_max_no_det = n_max_no_det,
                     n_max_det = n_max_det)
