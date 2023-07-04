#### Model across species ####
library(cmdstanr)
library(dplyr)
library(sf)
library(bayesplot)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
options(mc.cores=6)

# Source internal functions
sapply(list.files("functions", full.names = T), source, verbose = F)
# Species to run model for.
# deer_species_all <- c("Cervus unicolor", "Dama dama", "Cervus elaphus")
deer_species_all <- "Cervus elaphus"
# Projects to select.
project_short_name <- c("hog_deer_2023", "StatewideDeer")
# Buffer for data extraction.
spatial_buffer <- 1000
# Covariate Rasters
raster_files <- "/Volumes/DeerVic\ Photos/Processed_Rasters"
prediction_raster <- "data/prediction_raster/vic_model_data_resampled.tif"
n_max_no_det <- 5
n_max_det <- 50

#### Model formulas ####

transect_formula <- ~Survey
ab_formula <- ~ scale(BIO12)  + scale(BIO04)  + scale(BIO01) + scale(TreeDensity) + scale(sqrt(PastureDistance)) + scale(TWIND) + scale(BIO15) + scale(SLOPE) + scale(MRVBF) + scale(sqrt(ForestEdge))
det_formula <- ~ scale(HerbaceousUnderstoryCover)

#### Database connection ####
con <- weda::weda_connect(password = keyring::key_get(service = "ari-dev-weda-psql-01",
                                                username = "psql_user"), username = "psql_user")

cams_curated <- tbl(con, dbplyr::in_schema("camtrap", "curated_camtrap_operation")) %>%
  dplyr::filter(ProjectShortName %in% !!project_short_name) %>% # Only retrieve cam records for hog deer 2023
  dplyr::collect() %>%
  dplyr::arrange(SiteID)

model_data <- lapply(deer_species_all,
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
                     n_max_det = n_max_det, transects = FALSE)

names(model_data) <- deer_species_all


#### Run models for all species ####
# STAN settings
ni <- 250
nw <- 250
nt <- 1
nb <- 300
nc <- 6

inits = lapply(1:nc, function(i) list(beta_det=runif(2),
                                      beta_trans_det = runif(1),
                                      beta_psi_det = runif(1),
                                      beta_psi = runif((1+length(labels(terms(ab_formula))))),
                                      activ = runif(1),
                                      alpha = 2, rho=1,
                                      phi = runif(1),
                                      site_sd = runif(1),
                                      eps_ngs = 1/3,
                                      site_raw = rnorm(317),
                                      eta = rnorm(317)))

model_negbin_co <- cmdstan_model(here::here("stan", "count_only_negbin_gp.stan"))

fitintegrn_red <- model_negbin_co$sample(data = model_data[["Cervus elaphus"]], chains = nc,
                                 parallel_chains = nc,
                                 show_messages = TRUE,
                                 save_warmup = FALSE,
                                 iter_sampling = ni,
                                 iter_warmup = nw)

fitintegrn_hog$save_object("outputs/fitintegrn_hog.rds")


#### View the point densities ####
rn_dens <- fitintegrn_fallow$summary("N_site")
# Bind the density to the camera information
density_at_sites_rn <- cbind(cams_curated, rn_dens) %>%
  dplyr::mutate(means_sqrt = sqrt(mean)) %>%
  sf::st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4283)

mapview::mapview(density_at_sites_rn, zcol = "mean")

#### View the covariate effects ####
beta_psi_draws <- fitintegrn_hog$draws("beta_psi", format = "matrix") %>% `colnames<-`(c("Intercept", labels(terms(ab_formula))))
mcmc_areas(beta_psi_draws)

#### View the gp effect ####
gp_draws <- fitintegrn_hog$draws(c("rho", "alpha"), format = "matrix")
mcmc_areas(gp_draws)


