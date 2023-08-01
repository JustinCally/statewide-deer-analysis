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
deer_species_all <- c("Cervus unicolor", "Dama dama", "Cervus elaphus", "Axis porcinus")
# Projects to select.
project_short_name <- c("hog_deer_2023", "StatewideDeer")
# Buffer for data extraction.
spatial_buffer <- 1000
# Covariate Rasters
raster_files <- "/Volumes/DeerVic\ Photos/Processed_Rasters"
prediction_raster <- "data/prediction_raster/vic_model_data_resampled.tif"
n_max_no_det <- 10
n_max_det <- 50

#### Model formulas ####
transect_formula <- ~Survey
ab_formula <- ~ scale(BIO12)  + scale(BIO04)  + scale(BIO01) + scale(TreeDensity) + scale(sqrt(PastureDistance)) + scale(TWIND) + scale(BIO15) + scale(SLOPE) + scale(MRVBF) + scale(sqrt(ForestEdge))
# Use an intercept only detection model to start with
det_formula <- ~ 1 #scale(HerbaceousUnderstoryCover)

#### Database connection ####
con <- weda::weda_connect(password = keyring::key_get(service = "ari-dev-weda-psql-01",
                                                username = "psql_user"), username = "psql_user")

cams_curated <- tbl(con, dbplyr::in_schema("camtrap", "curated_camtrap_operation")) %>%
  dplyr::filter(ProjectShortName %in% !!project_short_name) %>% # Only retrieve cam records for hog deer 2023
  dplyr::collect() %>%
  dplyr::arrange(SiteID)

evaluate_transects <- c(TRUE, TRUE, TRUE, FALSE)
model_data <- list()
for(i in 1:length(deer_species_all)) {
model_data[[i]] <- prepare_model_data(species = deer_species_all[i],
                                 projects = project_short_name,
                     buffer = spatial_buffer,
                     detection_formula = det_formula,
                     abundance_formula = ab_formula,
                     transect_formula = transect_formula,
                     con = con,
                     raster_dir = raster_files,
                     prediction_raster = prediction_raster,
                     n_max_no_det = n_max_no_det,
                     n_max_det = n_max_det,
                     evaltransects = evaluate_transects[i],
                     snapshot_interval = 2)

}

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
model_negbin <- cmdstan_model(here::here("stan", "count_det_nondet_negbin_gp.stan"))
# model_poisson <- cmdstan_model(here::here("stan", "count_det_nondet_poisson_gp.stan"))
# model_poisson_co <- cmdstan_model(here::here("stan", "count_only_poisson_gp.stan"))

model_fits <- list()

for(i in 1:length(deer_species_all)) {

  if(evaluate_transects[i]) {
    model_to_fit <- model_negbin
  } else {
    model_to_fit <- model_negbin_co
  }

  model_fits[[i]] <- model_to_fit$sample(data = model_data[[i]],
                                         chains = nc,
                                 parallel_chains = nc,
                                 show_messages = TRUE,
                                 save_warmup = FALSE,
                                 iter_sampling = ni,
                                 iter_warmup = nw)

  model_fits[[i]]$save_object(paste0("outputs/fit_negbin_",deer_species_all[i],".rds"))

}

#### Generate Statewide Predictions ####
#### Make predictions ####
gp_preds_draws_all <- list()
pred_raster_means <- list()
for(i in 1:length(model_fits)) {

m_psi <- 1 + length(labels(terms(ab_formula))) # add intercept
n_site <- model_data[[i]]$n_site

# coordinates
coords_pred_scaled <- model_data[[i]]$coords_pred
coords <- model_data[[i]]$coords

obs_gp_function_posterior <-  model_fits[[i]]$draws('gp_predict',format = "df")[,1:n_site]  # gp covariance estimates
rho_posterior <- model_fits[[i]]$draws('rho', format = "df")$rho        # length scale
alpha_posterior <- model_fits[[i]]$draws('alpha', format = "df")$alpha  # marginal variance
beta_posterior<- model_fits[[i]]$draws('beta_psi', format = "df")[,1:m_psi] |> as.matrix()  # fixed effects
phi_posterior <- model_fits[[i]]$draws('phi', format = "df")$phi

nsims <- 50
iters <- sample(nrow(beta_posterior),
                size = nsims)

pred_raster_full <- terra::rast(prediction_raster)

Xpred_full <- model_data[[i]]$X_pred_psi

#Split predictions into chunks
chunks <- 10
pred_rows <- nrow(Xpred_full)

split_chunks <- function(x,n) split(1:x, cut(seq(1:x), n, labels = FALSE))

pred_chunks <- split_chunks(pred_rows, chunks)

gp_preds_draws <- list()

for(j in 1:chunks) {

  cat("Starting chunk", i, "at", as.character(Sys.time()))

  pred_grid <- pred_raster_full[pred_chunks[[j]],]

  # pred_grid # prediction grid including covariates
  OS<- rep(4, length(pred_chunks[[j]]))   # offset (cell size)
  Xpred <- Xpred_full[pred_chunks[[j]],]

  obs_distances <- calc_point_distances(coords)
  pred_distances <- calc_point_distances(coords_pred_scaled[pred_chunks[[j]],]) #, upper=TRUE, diag=TRUE)
  obs_to_pred_distances <- calc_point_distances(coords, coords_pred_scaled[pred_chunks[[j]],])

  gp_preds <- bettermc::mclapply(X = as.list(iters),
                                 FUN = predict_gp, obs_distances = obs_distances,
                                 pred_distances = pred_distances,
                                 obs_to_pred_distances = obs_to_pred_distances,
                                 obs_gp_function = obs_gp_function_posterior,
                                 alpha = alpha_posterior,
                                 rho_posterior,
                                 beta = beta_posterior,
                                 Xpred = Xpred,
                                 OS = OS,
                                 dist = "poisson",
                                 phi = phi_posterior,
                                 mc.progress = TRUE,
                                 mc.cores = 4,
                                 kernel = "quad")

  gp_preds_draws[[j]] <- do.call(cbind, gp_preds)
}

gp_preds_draws_all[[i]] <- dplyr::bind_rows(gp_preds_draws)
saveRDS(gp_preds_draws_all[[i]], paste0("outputs/gp_preds_draws_",deer_species_all[i],".rds"))

pred_raster_means[[i]] <- pred_raster_full[[1]]
terra::values(pred_raster_means[[i]])[!is.na(pred_raster_full[[1]])] <- apply(gp_preds_draws_all[[i]], 1, mean, na.rm = T)
terra::writeRaster(pred_raster_means[[i]], paste0("outputs/rasters/pred_mean_",deer_species_all[i],".tif"))

}



