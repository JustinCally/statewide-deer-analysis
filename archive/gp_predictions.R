##
##----- prediction for GP process ---------------------------
library(tidyverse)
library(dbplyr)
library(sf)
library(DBI)
library(terra)
library(weda)
library(bettermc)
library(VicmapR)
source("functions/gp_functions.r")

# Ensure GoConnect is active and keyring is stored with password
con <- weda_connect(password = keyring::key_get(service = "ari-dev-weda-psql-01",
                                                username = "psql_user"), username = "psql_user")

### coords data ####
# Select deer species of interest
deer_species <- c("Cervus unicolor")
# Project of interest (2018 or 2023)
project_short_name <- c("hog_deer_2023", "StatewideDeer")
theta <- 40 * pi / 180 # camera angle in radians

cams_curated <- tbl(con, in_schema("camtrap", "curated_camtrap_operation")) %>%
  dplyr::filter(ProjectShortName %in% !!project_short_name) %>% # Only retrieve cam records for hog deer 2023
  collect() %>%
  mutate(DateTimeDeploy = as.POSIXct(DateTimeDeploy),
         DateTimeRetrieve = as.POSIXct(DateTimeRetrieve),
         Tk = as.numeric(DateTimeRetrieve - DateTimeDeploy, units = "secs"), #seconds
         Tk_prob = coalesce(as.numeric(as.POSIXct(Problem1_to,
                                                  format = "%Y-%m-%d %H:%M:%OS") -
                                         as.POSIXct(Problem1_from, format = "%Y-%m-%d %H:%M:%OS"),
                                       units = "secs"), 0),
         Tk_adj = Tk-Tk_prob,
         Tkt = Tk_adj / 2, # snapshot moments: every second second
         Effort = (Tk_adj*theta)/(2 * pi * 2)) %>%
  arrange(SiteID)

n_site <- nrow(cams_curated)

# coordinates
#### Prediction Data ####
vic_model_data_resampled <- rast("data/vic_model_data_resampled.tif")
# vic_model_data_resampled <- aggregate(vic_model_data_resampled, 10, na.rm = T)
# Add forest edges
woody_forest_edges <- rast("data/woody_forest_edges.tif") %>%
  `names<-`("ForestEdge") %>%
  terra::`crs<-`("epsg:3111")
woody_forest_edges_rp2 <- project(woody_forest_edges, vic_model_data_resampled)
vic_model_data_resampled_add <- c(vic_model_data_resampled, woody_forest_edges_rp2)

ab_formula <- ~ scale(BIO12)  + scale(BIO04)  + scale(BIO01) + scale(TreeDensity) + scale(sqrt(PastureDistance)) + scale(TWIND) + scale(BIO15) + scale(SLOPE) + scale(MRVBF) + scale(ForestEdge)

vic_model_data_resampled_df <- as.data.frame(vic_model_data_resampled_add, xy = TRUE, cell = TRUE) #%>%
ab_model_pred_matrix <- model.matrix(ab_formula, data = vic_model_data_resampled_df)
prop_pred <- rep(1, nrow(ab_model_pred_matrix))

m_psi <- ncol(ab_model_pred_matrix)

# coordinates
coords_pred <- data.frame(X = vic_model_data_resampled_df$x, Y = vic_model_data_resampled_df$y)

coords_pred_scaled <- coords_pred
coords_pred_scaled$X <- scales::rescale(coords_pred$X, to = c(0,1))
coords_pred_scaled$Y <- scales::rescale(coords_pred$Y, to = c(0,1))


coords <- st_as_sf(cams_curated, coords = c("Longitude", "Latitude"), crs = 4283) %>%
  st_transform(3111) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(X = scales::rescale(X, to = c(0,1), from = range(vic_model_data_resampled_df$x)),
         Y = scales::rescale(Y, to = c(0,1), from = range(vic_model_data_resampled_df$y)))

#crop to eastern_vic
east_vic <- vicmap_query("open-data-platform:ad_vicgov_region") %>%
  filter(vicgov_region %in% c("EASTERN METROPOLITAN")) %>%
  collect()  %>%
  st_transform(3111) %>%
  terra::vect()

pred_grid <- vic_model_data_resampled_add %>% terra::mask(east_vic)

pred_grid # prediction grid including covariates
OS<- as.numeric(na.omit(as.vector(terra::cellSize(pred_grid, mask = T))/1e6))    # offset (cell size)
Xpredall<- model.matrix(ab_formula, data = model.frame(ab_formula, vic_model_data_resampled_df, na.action=na.pass))  # Fixed effects
Xpred <- Xpredall[cells(pred_grid),]

library(ModelMatrixModel)
mm=ModelMatrixModel(ab_formula,combined_spatial_data_fix, center = T, scale = T)
mm_pred=predict(mm,vic_model_data_resampled_df %>% na.omit())

obs_distances <- calc_point_distances(coords)
pred_distances <- calc_point_distances(coords_pred_scaled[cells(pred_grid),]) #, upper=TRUE, diag=TRUE)
# pred_distances <- calc_point_distances(coords_pred_scaled)
obs_to_pred_distances <- calc_point_distances(coords, coords_pred_scaled[cells(pred_grid),])

# fitintegrn <- readRDS("outputs/fitintegrn.rds")

obs_gp_function_posterior <- fitintegrn$draws('gp_predict',format = "df")[,1:n_site]  # gp covariance estimates
rho_posterior <- fitintegrn$draws('rho', format = "df")$rho        # length scale
alpha_posterior <- fitintegrn$draws('alpha', format = "df")$alpha  # marginal variance
beta_posterior<- fitintegrn$draws('beta_psi', format = "df")[,1:m_psi] |> as.matrix()  # fixed effects
phi_posterior <- fitintegrn$draws('phi', format = "df")$phi

nsims<- 20
iters <- sample(nrow(beta_posterior), size = nsims)

gp_preds <- lapply(X = as.list(iters),
                     FUN = predict_gp, obs_distances = obs_distances,
                       pred_distances = pred_distances,
                       obs_to_pred_distances = obs_to_pred_distances,
                       obs_gp_function = obs_gp_function_posterior,
                       alpha = alpha_posterior,
                       rho_posterior,
                       beta = beta_posterior,
                       Xpred = Xpred,
                       OS = OS,
                       dist = "negbin",
                       phi = phi_posterior,
                       kernel = "exp")

gp_preds_draws<- do.call(cbind, gp_preds)
gp_preds_mean <- apply(gp_preds_draws, 1, mean, na.rm = T)

pred_grid_vals <- pred_grid[[1]]
values(pred_grid_vals)[cells(pred_grid)] <- gp_preds_mean
plot(pred_grid_vals)
# plan(sequential)



saveRDS(gp_preds, "outputs/gp_preds.rds")
