library(tidyverse)
library(dbplyr)
library(sf)
library(DBI)
library(kableExtra)
library(data.table)
# library(Distance)
library(camtrapR)
library(cmdstanr)
library(activity)
library(bayesplot)
library(cowplot)
library(terra)
# remotes::install_github("JustinCally/weda")
library(weda)
library(VicmapR)
source("functions/gp_functions.r")
source("functions/availability.R")

# STAN settings
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
options(mc.cores=4)

# Select deer species of interest
deer_species <- c("Cervus unicolor")
# Project of interest (2018 or 2023)
project_short_name <- c("hog_deer_2023", "StatewideDeer")
# Camera information
theta <- 40 * pi / 180 # camera angle in radians

# Connect to database
# Ensure GoConnect is active and keyring is stored with password
con <- weda_connect(password = keyring::key_get(service = "ari-dev-weda-psql-01",
                                                username = "psql_user"), username = "psql_user")

# View what camera trap data is available (query the project)
# use the curated data
camtrap_projects_on_database <- tbl(con, in_schema("camtrap", "curated_project_information")) %>%
  collect()

# Can look at a data dictionary of all tables in the 'camtrap schema' with weda
# View(weda::data_dictionary)


### Get camera information for the 64 camera sites
### For cameras that have problems you need to use the date from when the problem started
cams_curated <- tbl(con, in_schema("camtrap", "curated_camtrap_operation")) %>%
  filter(ProjectShortName %in% !!project_short_name) %>% # Only retrieve cam records for hog deer 2023
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



### Get the deer records (hog deer)
camtrap_records_deer <- tbl(con, dbplyr::in_schema(schema = "camtrap", table = "curated_camtrap_records")) %>%
  filter(scientific_name %in% deer_species & ProjectShortName %in% !!project_short_name) %>%
  dplyr::select(Species = scientific_name, SiteID, Distance = metadata_Distance, size = metadata_Multiples, Date, Time, Behaviour = metadata_Behaviour) %>%
  collect() %>%
  mutate(Distance = case_when(Distance == "NA" | is.na(Distance) ~ "999",
                              TRUE ~ Distance)) %>%
  mutate(Time = as.POSIXct(Time, format = "%H:%M:%OS"),
         Time_n = as.numeric(Time, units = "secs"),
         Behaviour = na_if(x = Behaviour, y = "NA")) %>%
  rowwise() %>%
  mutate(DistanceMod = list(stringr::str_split(Distance, pattern = "_&_")[[1]])) %>%
  mutate(Distance = DistanceMod[which.min(as.numeric(stringr::str_extract(DistanceMod, pattern = "[0-9]+")))]) %>%
  mutate(Distance = case_when(Distance == "0-2.5" ~ "0 - 2.5",
                              Distance == "999" ~ NA_character_,
                              TRUE ~ Distance)) %>%
  ungroup() %>%
  filter(Time_n %% 2 == 0 & #& #snapshot moment interval of 2s
        is.na(Behaviour)) %>% # filter out behaviors such as camera or marker interaction
  group_by(SiteID, Time_n) %>%
  slice(1) %>% # if two photos occur in one second take only one (snapshot moment = 2)
  ungroup()

### Tidy format for the records
dcount<- camtrap_records_deer %>%
  dplyr::select(Species, SiteID, Distance, size, Date, Time, Behaviour) %>%
  full_join(cams_curated %>%
              dplyr::select(SiteID, Tkt, Effort), by="SiteID") %>%
  mutate(object=1:nrow(.)) %>%
  mutate(size = if_else(is.na(size),0L, size)) %>%
  arrange(SiteID)




### Format data in a time format for availability
summarised_count <- dcount
summarised_count$Distance <- factor(summarised_count$Distance, levels = c("0 - 2.5", "2.5 - 5", "5 - 7.5", "7.5 - 10", "10+"))

summarised_time <- dcount %>%
  mutate(rtime = as.integer(format(Time, "%H")),
         time_cut = cut(rtime, breaks =  seq(from = 0, to = 24, by = 2)),
         time_midpoint = 2*as.integer(time_cut)-1)

### Get herbaceous understorey for detection submodel
### This site data is in the deervic schema (for statewide project still in development)
site_vars <- tbl(con, in_schema("deervic", "curated_site_data")) %>%
  filter(SiteID %in% !!c(cams_curated$SiteID, "47191")) %>%
  collect() %>%
  mutate(HerbaceousUnderstoryCover = NNWHUCover + ENWHUCover,
         SiteID = case_when(SiteID == "47191" & CameraID == "HO04101053" ~ "47191A",
                                   SiteID == "47191" & CameraID != "HO04101053" ~ "47191B",
                                   TRUE ~ SiteID)) %>% # native + exotic herbaceous cover
  arrange(SiteID) %>%
  as.data.frame()

# abundance site variables #
site_locs <- site_vars %>%
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4283) %>%
  st_transform(3111)

# read in spatial data rasters
processed_tifs <- list.files("/Volumes/DeerVic\ Photos/Processed_Rasters", full.names = TRUE)[stringr::str_detect(list.files("/Volumes/DeerVic\ Photos/Processed_Rasters"), ".tif$")]

processed_stack <- terra::rast(processed_tifs)

# slga_stack <- rast("/Volumes/Cally_Camtr/StatewideRasters/slga_stack.tif")

# Add forest edges
woody_forest_edges <- rast("data/woody_forest_edges.tif") %>%
  `names<-`("ForestEdge") %>%
  terra::`crs<-`("epsg:3111")

woody_forest_edges_rp <- project(woody_forest_edges, processed_stack)

writeRaster(woody_forest_edges_rp, "/Volumes/Cally_Camtr/StatewideRasters/Processed_Rasters/woody_edges_1km.tif")

processed_stack_add <- c(processed_stack, woody_forest_edges_rp)

combined_spatial_data <- bind_cols(site_locs %>%
                                     dplyr::select(SiteID),
                                   terra::extract(x = processed_stack_add, terra::vect(site_locs),
                                                    method = "bilinear", na.rm = T, fun = "mean") %>%
                                     dplyr::select(-ID)) %>%
  st_drop_geometry()

# For some sites use a more expansive buffer and replace NAs
combined_spatial_data_buffer <- bind_cols(site_locs %>%
                                     dplyr::select(SiteID),
                                   terra::extract(x = processed_stack_add, terra::vect(site_locs %>% st_buffer(1000)),
                                                  method = "bilinear", na.rm = T, fun = "mean") %>%
                                     dplyr::select(-ID)) %>%
  st_drop_geometry()

combined_spatial_data_fix <- combined_spatial_data

combined_spatial_data_fix[is.na(combined_spatial_data)] <- combined_spatial_data_buffer[is.na(combined_spatial_data)]
combined_spatial_data_fix[is.na(combined_spatial_data_fix)] <- 0

combined_spatial_data_fix <- mutate(combined_spatial_data_fix, across(c(everything(), -SiteID), .fns = as.numeric))

# Basic detection model
det_formula <- ~ scale(HerbaceousUnderstoryCover)
det_model_matrix <- model.matrix(det_formula, data = site_vars)

# Intercept only abundance model: change to informative predictors
ab_formula <- ~ scale(BIO12)  + scale(BIO04)  + scale(BIO01) + scale(TreeDensity) + scale(sqrt(PastureDistance)) + scale(TWIND) + scale(BIO15) + scale(SLOPE) + scale(MRVBF) + scale(sqrt(ForestEdge))
# ab_model_matrix <- model.matrix(ab_formula, data = combined_spatial_data_fix %>% mutate(BIOREGION = factor(round(BIOREGION))))

#### Prediction Data ####
vic_model_data_resampled <- rast("/Volumes/DeerVic\ Photos/MaxentStack/vic_model_data_resampled.tif")
woody_forest_edges_rp2 <- project(woody_forest_edges, vic_model_data_resampled)
vic_model_data_resampled_add <- c(vic_model_data_resampled, woody_forest_edges_rp2) %>%
  aggregate(2, na.rm = TRUE)

site_loc_cells <- cells(vic_model_data_resampled_add, vect(site_locs))[,"cell"]

vic_model_data_resampled_df <- as.data.frame(vic_model_data_resampled_add, xy = TRUE, cell = TRUE) #%>%

ab_model_pred_matrix_full <- model.matrix(ab_formula, data = bind_rows(combined_spatial_data_fix, vic_model_data_resampled_df))
ab_model_pred_matrix <- ab_model_pred_matrix_full[(nrow(combined_spatial_data_fix)+1):nrow(ab_model_pred_matrix_full),]
ab_model_matrix <- ab_model_pred_matrix_full[1:nrow(combined_spatial_data_fix),]

prop_pred <- rep(1, nrow(ab_model_pred_matrix))

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

n_site <- length(unique(dcount$SiteID))
n_distance_bins <- 5
delta <- 2.5
midpts <- c(1.25, 3.75, 6.25, 8.75, 11.25)
max_distance <- 12.5
theta_frac <- 42/360
Tkt <- dcount %>% group_by(SiteID) %>% slice(1) %>% pull(Tkt)

# Get number of observations at each group size
n_obs <- table(dcount$SiteID, dcount$size)[,-1] # remove zeros
n_gs <- ncol(n_obs)
gs <- as.integer(colnames(n_obs))

yl <- list()

for(i in 1:n_gs) {
  yl[[i]] <- summarised_count %>%
    mutate(SiteID = factor(SiteID, levels = unique(summarised_count$SiteID))) %>%
    group_by(SiteID, Distance, .drop = FALSE) %>%
    filter(size == gs[i]) %>%
    summarise(n = length(size)) %>%
    reshape2::dcast(SiteID ~ Distance, value.var = "n") %>%
    dplyr::select(-SiteID, -`NA`)

  yl[[i]][is.na(yl[[i]])] <- 0
}

y <- str2str::ld2a(yl)


### Adjusted availability function

pa <- matrix(data = NA, nrow = max(gs), ncol = length(midpts))

for(j in 1:max(gs)) {

  pa[j, ] <- sapply(midpts, function(x) {pr_closest(x-(0.5*delta),
                                                    x+(0.5*delta),
                                                    max_distance = max_distance,
                                                    n = j)})


}


#### Next section only useful if you are using transects/multiple observations ####
# NOT USED FOR HOG DEER: only 1 site did transects #
### get any observed across survey types ###
# Load in transects that have been assifgned to species
# source("functions/species_assign.R")
# source("species_assign_wrangle.R")
Rusaunicolor_Detection <- readRDS("data/Rusaunicolor_Detection.rds")

presence_absence <-  tbl(con, in_schema("camtrap", "processed_site_substation_presence_absence")) %>%
              filter(scientific_name %in% deer_species & ProjectShortName %in% !!project_short_name) %>%
              dplyr::transmute(SiteID, Survey = "Camera", Presence, Count = 1) %>%
  collect()

transects <- Rusaunicolor_Detection %>%
  bind_rows(presence_absence) %>%
  arrange(SiteID, Survey) %>%
  mutate(surveyed = 1,
         row = 1:n()) %>%
  as.data.frame()

transects_ne <- site_vars %>%
  left_join(transects) %>%
  mutate(surveyed = coalesce(surveyed, 0),
         site = as.numeric(factor(SiteID)))

y_pel <- transects_ne %>%
  # filter(Survey == "Pellet") %>%
  group_by(SiteID) %>%
  summarise(Count = sum(Count, na.rm = T)) %>%
  pull(Count)

transects_ne_f <- transects_ne %>%
  filter(surveyed > 0)

n_survey <- transects_ne %>%
  group_by(SiteID) %>%
  summarise(n = sum(surveyed)) %>%
  pull(n)

# get start and end indices to extract slices of y for each site
start_idx <- rep(0, n_site)
end_idx <- rep(0, n_site)
for (i in 1:n_site) {
  if (n_survey[i] > 0) {
    site_indices <- which(transects_ne_f$site == i)
    start_idx[i] <- site_indices[1]
    end_idx[i] <- site_indices[n_survey[i]]
  }
}

transect_formula <- ~Survey

# Make sure camera is the intercept
transect_mm <- model.matrix(object = transect_formula, data = transects)

# If include cameras in the transect this bit is duplicated
cam_max <- apply(n_obs, MARGIN = 1, FUN = function(x) {
  m <- max(x)
  if(m > 0) {
    m <- 1
  }
  return(m)
  })

any_seen <- rep(0, n_site)
for (i in 1:n_site) {
  if (n_survey[i] > 0) {
    any_seen[i] <- max(transects_ne$Presence[start_idx[i]:end_idx[i]], na.rm = T)
  }
  any_seen[i] <- max(c(any_seen[i], cam_max[i]))
}




#### Determine availability prior ####
trigger.events <- camtrap_records_deer %>%
  # filter(is.na(Behaviour)) %>% ### Remove biased behaviour
  mutate(Sample.Label = as.character(SiteID))

trigger.events$rtime <- gettime(x = trigger.events$Time, tryFormats = "%Y-%m-%d %H:%M:%OS")

act_result <- fitact(trigger.events$rtime, sample="model", reps=50, bw = 30)

avail <- list(creation=data.frame(rate = act_result@act[1],
                                  SE   = act_result@act[2]))

# View activity and avialbility
plot(act_result)
avail # looks to be a bit higher than other deer and more variable due to lower counts

get.beta.param<- function(mu,sdev){
  # get parameters of beta by method of moments
  a<- mu*((mu*(1-mu))/sdev^2 - 1)
  b<- (1-mu)*((mu*(1-mu))/sdev^2 - 1)
  if(is.na(a) || is.na(b)) stop("Invalid Beta parameters")
  else return(list(a=a,b=b))
}

beta.pars<- get.beta.param(avail$creation$rate, 2*avail$creation$SE)

### Prepare data for STAN model

# Restrict n_max
n_max <- cam_max
n_max[n_max == 0] <- 3
n_max[n_max == 1] <- 20

data = list(N=sum(dcount$size, na.rm = T),
               delta = 2.5,
               n_site = n_site,
               n_gs = n_gs,
               gs = gs,
               n_distance_bins = n_distance_bins,
               midpts = midpts,
               max_distance = max_distance,
               max_int_dist = as.integer(round(max_distance)),
               theta_frac = theta_frac,
               effort = Tkt,
               n_obs = n_obs,
               y = y,
               pa = pa, # availability varies across sites based on number of individuals
               det_model_matrix = det_model_matrix,
               det_ncb = ncol(det_model_matrix),
               bshape = beta.pars[[1]],
               bscale = beta.pars[[2]],
               n_max = n_max,
               # transect based data: Can provide it bu not used in hog deer
               trans = nrow(transects),
               y2 = transects$Presence,
               start_idx = start_idx,
               end_idx = end_idx,
               trans_det_ncb = ncol(transect_mm),
               trans_pred_matrix = transect_mm,
               any_seen = any_seen,
               n_survey = n_survey,
               m_psi =  ncol(ab_model_matrix),
               X_psi = ab_model_matrix,
               npc = nrow(ab_model_pred_matrix),
               X_pred_psi = ab_model_pred_matrix,
               prop_pred = prop_pred,
               coords = coords,
               # y_pel = y_pel,
               reciprocal_phi_scale = 1,
               hs_df = 1,
               hs_df_global = 1,
               hs_scale_global = 2/sqrt(n_site), # ratio of expected non-zero to zero divided by total observation as per brms convention
               hs_scale_slab = 1,
               hs_df_slab = 4)

ni <- 250
nw <- 250
nt <- 1
nb <- 300
nc <- 6

inits = lapply(1:nc, function(i) list(beta_det=runif(2),
                                      beta_trans_det = runif(1),
                                      beta_psi_det = runif(1),
                                      beta_psi = runif(ncol(ab_model_matrix)),
                                      activ = runif(1),
                                      alpha = 2, rho=1,
                                      phi = rep(1, n_site),
                                      site_sd = runif(1),
                                      eps_ngs = rep(1/n_gs, n_gs),
                                      site_raw = rnorm(n_site),
                                      eta = rnorm(n_site)))


# mimimum distance
coords_m <-  st_as_sf(cams_curated, coords = c("Longitude", "Latitude"), crs = 4283) %>%
  st_transform(3111)

st_dist <- st_distance(coords_m) %>% units::set_units("km")
st_dist <- st_dist[as.numeric(st_dist) > 0]
st_dist %>% min()
st_dist %>% quantile(probs = seq(0, 1, 0.1))

distance_scaled <- dist(coords) # was messing up

distance_scaled %>% quantile(probs = seq(0, 1, 0.1))

# Upper limit will be ~ 60th quantile
# Lower limit will be ~ 30th quantile

# prior_tune <- cmdstan_model(here::here("stan", "gp_prior_tune.stan"))

# priorfit<- prior_tune$sample(iter_sampling=1, iter_warmup=0, chains=1,
#                                seed=5838298, fixed_param = TRUE)


model_negbin <- cmdstan_model(here::here("stan", "count_det_nondet_negbin_gp.stan"))
fitintegrn<- model_negbin$sample(data = data, chains = nc,
                                  parallel_chains = nc,
                                  show_messages = TRUE,
                                  save_warmup = FALSE,
                                  iter_sampling = ni,
                                  iter_warmup = nw)
Sys.time()

fitintegrn$save_object("outputs/fitintegrn.rds")

Sys.time()

#### Make predictions ####
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


# fitintegrn <- readRDS("outputs/fitintegrn.rds")

obs_gp_function_posterior <- fitintegrn$draws('gp_predict',format = "df")[,1:n_site]  # gp covariance estimates
rho_posterior <- fitintegrn$draws('rho', format = "df")$rho        # length scale
alpha_posterior <- fitintegrn$draws('alpha', format = "df")$alpha  # marginal variance
beta_posterior<- fitintegrn$draws('beta_psi', format = "df")[,1:m_psi] |> as.matrix()  # fixed effects
phi_posterior <- fitintegrn$draws('phi', format = "df")$phi

nsims<- 12
iters <- sample(nrow(beta_posterior), size = nsims)

# Regions
vic_regions <- vicmap_query("open-data-platform:delwp_region") %>%
  collect()  %>%
  st_transform(3111)

ab_model_pred_matrix_full_NAs <- model.matrix(ab_formula,
                                              data = model.frame(ab_formula, bind_rows(combined_spatial_data_fix, vic_model_data_resampled_df), na.action = na.pass))

gp_preds_draws <- list()

for(i in 1:nrow(vic_regions)) {

  cat("Starting region", i, "at", as.character(Sys.time()))

pred_grid <- vic_model_data_resampled_add %>% terra::mask(vect(vic_regions[i,]))

# pred_grid # prediction grid including covariates
OS<- rep(4, length(cells(pred_grid)))   # offset (cell size)
Xpred <- ab_model_pred_matrix_full_NAs[(n_site + cells(pred_grid)),]

obs_distances <- calc_point_distances(coords)
pred_distances <- calc_point_distances(coords_pred_scaled[cells(pred_grid),]) #, upper=TRUE, diag=TRUE)
obs_to_pred_distances <- calc_point_distances(coords, coords_pred_scaled[cells(pred_grid),])

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
                   dist = "negbin",
                   phi = phi_posterior, mc.progress = TRUE, mc.cores = 4,
                   kernel = "quad")

gp_preds_draws[[i]] <- do.call(cbind, gp_preds)
}

names(gp_preds_draws) <- vic_regions$delwp_region

saveRDS(gp_preds_draws, "outputs/gp_preds_draws.rds")

Sys.time()


# Make Raster of each region (mean)
pred_grid_vals <- list()

for(i in 1:nrow(vic_regions)) {

gp_preds_mean <- apply(gp_preds_draws[[i]], 1, mean, na.rm = T)
pred_grid <- vic_model_data_resampled_add %>% terra::mask(vect(vic_regions[i,]))
pred_grid_vals[[i]] <- pred_grid[[1]]
values(pred_grid_vals[[i]]) <- NA
values(pred_grid_vals[[i]])[cells(pred_grid)] <- gp_preds_mean
names(pred_grid_vals[[i]]) <- "Abundance"
}

statewide_pred_rast <- mosaic(sprc(pred_grid_vals))
plot(statewide_pred_rast)
# inspect the length scale
ar <- fitintegrn$summary(c("rho", "alpha"))
gp_draws <- fitintegrn$draws(c("rho", "alpha"), format = "matrix")
mcmc_areas(gp_draws)

gp_predict <- fitintegrn$summary("gp_predict")

rn_dens <- fitintegrn$summary("N_site")
# Bind the density to the camera information
density_at_sites_rn <- cbind(cams_curated, rn_dens) %>%
  mutate(means_sqrt = sqrt(mean)) %>%
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4283)

mapview::mapview(density_at_sites_rn, zcol = "mean")
# rn_lamb <- fitintegrn$summary("Site_lambda")
# rn_r <- fitintegrn$summary("r")
# rn_psi <- fitintegrn$summary("psi")
fitintegrn$summary("Nhat")
survey_det <- fitintegrn$summary("beta_trans_det")
rn_eps_site <- fitintegrn$summary("eps_site")
rn_beta_psi <- fitintegrn$summary("beta_psi")
beta_psi_draws <- fitintegrn$draws("beta_psi", format = "matrix") %>% `colnames<-`(colnames(ab_model_matrix))
mcmc_areas(beta_psi_draws)

#### Spatial predictions ####
rp <- fitintegrn$summary("pred", c("mean"))
spatial_preds <- bind_cols(vic_model_data_resampled_df, rp)
PredRast <- rast(vic_model_data_resampled, nlyr=1)

PredRast[spatial_preds$cell] <- spatial_preds[,"mean"]
# PredRast[!(PredRast %in% spatial_preds$cell)] <- NA
plot(PredRast)

#### Spatial interpolation from points ####
library(gstat)
gs <- gstat(formula=mean~1, locations=~x+y, data=cbind(density_at_sites_rn %>%
                                                         st_transform(3111) %>%
                                                         st_coordinates(), rn_dens) %>%
              rename(x = X, y = Y), nmax=5, set=list(idp = 0))
nn <- interpolate(PredRast, gs, debug.level=0)
nnmsk <- mask(nn, PredRast)
sum(values(nnmsk[[1]], na.rm = T))
plot(nnmsk[[1]])

#### Detection Function ####
det_curve <- fitintegrn$draws("DetCurve", format = "draws_matrix") %>%
  as.data.frame() %>%
  head(250) %>%
  pivot_longer(cols = everything())

det_curve_wr <- det_curve %>%
  mutate(var = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) %>%
  separate(var, into = c("s", "Distance"), sep = ",")

det_vars_pred <- site_vars %>%
  mutate(s = as.character(1:nrow(.)),
         herbaceouslvl = cut(HerbaceousUnderstoryCover,
                             breaks = c(0, 25, 50, 75, 100), include.lowest = T, right = FALSE))

det_curve_sum <- det_curve_wr %>%
  mutate(Distance = as.numeric(Distance)-1) %>%
  left_join(det_vars_pred) %>%
  group_by(herbaceouslvl, Distance) %>%
  summarise(median = quantile(value, 0.5),
            q5 = quantile(value, 0.05),
            q95 = quantile(value, 0.95))


y_combined <- colSums(y[,,1]) %>% # just for group size 1
  as.data.frame() %>%
  `colnames<-`("Count") %>%
  mutate(Distance = midpts,
         Prop = Count/(max(Count)),
         CountS = Count/(2 * 2.5 * Distance)/max_distance^2,
         PropS = CountS/(max(CountS)))

detection_plot_HN <- ggplot(aes(x = Distance), data = det_curve_sum) +
  geom_col(aes(y = PropS), fill = "grey70", colour = "grey20", width = 2.5, data = y_combined, alpha = 0.7) +
  # geom_errorbar(aes(ymin = PropSq5, ymax = PropSq95),  data = y_combined) +
  # geom_line(aes(y = HNS)) +
  geom_ribbon(aes(ymin = q5, ymax = q95, fill = herbaceouslvl), alpha = 0.5) +
  geom_line(aes(y = median, colour = herbaceouslvl)) +
  ggtitle("MRDS Model", subtitle = "Red = p(dist)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  ylab("Detection probability (p)") +
  theme_classic()

detection_plot_HN
