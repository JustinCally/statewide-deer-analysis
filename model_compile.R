library(tidyverse)
library(dbplyr)
library(sf)
library(DBI)
library(kableExtra)
library(data.table)
library(Distance)
library(camtrapR)
library(cmdstanr)
library(activity)
library(bayesplot)
library(cowplot)
library(terra)
# remotes::install_github("JustinCally/weda")
library(weda)

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

# coordinates
coords <- st_as_sf(cams_curated, coords = c("Longitude", "Latitude"), crs = 4283) %>%
  st_transform(3111) %>%
  st_coordinates() %>%
  scale()

### Get the deer records (hog deer)
camtrap_records_deer <- tbl(con, dbplyr::in_schema(schema = "camtrap", table = "curated_camtrap_records")) %>%
  filter(scientific_name %in% deer_species & ProjectShortName %in% !!project_short_name) %>%
  dplyr::select(Species = scientific_name, SiteID, Distance = metadata_Distance, size = metadata_Multiples, Date, Time, Behaviour = metadata_Behaviour) %>%
  collect() %>%
  mutate(Time = as.POSIXct(Time, format = "%H:%M:%OS"),
         Time_n = as.numeric(Time, units = "secs"),
         Behaviour = na_if(x = Behaviour, y = "NA")) %>%
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


#' Availability calculation for a point transect
#'
#' @param delta bin width
#' @param midpoint bin midpoint
#' @param max_distance maximum distance of survey (truncation distance)
#'
#' @return numeric
availability_fn <- function(delta, midpoint, max_distance) {
  (2 * delta * midpoint)/max_distance^2
}

#' Probability closest animal is in a given bin, given a provided bin and number of individuals in a photo
#'
#' @param bin_start start point of bin
#' @param bin_end end point of bin
#' @param max_distance maximum distance of survey (truncation distance)
#' @param n number of individuals in photo
#'
#' @return numeric
pr_closest <- function(bin_start, bin_end, max_distance, n) {
  if(bin_start == 0) {
    prob_closer <- 0
  } else {
    closer_midpoint <- (bin_start)/2
    prob_closer <- 1 - (1 - availability_fn(bin_start, closer_midpoint, max_distance = max_distance))^n
  }
  if(bin_end == max_distance) {
    prob_further <- 0
  } else {
    further_midpoint <- (bin_end+max_distance)/2
    further_delta <- max_distance-bin_end
    prob_further <- availability_fn(further_delta, further_midpoint, max_distance = max_distance)
  }

  # Combined probability
  # probability that the closest indiv is in this bin is 1 minus the probability that the closest was in a closer bin, minus the probability that all individuals are in further bins
  pi_m <- 1 - (prob_closer + prob_further^n)

  return(pi_m)
}

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
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4283)

# read in spatial data rasters
processed_tifs <- list.files("/Volumes/DeerVic\ Photos/Processed_Rasters", full.names = TRUE)[stringr::str_detect(list.files("/Volumes/DeerVic\ Photos/Processed_Rasters"), ".tif$")]

processed_stack <- terra::rast(processed_tifs)

combined_spatial_data <- bind_cols(site_locs %>%
                                     dplyr::select(SiteID),
                                   terra::extract(x = processed_stack, terra::vect(site_locs %>%
                                                                                     sf::st_transform(3111)),
                                                  method = "bilinear", na.rm = T) %>%
                                     dplyr::select(-ID)) %>%
  # left_join(combined_spatial_data %>% dplyr::select(SiteID, WatercourseDistance) %>% st_drop_geometry()) %>%
  mutate(DeerHunt = as.factor(DeerHunt),
         SambarHunt = as.factor(DeerHunt)) %>%
  st_drop_geometry()


# Basic detection model
det_formula <- ~ scale(HerbaceousUnderstoryCover)
det_model_matrix <- model.matrix(det_formula, data = site_vars)

# Intercept only abundance model: change to informative predictors
ab_formula <- ~ scale(BIO12)  + scale(BIO04)  + scale(BIO01)  + scale(BIO12) + scale(sqrt(DistancetoWater)) + scale(BIO12*BIO04) + scale(BIO01*BIO12) + scale(BIO04*sqrt(DistancetoWater)) + scale(TreeDensity) + scale(sqrt(SLOPE)*BIO04) + scale(sqrt(PastureDistance)) + scale(sqrt(CrownGrazingDistance)) + scale(TWIND) + scale(BIO15)
ab_formula <- ~ scale(BIO12)  + scale(BIO04)  + scale(BIO01)  + scale(BIO12) + scale(BIO12*BIO04) + scale(TreeDensity) + scale(sqrt(PastureDistance)) + scale(TWIND) + scale(BIO15)
ab_model_matrix <- model.matrix(ab_formula, data = combined_spatial_data)

#### Prediction Data ####
vic_model_data_resampled <- rast("/Volumes/DeerVic\ Photos/MaxentStack/vic_model_data_resampled.tif")
vic_model_data_resampled_df <- as.data.frame(vic_model_data_resampled, xy = TRUE, na.rm = TRUE, cell = TRUE)
ab_model_pred_matrix <- model.matrix(ab_formula, data = vic_model_data_resampled_df)


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

# y_avail <- summarised_time %>%
#   group_by(SiteID, time_cut) %>%
#   summarise(n = sum(size, na.rm = T)) %>%
#   reshape2::dcast(SiteID ~ time_cut, value.var = "n") %>%
#   dplyr::select(-SiteID, -`NA`)
#
# y_avail[is.na(y_avail)] <- 0

# y1_df <- summarised_time %>%
#   filter(!is.na(time_midpoint)) %>%
#   group_by(time_midpoint) %>%
#   summarise(n = sum(size, na.rm = T))

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
              dplyr::transmute(SiteID, Survey = "Camera", Presence) %>%
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
               n_max = 20,
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
               coords = coords,
            y_pel = )

ni <- 400
nw <- 400
nt <- 1
nb <- 300
nc <- 6

inits = lapply(1:nc, function(i) list(beta_det=runif(2),
                                      beta_trans_det = runif(1),
                                      beta_psi_det = runif(1),
                                      beta_psi = runif(ncol(ab_model_matrix)),
                                      activ = runif(1),
                                      alpha = 1, rho=2,
                                      site_sd = runif(1),
                                      eps_ngs = rep(1/n_gs, n_gs),
                                      site_raw = rnorm(n_site),
                                      eta = rnorm(n_site)))

### Read in STAN model
# count only model
# Use relative file path
# count_only_model <- cmdstan_model(here::here("stan", "count_only.stan"))
# count + royle-nichols model
# modelRN <- cmdstan_model(here::here("_posts", "2022-06-29-sambards", "count_det_non_det_rn.stan"))

# Runs in 1-2 minutes
# fitcount <- count_only_model$sample(data = data,
#                                     chains = nc,
#                                     parallel_chains = nc,
#                                     show_messages = FALSE,
#                                     save_warmup = FALSE,
#                                     iter_sampling = ni,
#                                     iter_warmup = nw)
#
# fitcount$save_object("outputs/fitcount.rds")
#
# dens <- fitcount$summary("N_site")
#
#
# # Bind the density to the camera information
# density_at_sites <- cbind(cams_curated, dens) %>%
#   st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4286)
#
# # Currently looking at this map it seems like quite a few site coords are long (following up with Wildlife Unlimited)
# mapview::mapview(density_at_sites, zcol = "mean")
#### Integrated model (not used for hog deer) ####
model_hurdle <- cmdstan_model(here::here("stan", "count_det_non_det_rn_hurdle.stan"))
fit_hurdle <- model_hurdle$sample(data = data, chains = nc,
                             parallel_chains = nc,
                             show_messages = TRUE,
                             save_warmup = FALSE,
                             iter_sampling = ni,
                             iter_warmup = nw)


#### Integrated model (not used for hog deer) ####
modelRN <- cmdstan_model(here::here("stan", "count_det_non_det_rn.stan"))
fitintegrn <- modelRN$sample(data = data, chains = nc,
                             parallel_chains = nc,
                             show_messages = TRUE,
                             save_warmup = FALSE,
                             iter_sampling = ni,
                             iter_warmup = nw)

modelCO <- cmdstan_model(here::here("stan", "count_only_re.stan"))
fitCO <- modelCO$sample(data = data, chains = nc,
                             parallel_chains = nc,
                             show_messages = TRUE,
                             save_warmup = FALSE,
                             iter_sampling = ni,
                             iter_warmup = nw)

fitintegrn$save_object("outputs/fitintegrn.rds")

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
beta_psi_draws <- fitintegrn$draws("beta_psi", format = "matrix") %>% `colnames<-`(c("Intercept", labels(terms(ab_formula))))
mcmc_areas(beta_psi_draws)
#### Spatial predictions ####
rp <- fitintegrn$summary("pred", c("mean"))
spatial_preds <- bind_cols(vic_model_data_resampled_df, rp)
PredRast <- rast(vic_model_data_resampled, nlyr=1)

PredRast[spatial_preds$cell] <- spatial_preds[,"mean"]

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
