# Prepare Data Function

prepare_model_data_multispecies <- function(species,
                               projects,
                               buffer,
                               detection_formula,
                               abundance_formula,
                               transect_formula,
                               con,
                               raster_dir,
                               prediction_raster,
                               n_max_no_det,
                               n_max_det,
                               evaltransects = TRUE,
                               snapshot_interval = 2) {

  if(species == "All deer") {
    species <- c("Cervus unicolor", "Dama dama", "Cervus elaphus", "Axis porcinus")
  }

  # Camera information
  theta <- 40 * pi / 180 # camera angle in radians

  ### Get camera information for the camera sites
  ### For cameras that have problems you need to use the date from when the problem started
  cams_curated <- dplyr::tbl(con, dbplyr::in_schema("camtrap", "curated_camtrap_operation")) %>%
    dplyr::filter(ProjectShortName %in% !!projects) %>% # Only retrieve cam records for hog deer 2023
    dplyr::collect() %>%
    dplyr::mutate(DateTimeDeploy = as.POSIXct(DateTimeDeploy),
           DateTimeRetrieve = as.POSIXct(DateTimeRetrieve),
           Tk = as.numeric(DateTimeRetrieve - DateTimeDeploy, units = "secs"), #seconds
           Tk_prob = dplyr::coalesce(as.numeric(as.POSIXct(Problem1_to,
                                                    format = "%Y-%m-%d %H:%M:%OS") -
                                           as.POSIXct(Problem1_from, format = "%Y-%m-%d %H:%M:%OS"),
                                         units = "secs"), 0),
           Tk_adj = Tk-Tk_prob,
           Tkt = Tk_adj / snapshot_interval, # snapshot moments: every second second
           Effort = (Tk_adj*theta)/(snapshot_interval * pi * 2)) %>%
    dplyr::arrange(SiteID)


  ### Get the deer records (hog deer)
  camtrap_records_deer <- dplyr::tbl(con, dbplyr::in_schema(schema = "camtrap", table = "curated_camtrap_records")) %>%
    dplyr::filter(scientific_name %in% species & ProjectShortName %in% !!projects) %>%
    dplyr::select(Species = scientific_name, SiteID, Distance = metadata_Distance, size = metadata_Multiples, Date, Time, Behaviour = metadata_Behaviour) %>%
    dplyr::collect() %>%
    dplyr::mutate(Distance = dplyr::case_when(Distance == "NA" | is.na(Distance) ~ "999",
                                TRUE ~ Distance)) %>%
    dplyr::mutate(Time = as.POSIXct(Time, format = "%H:%M:%OS"),
           Time_n = as.numeric(Time, units = "secs"),
           Behaviour = na_if(x = Behaviour, y = "NA")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(DistanceMod = list(stringr::str_split(Distance, pattern = "_&_")[[1]])) %>%
    dplyr::mutate(Distance = DistanceMod[which.min(as.numeric(stringr::str_extract(DistanceMod, pattern = "[0-9]+")))]) %>%
    dplyr::mutate(Distance = dplyr::case_when(Distance == "0-2.5" ~ "0 - 2.5",
                                Distance == "999" ~ NA_character_,
                                TRUE ~ Distance)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Time_n %% snapshot_interval == 0 & #& #snapshot moment interval of 2s
             is.na(Behaviour)) %>% # filter out behaviors such as camera or marker interaction
    dplyr::group_by(SiteID, Time_n, Species) %>%
    dplyr::slice(1) %>% # if two photos occur in one second take only one (snapshot moment = 2)
    dplyr::ungroup()

  ### Tidy format for the records
  dcount<- camtrap_records_deer %>%
    dplyr::select(Species, SiteID, Distance, size, Date, Time, Behaviour) %>%
    dplyr::full_join(cams_curated %>%
                dplyr::select(SiteID, Tkt, Effort), by="SiteID") %>%
    dplyr::mutate(object=1:nrow(.)) %>%
    dplyr::mutate(size = if_else(is.na(size),0L, size)) %>%
    dplyr::arrange(SiteID)

  ### Format data in a time format for availability
  summarised_count <- dcount
  summarised_count$Distance <- factor(summarised_count$Distance, levels = c("0 - 2.5", "2.5 - 5", "5 - 7.5", "7.5 - 10", "10+"))

  summarised_time <- dcount %>%
    dplyr::mutate(rtime = as.integer(format(Time, "%H")),
           time_cut = cut(rtime, breaks =  seq(from = 0, to = 24, by = 2)),
           time_midpoint = 2*as.integer(time_cut)-1)

  site_vars <- dplyr::tbl(con, dbplyr::in_schema("deervic", "curated_site_data")) %>%
    dplyr::filter(SiteID %in% !!c(cams_curated$SiteID, "47191")) %>%
    dplyr::collect() %>%
    dplyr::mutate(HerbaceousUnderstoryCover = NNWHUCover + ENWHUCover,
           SiteID = dplyr::case_when(SiteID == "47191" & CameraID == "HO04101053" ~ "47191A",
                              SiteID == "47191" & CameraID != "HO04101053" ~ "47191B",
                              TRUE ~ SiteID)) %>% # native + exotic herbaceous cover
    dplyr::arrange(SiteID) %>%
    as.data.frame()

  # abundance site variables #
  site_locs <- site_vars %>%
    sf::st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4283) %>%
    sf::st_transform(3111)

  # Read in rasters
  processed_tifs <- list.files(raster_dir, full.names = TRUE)[stringr::str_detect(list.files(raster_dir), ".tif$")]

  processed_stack <- terra::rast(processed_tifs)

  # For some sites use a more expansive buffer and replace NAs
  combined_spatial_data <- dplyr::bind_cols(site_locs %>%
                                       dplyr::select(SiteID),
                                            exactextractr::exact_extract(x = processed_stack,
                                                                         site_locs %>%
                                                                           sf::st_buffer(buffer),
                                                                         fun = "mean",
                                                                         colname_fun = function(values, fun_name, weights, fun_value, nvalues, nweights)  values)) %>%
    sf::st_drop_geometry()

  # Elevation is an issue
  combined_spatial_data$Elevation[is.na(combined_spatial_data$Elevation)] <- 0

  # No categorical predictors
  combined_spatial_data_fix <- dplyr::mutate(combined_spatial_data,
                                             dplyr::across(c(dplyr::everything(), -SiteID), .fns = as.numeric))

  # Basic detection model
  det_model_matrix <- model.matrix(detection_formula, data = site_vars)

  #### Prediction Data ####
  vic_model_data_resampled <- terra::rast(prediction_raster)

  site_loc_cells <- terra::cells(vic_model_data_resampled, terra::vect(site_locs))[,"cell"]

  vic_model_data_resampled_df <- terra::as.data.frame(vic_model_data_resampled, xy = TRUE, cell = TRUE)

  ab_model_pred_matrix_full <- model.matrix(abundance_formula,
                                            data = bind_rows(combined_spatial_data_fix, vic_model_data_resampled_df))

  ab_model_pred_matrix <- ab_model_pred_matrix_full[(nrow(combined_spatial_data_fix)+1):nrow(ab_model_pred_matrix_full),]
  ab_model_matrix <- ab_model_pred_matrix_full[1:nrow(combined_spatial_data_fix),]
  prop_pred <- rep(1, nrow(ab_model_pred_matrix))

  #### Coordinate data ####
  # coordinates
  coords_pred <- data.frame(X = vic_model_data_resampled_df$x,
                            Y = vic_model_data_resampled_df$y)

  coords_pred_scaled <- coords_pred
  coords_pred_scaled$X <- scales::rescale(coords_pred$X, to = c(0,1))
  coords_pred_scaled$Y <- scales::rescale(coords_pred$Y, to = c(0,1))


  coords <- sf::st_as_sf(cams_curated,
                         coords = c("Longitude", "Latitude"), crs = 4283) %>%
    sf::st_transform(3111) %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::mutate(X = scales::rescale(X, to = c(0,1),
                                      from = range(vic_model_data_resampled_df$x)),
                  Y = scales::rescale(Y, to = c(0,1),
                                      from = range(vic_model_data_resampled_df$y)))

  #### Model data boundries ####
  n_site <- length(unique(dcount$SiteID))
  n_distance_bins <- 5
  delta <- 2.5
  midpts <- c(1.25, 3.75, 6.25, 8.75, 11.25)
  max_distance <- 12.5
  theta_frac <- 42/360
  Tkt <- dcount %>%
    dplyr::group_by(SiteID) %>%
    dplyr::slice(1) %>%
    dplyr::pull(Tkt)

  ### Variable group size ####
  # Get number of observations at each group size
  n_obs <- table(dcount$SiteID, dcount$size, dcount$Species)[,-1,] # remove zeros
  n_gs <- ncol(n_obs)
  gs <- as.integer(colnames(n_obs))

  yl <- list()

  for(i in 1:n_gs) {
    yl[[i]] <- summarised_count %>%
      dplyr::mutate(SiteID = factor(SiteID,
                                    levels = unique(summarised_count$SiteID))) %>%
      dplyr::group_by(SiteID, Distance, .drop = FALSE) %>%
      dplyr::filter(size == gs[i]) %>%
      dplyr::summarise(n = length(size)) %>%
      reshape2::dcast(SiteID ~ Distance, value.var = "n") %>%
      dplyr::select(-SiteID, -`NA`)

    yl[[i]][is.na(yl[[i]])] <- 0
  }
  # observations per group
  y <- str2str::ld2a(yl)

  #### Adjusted availability function for group sizes ####
  pa <- matrix(data = NA, nrow = max(gs), ncol = length(midpts))
  for(j in 1:max(gs)) {
    pa[j, ] <- sapply(midpts, function(x) {pr_closest(x-(0.5*delta),
                                                      x+(0.5*delta),
                                                      max_distance = max_distance,
                                                      n = j)})
  }

  if(evaltransects) {
  #### Transect-based detections ####
  #### Next section only useful if you are using transects/multiple observations ####
  # NOT USED FOR HOG DEER: only 1 site did transects #
  ### get any observed across survey types ###
  # Load in transects that have been assifgned to species
  # source("functions/species_assign.R")
  # source("species_assign_wrangle.R")

  filepath_det <- sapply(species, function(x) paste0("data/transects/",stringr::str_replace_all(x, " ", "_"), "_Detection.rds"))
  filepath_det<- filepath_det[sapply(filepath_det, file.exists)]
  Deer_Detection <- lapply(filepath_det, readRDS)
  for(i in 1:length(species)) {
    Deer_Detection[[i]]$Species <- species[i]
  }
  Deer_Detection <- dplyr::bind_rows(Deer_Detection) %>%
    dplyr::group_by(Species, SiteID, Distance, Survey) %>%
    dplyr::mutate(Transect = 1:n()) %>%
    dplyr::group_by(Species, SiteID, Distance, Survey, Transect) %>%
    dplyr::summarise(Count = max(Count, na.rm = T),
              Presence = max(Presence, na.rm = T)) %>%
    dplyr::ungroup()

  # filter deer detections to only include methods that detected deer
  # methods_of_det <- Deer_Detection %>%
  #   dplyr::filter(Presence == 1) %>%
  #   dplyr::pull(Survey) %>%
  #   unique()
  #
  # Deer_Detection <- Deer_Detection %>%
  #   dplyr::filter(Survey %in% methods_of_det)

  pa_empty <- expand.grid(Species = species,
                          SiteID = combined_spatial_data$SiteID,
                          Survey = "Camera")

  presence_absence_ic <-  dplyr::tbl(con,
                                  dbplyr::in_schema("camtrap", "processed_site_substation_presence_absence")) %>%
    dplyr::filter(scientific_name %in% species & ProjectShortName %in% !!projects) %>%
    dplyr::group_by(Species = scientific_name, SiteID, Presence) %>%
    dplyr::summarise(Survey = "Camera", Count = Presence) %>%
    dplyr::collect()

  presence_absence <- pa_empty %>%
    dplyr::left_join(presence_absence_ic, by = c("Species", 'SiteID', 'Survey')) %>%
    dplyr::mutate(Presence = case_when(is.na(Presence) ~ 0,
                                       TRUE ~ Presence),
                  Count = case_when(is.na(Count) ~ 0,
                                       TRUE ~ Count))

  transects <- Deer_Detection %>%
    dplyr::bind_rows(presence_absence) %>%
    dplyr::arrange(Species, SiteID, Survey) %>%
    dplyr::mutate(surveyed = 1,
           row = 1:n()) %>%
    as.data.frame()

  transects_grouped <- transects %>%
    dplyr::group_by(SiteID, Transect, Distance, Survey, surveyed) %>%
    dplyr::summarise(Count = max(Count, na.rm = T),
                     Presence = max(Presence, na.rm = T)) %>%
    dplyr::ungroup()

  transects_ne <- site_vars %>%
    dplyr::left_join(transects_grouped) %>%
    dplyr::mutate(surveyed = dplyr::coalesce(surveyed, 0),
                  site = as.numeric(factor(SiteID)))

  transects_ne_all <- site_vars %>%
    dplyr::left_join(transects) %>%
    dplyr::mutate(surveyed = dplyr::coalesce(surveyed, 0),
                  site = as.numeric(factor(SiteID)))

  transects_ne_f <- transects_ne %>%
    dplyr::filter(surveyed > 0)

  n_survey <- transects_ne %>%
    dplyr::group_by(SiteID) %>%
    dplyr::summarise(n = sum(surveyed)) %>%
    dplyr::pull(n)

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

  # Make sure camera is the intercept
  transect_mm <- model.matrix(object = transect_formula, data = transects_grouped)

  # If include cameras in the transect this bit is duplicated
  cam_max <- apply(n_obs, MARGIN = c(1,3), FUN = function(x) {
    m <- max(x)
    if(m > 0) {
      m <- 1
    }
    return(m)
  })

  any_seen <- matrix(0, nrow = length(species), ncol = n_site)
  for(k in 1:length(species)) {
  for (i in 1:n_site) {
    if (n_survey[i] > 0) {
      tran_sp <- transects_ne_all %>%
        dplyr::filter(Species == species[k])
      any_seen[k,i] <- max(tran_sp$Presence[start_idx[i]:end_idx[i]], na.rm = T)
    }
    any_seen[k,i] <- max(c(any_seen[k,i], cam_max[i, k]))
  }
  }

  # Restrict n_max
  n_max <- cam_max
  n_max[n_max == 0] <- n_max_no_det
  n_max[n_max == 1] <- n_max_det

  y2 <- matrix(transects$Presence, nrow = length(species), byrow = TRUE)

  }

  #### Availability ####
  #### Determine availability prior ####
  trigger.events <- camtrap_records_deer %>%
    # filter(is.na(Behaviour)) %>% ### Remove biased behaviour
    dplyr::mutate(Sample.Label = as.character(SiteID))

  trigger.events$rtime <- activity::gettime(x = trigger.events$Time,
                                            tryFormats = "%Y-%m-%d %H:%M:%OS")

  act_result <- activity::fitact(trigger.events$rtime,
                                 sample="model",
                                 reps=50,
                                 bw = 30)

  avail <- list(creation=data.frame(rate = act_result@act[1],
                                    SE   = act_result@act[2]))

  beta.pars<- get.beta.param(avail$creation$rate,
                             2*avail$creation$SE)

  data = list(N=sum(dcount$size, na.rm = T),
              S = length(species),
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
              m_psi =  ncol(ab_model_matrix),
              X_psi = ab_model_matrix,
              npc = nrow(ab_model_pred_matrix),
              X_pred_psi = ab_model_pred_matrix,
              prop_pred = prop_pred,
              coords = coords,
              coords_pred = coords_pred_scaled,
              reciprocal_phi_scale = 1)
  if(evaltransects) {
  data_trans <- list(trans = nrow(transects_grouped),
                     y2 = y2,
                     start_idx = start_idx,
                     end_idx = end_idx,
                     trans_det_ncb = ncol(transect_mm),
                     trans_pred_matrix = transect_mm,
                     any_seen = any_seen,
                     n_survey = n_survey,
                     n_max = n_max)
  data <- c(data, data_trans)
  }

  return(data)

}
