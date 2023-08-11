con <- weda::weda_connect(password = keyring::key_get(service = "ari-dev-weda-psql-01",
                                                      username = "psql_user"), username = "psql_user")

cams_curated_buffer <- tbl(con, dbplyr::in_schema("camtrap", "curated_camtrap_operation")) %>%
  dplyr::filter(ProjectShortName %in% !!project_short_name) %>% # Only retrieve cam records for hog deer 2023
  dplyr::collect() %>%
  dplyr::arrange(SiteID) %>%
  sf::st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4283)

evc_data <- sf::st_read("data/nv2005_evcbcs")

evc_data_valid <- st_make_valid(evc_data) %>%
  filter()

evc_joined <- cams_curated_buffer %>%
  st_join(evc_data_valid)

# bioreg - evc codes
evc_joined_mutate <- evc_joined %>%
  mutate(BIOEVC = paste(XGROUPNAME, sep = " || ")) %>%
  arrange(BIOEVC) %>%
  mutate(BIOEVCLVL = as.integer(as.factor(BIOEVC)))

# evc data all
evc_data_all <- evc_data_valid %>%
  mutate(BIOEVC = paste(XGROUPNAME, sep = " || "))

# filter to sampled areas
evc_data_filtered <- evc_data_all %>%
  filter(BIOEVC %in% evc_joined_mutate$BIOEVC) %>%
  arrange(BIOEVC) %>%
  mutate(BIOEVCLVL = as.integer(as.factor(BIOEVC)))

pred_raster_full <- terra::rast(prediction_raster)

pred_raster <- terra::app(pred_raster_full[[stringr::str_subset(
  stringr::str_remove_all(labels(terms(ab_formula_4)),
                          "scale[(]|[)]|log[(]|sqrt[(]"),
  pattern = "[*]", negate = T)]], mean)

vic_model_data_resampled_df <- terra::as.data.frame(pred_raster, xy = TRUE, cell = TRUE, na.rm = TRUE)

pred_evc_lvl_i <-  vic_model_data_resampled_df %>%
  sf::st_as_sf(., coords = c("x", "y"), crs = 3111) %>%
  sf::st_nearest_feature(evc_data_filtered %>% st_transform(3111))

evc_groups <- list()
evc_groups[["pred_evc"]] <- evc_data_filtered$BIOEVCLVL[pred_evc_lvl_i]
evc_groups[["site_evc"]] <- evc_joined_mutate$BIOEVCLVL
evc_groups[["np_evc"]] <- length(unique(evc_groups[["site_evc"]]))
saveRDS(evc_groups, "data/evc_groupnames.rds")

# for(i in 1:length(model_data)) {
#   model_data[[i]]$np_evc <- evc_groups[["np_evc"]]
#   model_data[[i]]$pred_evc <- evc_groups[["pred_evc"]]
#   model_data[[i]]$site_evc <- evc_groups[["site_evc"]]
# }
