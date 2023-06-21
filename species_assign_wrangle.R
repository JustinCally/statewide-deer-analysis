library(VicmapR)
library(sf)
library(galah)
library(terra)
library(raster)
galah_config(email = "justin.cally@delwp.vic.gov.au")

con <- weda_connect(password = keyring::key_get(service = "ari-dev-weda-psql-01",
                                                username = "psql_user"))
#### Make KDEs ####
vic_bound <- vicmap_query("open-data-platform:vmlite_victoria_polygon_su5") %>%
  filter(state == "VIC") %>% # Philip Island and French Island
  collect() %>%
  st_transform(3111)

area <- vic_bound %>%
  sf::st_union() %>%
  sf::st_buffer(25000) %>%
  sf::st_simplify(dTolerance = 10000) %>%
  sf::st_transform(4283) %>%
  sf::st_as_sf()

if(!file.exists("data/historic_deer_ala_records.rds")) {
  historic_deer_ala_records <- galah_call() %>%
    galah_identify(c("Cervus unicolor", "Dama dama", "Cervus elaphas")) %>%
    galah::galah_geolocate(area) %>%
    galah_filter(year >= 2003) %>%
    atlas_occurrences() %>%
    filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4283) %>%
    mutate(scientificName = case_when(scientificName == "Dama dama dama" ~ "Dama dama",
                                      TRUE ~ scientificName))
  saveRDS(historic_deer_ala_records, "data/historic_deer_ala_records.rds")
} else {
  historic_deer_ala_records <- readRDS("data/historic_deer_ala_records.rds")
}


smoothed_kernel_species <- function(sr, all_area = vic_bound) {

  hist_3111 <- sr %>%
    sf::st_transform(3111)

  points_kde <- st_kde(hist_3111,
                       cellsize = 1000,
                       bandwith = nrd_new(hist_3111 %>% sf::st_coordinates()),
                       extent = terra::ext(all_area))

  pred_vals <- raster::extract(points_kde, hist_3111)
  q995 <- quantile(pred_vals, 0.001, na.rm = T)

  points_kde_cut <- points_kde
  terra::values(points_kde_cut)[terra::values(points_kde_cut) >= q995] <- as.factor("1")
  terra::values(points_kde_cut)[terra::values(points_kde_cut) < q995] <- NA
  points_kde_cut <- terra::rast(points_kde_cut)
  points_kde_cut <- terra::crop(points_kde_cut, terra::vect(vic_bound), mask = TRUE)

  points_kde_cut
}

split_deer_records <- split(historic_deer_ala_records, historic_deer_ala_records$scientificName)

deer_kde <- lapply(split_deer_records, smoothed_kernel_species) %>%
  terra::rast()

names(deer_kde) <- paste(names(deer_kde), "distribution")

#### For the transect data assign it to the three species ####
deer_species <- c("Sambar Deer", "Fallow Deer", "Red Deer")
project_short_name <- c("hog_deer_2023", "StatewideDeer")

raw_transects <- tbl(con, dbplyr::in_schema(schema = "deervic", table = "raw_site_transects")) %>%
  collect() %>%
  na.omit()

presence_absence <-  tbl(con, in_schema("camtrap", "curated_camtrap_operation")) %>%
  filter(ProjectShortName %in% !!project_short_name) %>%
  left_join(tbl(con, in_schema("camtrap", "processed_site_substation_presence_absence")) %>%
              filter(common_name %in% deer_species) %>%
              dplyr::select(-SubStation)) %>%
  collect() %>%
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4283)

presence_absence_wide <- presence_absence %>%
  pivot_wider(id_cols = "SiteID", values_from = "Presence", names_from ="scientific_name") %>%
  st_drop_geometry() %>%
  left_join(presence_absence %>% dplyr::select(SiteID) %>% distinct()) %>%
  st_as_sf()

species_overlap <- bind_cols(presence_absence_wide,
                             terra::extract(deer_kde, vect(presence_absence_wide %>%
                                                             sf::st_transform(3111)))) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(-ID)

species_overlap[is.na(species_overlap)] <- 0

site_vars <- tbl(con, in_schema("deervic", "curated_site_data")) %>%
  filter(SiteID %in% !!c(presence_absence$SiteID, "47191")) %>%
  collect() %>%
  mutate(HerbaceousUnderstoryCover = NNWHUCover + ENWHUCover,
         SiteID = case_when(SiteID == "47191" & CameraID == "HO04101053" ~ "47191A",
                            SiteID == "47191" & CameraID != "HO04101053" ~ "47191B",
                            TRUE ~ SiteID)) %>% # native + exotic herbaceous cover
  arrange(SiteID) %>%
  as.data.frame()

transects_long <- tbl(con, dbplyr::in_schema(schema = "deervic", table = "raw_site_transects")) %>%
  filter(Data_File_Id %in% !!site_vars$Data_File_Id) %>%
  collect() %>%
  dplyr::select(Data_File_Id,
         Transect = I0_TransectNo,
         Distance = I0_TransDistance,
         Pellet = I0_Pellets,
         Rubbing = I0_Rubbings,
         Footprint = I0_Footprints,
         Wallow = I0_Wallows) %>%
  mutate(Footprint = case_when(Footprint == "Observed" ~ 1L,
                                TRUE ~ 0L)) %>%
  left_join(site_vars %>%
              dplyr::select(Data_File_Id, SiteID) %>%
              st_drop_geometry(), join_by(Data_File_Id)) %>%
  dplyr::select(SiteID, everything(), -Data_File_Id) %>%
  # pivot_longer(cols = Pellets:Wallows, names_to = c("Survey"), values_to = "Presence") %>%
  na.omit()

transects_sp <- left_join(transects_long, species_overlap, by = "SiteID")

Damadama_Detection <- detection_assign(transects_sp, "Dama dama", othersp = c("Cervus unicolor", "Cervus elaphus"))
Rusaunicolor_Detection <- detection_assign(transects_sp, "Cervus unicolor", othersp = c("Dama dama", "Cervus elaphus"))
Cervuselaphus_Detection <- detection_assign(transects_sp, "Cervus elaphus", othersp = c("Dama dama", "Cervus unicolor"))

saveRDS(Damadama_Detection, "data/Damadama_Detection.rds")
saveRDS(Rusaunicolor_Detection, "data/Rusaunicolor_Detection.rds")
saveRDS(Cervuselaphus_Detection, "data/Cervuselaphus_Detection.rds")
