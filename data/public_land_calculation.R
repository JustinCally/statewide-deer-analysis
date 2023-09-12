# creating the public land area to use
sf_use_s2(FALSE)
plm <- vicmap_query("open-data-platform:plm25") %>%
  filter(!(gener_desc %in% c("INLAND WATER BODY", "MARINE")),
         rec_cat != "SERVICES AND UTILITIES") %>%
  collect()

plm_f <- plm %>%
  st_transform(3111) %>%
  st_make_valid()

wa <- vicmap_query("open-data-platform:hy_water_area_polygon") %>%
  filter(feature_type_code %in% c("wb_lake", "wb_lake_salt", "reservoir", "pondage_sewerage", "pondage", "pondage_saltpan", "watercourse_area_drain", "watercourse_area_river", "watercourse_area_channel", "watercourse_area_channel_drain")) %>%
  collect()

wa_v <- wa %>%
  st_transform(3111) %>%
  st_filter(plm_f)

plm_to_filter <- st_filter(plm_f, wa_v)

split_vector <- rep(1:400, each = nrow(plm_to_filter) / nrow(plm_to_filter), length.out = nrow(plm_to_filter))

split_results <- split(plm_to_filter, split_vector)

plm_nw_t <- list()
t1 <- Sys.time()
for(i in 1:length(split_results)) {
  cat(i)
  sp1 <- wa_v %>%
    st_filter(split_results[[i]]) %>%
    st_make_valid()

  gc(verbose = F)
  try({
    plm_nw_t[[i]] <- rmapshaper::ms_erase(split_results[[i]], sp1)
  })
}
Sys.time()-t1

plm_nw <- bind_rows(plm_nw_t)

plm_old_f <- plm_f %>%
  filter(!(id %in% plm_nw$id))

plm_c <- bind_rows(plm_nw, plm_old_f)

crown_land_by_reserve <- plm_c %>%
  group_by(label) %>%
  summarise(geometry = st_union(geometry),
            area = sum(area_ha))

saveRDS(crown_land_by_reserve, "data/crown_land_by_reserve.rds", compress = "xz")
saveRDS(plm_c, "data/plm_c.rds", compress = "xz")

crown_land <- plm_c %>% st_union()

# filter islands

saveRDS(crown_land, "data/crown_land.rds", compress = "xz")

# Read in resolved rasters
processed_tifs <- list.files(raster_files, full.names = TRUE)[stringr::str_detect(list.files(raster_files), ".tif$")]

processed_stack <- terra::rast(processed_tifs)

# project to 1km resolution
projected_sw <- project(processed_stack, rast(prediction_raster))

# crop to public land
cropped_prediction <- mask(projected_sw,
                           vect(crown_land),
                           touches=TRUE)

writeRaster(cropped_prediction, "data/prediction_raster/statewide_raster.tif", overwrite=T)
