# Power analysis and sample size selection
library(terra)
library(lubridate)
library(VicmapR)
library(sf)
library(dplyr)
library(mapview)

# Camera trap availability
avail <- rbeta(1, 215, 181)

# Survey effort: snapshot moments
dep_times <- seconds(weeks(c(2,4,6,8,10,12)))/2
survey_effort <- as.numeric(42/360 * dep_times * pi * (12.5/1000)^2)

# Average detection probability
det_p <- 0.297

# Simulate detections at sites
buffalo <- vicmap_query("parkres") %>%
  filter(name_short == "Mount Buffalo NP") %>%
  collect(quiet = TRUE) %>%
  st_transform(3111) %>%
  st_combine()

site_ns <- c(10, 20, 30, 40, 50, 60)

sample_locs <- sapply(site_ns, function(x) st_sample(buffalo, x, type = "hexagonal"))

# mapview(buffalo) + mapview(sample_locs[[3]])

combined_raster <- rast("outputs/rasters/combined_deer_average_density.tif")
combined_sd_raster <- rast("outputs/rasters/combined_deer_sd.tif")

abundance_est <- terra::extract(combined_raster, vect(sample_locs[[6]]), ID = F)

count_est <- abundance_est*det_p*avail*survey_effort[4]
