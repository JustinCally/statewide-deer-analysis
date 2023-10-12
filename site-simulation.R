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

survey_area <- st_area(buffalo) %>% units::set_units("km2")

site_ns <- c(10, 20, 30, 40, 50, 60, 80, 100, 120, 200, 300)

site_idx <- 11

sample_locs <- sapply(site_ns, function(x) st_sample(buffalo, x, type = "hexagonal"))

# mapview(buffalo) + mapview(sample_locs[[3]])

combined_raster <- rast("outputs/rasters/combined_deer_average_density.tif")
combined_sd_raster <- rast("outputs/rasters/combined_deer_sd.tif")

abundance_area <- terra::extract(combined_raster, vect(buffalo), ID = F)

abundance_est <- terra::extract(combined_raster, vect(sample_locs[[site_idx]]), ID = F)
sd_est <- terra::extract(combined_sd_raster, vect(sample_locs[[site_idx]]), ID = F)

species_idx <- 1

raw_counts <- matrix(nrow = length(sample_locs[[site_idx]]), ncol = 1000)
dens_est <- matrix(nrow = length(sample_locs[[site_idx]]), ncol = 1000)

for(i in 1:length(sample_locs[[site_idx]])) {
raw_counts[i,] <- rnbinom(1000,mu=abundance_est[i,species_idx]*det_p*avail*survey_effort[4],size=abundance_est[i,species_idx]^2/((sd_est[i,species_idx]^2)-abundance_est[i,species_idx]))
dens_est[i,] <- raw_counts[i,] /(det_p*avail*survey_effort[4])
}

dens_mean <- vector()
ab_mean <- vector()
for(j in 1:1000) {
  dens_mean[j] <- mean(dens_est[,j])
  ab_mean[j] <- dens_mean[j]*survey_area
}

grand_mean <- mean(ab_mean)
grand_cv <- sd(ab_mean)/grand_mean
grand_mean
grand_cv

#true mean
sum(abundance_area[,species_idx])
