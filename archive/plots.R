library(ggplot2)
library(ggtext)
library(VicmapR)
library(sf)
library(dplyr)
project_short_name <- c("hog_deer_2023", "StatewideDeer")
#### Database connection ####
con <- weda::weda_connect(password = keyring::key_get(service = "ari-dev-weda-psql-01",
                                                      username = "psql_user"), username = "psql_user")

cams_curated <- tbl(con, dbplyr::in_schema("camtrap", "curated_camtrap_operation")) %>%
  dplyr::filter(ProjectShortName %in% !!project_short_name) %>% # Only retrieve cam records for hog deer 2023
  dplyr::collect() %>%
  dplyr::arrange(SiteID)

# fitintegrn <- readRDS("outputs/fitintegrn.rds")

rn_dens <- fitintegrn_hog$summary("N_site")
# Bind the density to the camera information
# Regions
vic_regions <- vicmap_query("open-data-platform:delwp_region") %>%
  collect() %>%
  st_transform(3111) %>%
  st_simplify(dTolerance = 500)

density_at_sites_rn <- cbind(cams_curated, rn_dens) %>%
  mutate(means_sqrt = sqrt(mean)) %>%
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4283) %>%
  st_transform(3111)



density_at_sites_rn$density <- cut(density_at_sites_rn$mean,
                                   breaks = c(0, 1, 5, 10, max(density_at_sites_rn$mean)),
                                   labels = c("<1", "1-5", "5-10", "10+"), include.lowest = T, right = T)

ggplot(data = density_at_sites_rn) +
  geom_sf(data = vic_regions, alpha = 0.75, fill = "grey80") +
  geom_sf(aes(fill = density, alpha = mean), shape = 21, size = 4) +
  scale_fill_viridis_d(name = "", guide = guide_legend(override.aes = list(size = 6))) +
  scale_alpha_continuous(range = c(0.5,1), guide = "none") +
  labs(title = bquote('Density of Hog Deer per'~km^2)) +
  theme_bw() +
  theme(legend.text = element_text(size = 18), legend.key.size = unit(1, "cm"),
        title = element_text(size = 22))
